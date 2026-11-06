//! tcrdist-style cluster discovery: joint alpha/beta distance graph, clustered per threshold.

use indicatif::ProgressIterator;
mod clustering;
mod distance;
mod lookup_cdr_seqs;

pub use lookup_cdr_seqs::TcrDistDb;
use lookup_cdr_seqs::{ChainCdrs, VGeneCdrs, infer_record_cdrs};

use crate::progress::{Progress, ProgressSettings, dual_progress_bars};
use crate::tcr_cluster::blosum::AminoAcidSubMatrix;
use crate::tcr_cluster::model::*;
use crate::tcr_cluster::pipeline::ClusterFinder;
use crate::tcr_cluster::stats;
use anyhow::{Context, Result, bail};
pub use clustering::PairedCluster;
use clustering::build_clusters_for_threshold;
use distance::TcrDistanceMatrix;
use rayon::prelude::*;
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::sync::Mutex;

const THRESHOLDS: [u8; 6] = [4, 8, 12, 18, 24, 30];
/// Clusters containing at least this fraction of all records are considered uninformative.
const MAX_CLUSTER_FRACTION: f32 = 0.25;

/// The tcrdist cluster finder.
pub struct TcrDist {
    pub sub_matrix: AminoAcidSubMatrix,
    pub params: DistanceParams,
    /// Ignore `params.mode` and pick it from the input with [`DistanceMode::infer`].
    pub infer_mode: bool,
    /// tcrdist3 gene db; required by the modes that use germline CDRs.
    pub gene_db: Option<TcrDistDb>,
}

impl TcrDist {
    /// The distance parameters for `records`, resolving an inferred mode.
    fn params_for(&self, records: &[ChainRecord]) -> DistanceParams {
        let mode = if self.infer_mode {
            DistanceMode::infer(records)
        } else {
            self.params.mode
        };
        let cdr3_weight = if mode.uses_germline_cdrs() {
            self.params.cdr3_weight
        } else {
            1
        };
        DistanceParams {
            mode,
            cdr3_weight,
            ..self.params
        }
    }

    fn gene_db_for(&self, mode: DistanceMode) -> Result<Option<&TcrDistDb>> {
        if !mode.uses_germline_cdrs() {
            return Ok(None);
        }
        self.gene_db
            .as_ref()
            .map(Some)
            .with_context(|| format!("Distance mode '{mode}' needs the tcrdist3 gene database"))
    }
}

/// Scales `matrix` in place; weight 1 is a no-op.
fn weighted(mut matrix: TcrDistanceMatrix, weight: u8) -> TcrDistanceMatrix {
    if weight != 1 {
        matrix *= weight;
    }
    matrix
}

/// Joint tcrdist distance between all `records` for the chains and CDRs selected by
/// `params.mode`. Every record must carry the data the mode needs (see `TcrDist::prepare`).
fn joint_distance_matrix(
    records: &[ChainRecord],
    params: &DistanceParams,
    gene_db: Option<&TcrDistDb>,
    sub_matrix: &AminoAcidSubMatrix,
) -> Result<TcrDistanceMatrix> {
    let mode = params.mode;
    let mismatch_matrix = sub_matrix.build_mismatch_matrix(i32::from(params.max_mismatch));
    let cdr3 = |seqs: &[String], name: &str| {
        distance::compute_matrix_auto_detect(
            seqs,
            &mismatch_matrix,
            i32::from(params.gap_penalty),
            u8::MAX,
            name,
        )
        .with_context(|| format!("Could not calculate {name} distance matrix"))
        .map(|matrix| weighted(matrix, params.cdr3_weight))
    };

    let beta: Vec<String> = records.iter().map(|r| r.beta.clone()).collect();
    let mut total = cdr3(&beta, "CDR3β")?;
    if mode.uses_alpha() {
        let alpha: Vec<String> = records
            .iter()
            .map(|r| r.alpha.clone().unwrap_or_default())
            .collect();
        total += cdr3(&alpha, "CDR3α")?;
    }

    if let Some(db) = gene_db.filter(|_| mode.uses_germline_cdrs()) {
        let cdrs: Vec<ChainCdrs> = records
            .iter()
            .map(|record| infer_record_cdrs(db, record, true))
            .collect();
        // Germline distances are added straight into `total`: no extra matrix is built.
        let mut add_germline = |chain: fn(&ChainCdrs) -> Option<&VGeneCdrs>, name: &str| {
            let chain_cdrs: Vec<&VGeneCdrs> = cdrs
                .iter()
                .map(|c| chain(c).with_context(|| format!("Record without {name}")))
                .collect::<Result<_>>()?;
            distance::add_germline_distances(
                &mut total,
                &chain_cdrs,
                &mismatch_matrix,
                params.gap_penalty,
                params.germline_weight,
                name,
            )
        };
        add_germline(|c| c.beta.as_ref(), "germline CDRβ")?;
        if mode.uses_alpha() {
            add_germline(|c| c.alpha.as_ref(), "germline CDRα")?;
        }
    } else if mode.uses_germline_cdrs() {
        bail!("Distance mode '{mode}' needs the tcrdist3 gene database");
    }
    Ok(total)
}

impl ClusterFinder for TcrDist {
    type Cluster = PairedCluster;

    /// Drops records missing what the distance mode needs, so that matrix indices,
    /// clustering and scoring all refer to the same records.
    fn prepare<'a>(&self, records: &'a [ChainRecord]) -> Result<Cow<'a, [ChainRecord]>> {
        let params = self.params_for(records);
        let mode = params.mode;
        let db = self.gene_db_for(mode)?;
        let (mut missing_alpha, mut unknown_v) = (0usize, 0usize);
        let kept: Vec<ChainRecord> = records
            .iter()
            .filter(|record| {
                if mode.uses_alpha() && !record.has_alpha_cdr3() {
                    missing_alpha += 1;
                    return false;
                }
                if let Some(db) = db {
                    let cdrs = infer_record_cdrs(db, record, true);
                    if cdrs.beta.is_none() || (mode.uses_alpha() && cdrs.alpha.is_none()) {
                        unknown_v += 1;
                        return false;
                    }
                }
                true
            })
            .cloned()
            .collect();

        eprintln!("tcrdist distance mode: {mode}");
        if missing_alpha > 0 {
            eprintln!("Warning: {missing_alpha} records without a CDR3α will be ignored");
        }
        if unknown_v > 0 {
            eprintln!(
                "Warning: {unknown_v} records with a missing or unknown V gene for {} will be ignored",
                db.map_or_else(String::new, |db| db.organism().to_string())
            );
        }
        if kept.is_empty() {
            bail!("No records have the data needed for distance mode '{mode}'");
        }
        Ok(if kept.len() == records.len() {
            Cow::Borrowed(records)
        } else {
            Cow::Owned(kept)
        })
    }

    fn find(
        &self,
        records: &[ChainRecord],
        reference: &ReferenceData,
    ) -> Result<Vec<PairedCluster>> {
        let ref_seqs = reference.valid_sequences();
        let params = self.params_for(records);
        params.check_thresholds(&THRESHOLDS);
        let distance_matrix = joint_distance_matrix(
            records,
            &params,
            self.gene_db_for(params.mode)?,
            &self.sub_matrix,
        )?;
        let beta_seqs: Vec<String> = records.iter().map(|r| r.beta.clone()).collect();

        let fisher_lookup: HashMap<String, f64> =
            stats::motif_statistics_iter(&beta_seqs, &ref_seqs)
                .map(|(motif, stats)| (motif, stats.fisher))
                .collect();
        let log_fact_table = stats::LogFactorialTable::new(ref_seqs.len() + beta_seqs.len());
        let fisher_cache: Mutex<HashMap<String, f64>> = Mutex::new(HashMap::new());
        let fisher_for = |pattern: &str| -> f64 {
            if let Some(&v) = fisher_lookup.get(pattern) {
                return v;
            }
            if let Some(&v) = fisher_cache.lock().unwrap().get(pattern) {
                return v;
            }
            let v = stats::pattern_fisher(pattern, &beta_seqs, &ref_seqs, &log_fact_table);
            fisher_cache.lock().unwrap().insert(pattern.to_string(), v);
            v
        };

        let max_members = (records.len() as f32 * MAX_CLUSTER_FRACTION) as usize;

        let (opb, ipb) = dual_progress_bars(
            THRESHOLDS.len() as u64,
            None,
            ProgressSettings::Position,
            ProgressSettings::Position,
            Progress::Progress,
        )?;
        let discovered: HashSet<PairedCluster> = THRESHOLDS
            .into_iter()
            .map(|threshold| {
                opb.set_message(format!("Clustering with max distance {threshold}"));
                // Show that the inner bar is idle while the graph is built, instead of
                // leaving the previous threshold's count on screen.
                ipb.reset();
                ipb.set_length(0);
                ipb.set_message("Building graph");
                if distance_matrix.max() >= threshold {
                    eprintln!("Skipping threshold {threshold} because the resulting graph would contain 1 cluster");
                    return Ok(Vec::new());
                }
                let graph = distance_matrix
                    .to_igraph(params.scaled_threshold(threshold))
                    .context("Failed to convert distance matrix to graph")?;
                build_clusters_for_threshold(
                    threshold,
                    &graph,
                    records,
                    &self.sub_matrix,
                    true,
                    &fisher_for,
                    &ClusterKind::Global,
                    max_members,
                    &ipb,
                )
                .context("Failed to build clusters for threshold")
            })
            .progress_with(opb.clone())
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .flatten()
            .collect();

        ipb.finish_and_clear();
        opb.finish();

        // HashSet order changes per process; the pipeline numbers clusters (and breaks output
        // ties) by position, so return them in a fixed order.
        let mut clusters: Vec<PairedCluster> = discovered.into_iter().collect();
        clusters.sort_unstable_by(|a, b| {
            (
                a.threshold,
                a.cluster_id,
                &a.member_indices,
                &a.cluster.pattern,
            )
                .cmp(&(
                    b.threshold,
                    b.cluster_id,
                    &b.member_indices,
                    &b.cluster.pattern,
                ))
        });
        Ok(clusters)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const DB: &str =
        "id\torganism\tchain\tregion\tnucseq\tframe\taligned_protseq\tcdr_columns\tcdrs
TRBV19*01\thuman\tB\tV\tg\t+1\tX\t-\tMN.HEY;SVG;TDQ;CASS
TRBV9*01\thuman\tB\tV\tg\t+1\tX\t-\tSG.HVS;FQN;SEQ;CASS
TRAV12-2*01\thuman\tA\tV\tg\t+1\tX\t-\tDRG.SQS;IYS;NKG;CAVN
TRAV1-2*01\thuman\tA\tV\tg\t+1\tX\t-\tTSG.FNG;VVL;DGL;CAVR
";

    fn record(beta: &str, alpha: Option<&str>, vb: Option<&str>, va: Option<&str>) -> ChainRecord {
        ChainRecord {
            beta: beta.to_owned(),
            alpha: alpha.map(str::to_owned),
            v_gene: vb.map(str::to_owned),
            va_gene: va.map(str::to_owned),
            ..ChainRecord::default()
        }
    }

    fn finder(mode: DistanceMode) -> TcrDist {
        TcrDist {
            sub_matrix: AminoAcidSubMatrix::Blosum62,
            params: DistanceParams::new(mode),
            infer_mode: false,
            gene_db: Some(TcrDistDb::from_tsv(DB, Organism::Human).unwrap()),
        }
    }

    fn distances(mode: DistanceMode, records: &[ChainRecord]) -> TcrDistanceMatrix {
        let f = finder(mode);
        let params = f.params_for(records);
        joint_distance_matrix(
            records,
            &params,
            f.gene_db_for(mode).unwrap(),
            &f.sub_matrix,
        )
        .unwrap()
    }

    #[test]
    fn each_mode_adds_its_chains_and_cdrs() {
        // Same CDR3β; differ in CDR3α, TRBV and TRAV.
        let records = [
            record(
                "CASSIRSSYEQYF",
                Some("CAVRDSNYQLIW"),
                Some("TRBV19"),
                Some("TRAV12-2"),
            ),
            record(
                "CASSIRSSYEQYF",
                Some("CAVRDSNYQLVW"),
                Some("TRBV9"),
                Some("TRAV1-2"),
            ),
        ];
        let beta = distances(DistanceMode::Beta, &records)[(1, 0)];
        let alpha_beta = distances(DistanceMode::AlphaBeta, &records)[(1, 0)];
        let beta_all = distances(DistanceMode::BetaFull, &records)[(1, 0)];
        let all = distances(DistanceMode::AlphaBetaFull, &records)[(1, 0)];

        assert_eq!(beta, 0);
        assert!(alpha_beta > beta, "CDR3α mismatch counts");
        assert!(beta_all > beta, "TRBV germline CDRs count");
        // all = CDR3α + CDR3β + germline β + germline α (TRAVs differ, so germline α > 0).
        assert!(
            all > alpha_beta + (beta_all - beta),
            "germline α contributes on top"
        );
    }

    #[test]
    fn germline_distance_is_zero_for_same_gene_with_and_without_allele() {
        let records = [
            record("CASSIRSSYEQYF", None, Some("TRBV19"), None),
            record("CASSIRSSYEQYF", None, Some("TRBV19*01"), None),
        ];
        assert_eq!(distances(DistanceMode::BetaFull, &records)[(1, 0)], 0);
    }

    #[test]
    fn prepare_drops_records_missing_mode_data() {
        let records = vec![
            record(
                "CASSIRSSYEQYF",
                Some("CAVRDSNYQLIW"),
                Some("TRBV19"),
                Some("TRAV12-2"),
            ),
            record("CASSPRSSYEQYF", None, Some("TRBV19"), Some("TRAV12-2")),
            record(
                "CASSLRSSYEQYF",
                Some("CAVRDSNYQLIW"),
                Some("TRBV99"),
                Some("TRAV12-2"),
            ),
            record("CASSKRSSYEQYF", Some("CAVRDSNYQLIW"), Some("TRBV9"), None),
        ];
        let kept = |mode| finder(mode).prepare(&records).unwrap().len();

        assert_eq!(kept(DistanceMode::Beta), 4);
        assert_eq!(kept(DistanceMode::AlphaBeta), 3);
        assert_eq!(kept(DistanceMode::BetaFull), 3, "unknown TRBV99 dropped");
        assert_eq!(kept(DistanceMode::AlphaBetaFull), 1);
    }

    #[test]
    fn inferred_mode_uses_alpha_only_when_present() {
        let beta_only = [record("CASSIRSSYEQYF", None, None, None)];
        let paired = [record("CASSIRSSYEQYF", Some("CAVRDSNYQLIW"), None, None)];
        let f = TcrDist {
            infer_mode: true,
            gene_db: None,
            ..finder(DistanceMode::AlphaBetaFull)
        };
        assert_eq!(f.params_for(&beta_only).mode, DistanceMode::Beta);
        assert_eq!(f.params_for(&paired).mode, DistanceMode::AlphaBeta);
    }

    #[test]
    fn germline_mode_without_db_is_an_error() {
        let f = TcrDist {
            gene_db: None,
            ..finder(DistanceMode::BetaFull)
        };
        assert!(
            f.prepare(&[record("CASSIRSSYEQYF", None, Some("TRBV19"), None)])
                .is_err()
        );
    }
}
