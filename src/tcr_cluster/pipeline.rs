//! The shared find → score → format pipeline used by every clustering algorithm.

use crate::progress::{ProgressSettings, progress_bar_with_message};
use crate::tcr_cluster::model::*;
use crate::tcr_cluster::scoring::score_cluster;
use anyhow::{Result, bail};
use rayon::prelude::*;
use std::borrow::Cow;

/// A strategy that groups chains into clusters (GLIPH, tcrdist, ...).
pub trait ClusterFinder: Sync {
    type Cluster: HasCluster + Send + Sync;

    /// Optionally restricts the records the finder (and scoring) operate on.
    fn prepare<'a>(&self, records: &'a [ChainRecord]) -> Result<Cow<'a, [ChainRecord]>> {
        Ok(Cow::Borrowed(records))
    }

    fn find(
        &self,
        records: &[ChainRecord],
        reference: &ReferenceData,
    ) -> Result<Vec<Self::Cluster>>;
}

/// Parameters shared by all algorithms.
#[derive(Clone, Copy, Debug)]
pub struct ScoringParams {
    pub algorithm: Algorithm,
    pub sim_depth: usize,
    pub hla_threshold: f64,
}

/// Scores every cluster in parallel and builds the output rows for its members.
/// Clusters without any matching record are dropped.
pub fn score_clusters<C>(
    records: &[ChainRecord],
    clusters: Vec<C>,
    ctx: &ScoringContext,
) -> Result<Vec<ScoredCluster<C>>>
where
    C: HasCluster + Send + Sync,
{
    let progress = progress_bar_with_message(
        clusters.len(),
        "Scoring clusters",
        ProgressSettings::PositionSpeed,
    )?;
    let scored: Vec<ScoredCluster<C>> = clusters
        .into_par_iter()
        .enumerate()
        .filter_map(|(index, cluster)| {
            let result = {
                let members = ClusterMembers::new(records, cluster.cluster());
                if members.is_empty() {
                    None
                } else {
                    let score = score_cluster(ctx, &members);
                    let rows = members
                        .records
                        .iter()
                        .map(|record| OutputRow::new(index, &members, &score, record))
                        .collect();
                    Some((score, rows))
                }
            };
            progress.inc(1);
            result.map(|(score, rows)| ScoredCluster {
                cluster,
                score,
                rows,
            })
        })
        .collect();
    progress.finish_and_clear();
    eprintln!(
        "Cluster scoring complete in {}",
        indicatif::HumanDuration(progress.elapsed())
    );
    Ok(scored)
}

/// Runs `finder` on `input_records`, scores the clusters and returns the sorted output rows.
pub fn run<F: ClusterFinder>(
    finder: &F,
    input_records: &[ChainRecord],
    reference: &ReferenceData,
    params: ScoringParams,
) -> Result<Vec<OutputRow>> {
    if params.sim_depth == 0 {
        bail!("Simulation depth must be greater than zero");
    }
    let hla = HlaContext::build(input_records, params.hla_threshold);
    let expansion = ExpansionContext::build(input_records);

    let records = finder.prepare(input_records)?;
    let clusters = finder.find(&records, reference)?;
    let ctx = ScoringContext::new(
        reference,
        params.algorithm,
        params.sim_depth,
        &hla,
        &expansion,
        records.len(),
    );

    let scored = score_clusters(&records, clusters, &ctx)?;
    let mut rows: Vec<OutputRow> = scored.into_iter().flat_map(|c| c.rows).collect();
    rows.sort_by(|left, right| {
        left.fisher_score
            .total_cmp(&right.fisher_score)
            .then_with(|| left.pattern.cmp(&right.pattern))
            .then_with(|| left.index.cmp(&right.index))
    });
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tcr_cluster::blosum::AminoAcidSubMatrix;
    use crate::tcr_cluster::gliph::Gliph;
    use pretty_assertions::assert_eq;

    fn map_test_sample_to_chain_record(index: usize, beta: &str) -> ChainRecord {
        ChainRecord {
            subject: format!("S{}", index + 1),
            condition: "case".to_owned(),
            beta: (*beta).to_owned(),
            alpha: Some(String::new()),
            v_gene: Some("TRBV1".to_owned()),
            j_gene: Some("TRBJ1".to_owned()),
            va_gene: None,
            ja_gene: None,
            counts: None,
            hla: [
                "DPA1*01:03".to_owned(),
                "DPB1*04:01".to_owned(),
                "DQA1*05:01".to_owned(),
                "DQB1*02:01".to_owned(),
                "DRB1*15:01".to_owned(),
                String::new(),
                String::new(),
                String::new(),
            ],
        }
    }
    fn run_gliph(
        records: &[ChainRecord],
        reference: &ReferenceData,
        algorithm: Algorithm,
    ) -> Vec<OutputRow> {
        let finder = Gliph {
            algorithm,
            sim_depth: 20,
            sub_matrix: Some(AminoAcidSubMatrix::Blosum62),
        };
        let params = ScoringParams {
            algorithm,
            sim_depth: 20,
            hla_threshold: 0.05,
        };
        run(&finder, records, reference, params).unwrap()
    }

    #[test]
    fn gliph_pipeline_builds_annotated_cluster_rows() {
        let sample = ["CASSIRSSYEQYF", "CASSIRSSYEQAF", "CASSIRSSYEQVF"];
        let records = sample
            .iter()
            .enumerate()
            .map(|(index, beta)| map_test_sample_to_chain_record(index, beta))
            .collect::<Vec<_>>();
        let alphabet = b"ACDEFGHIKLMNPQTVWY";
        let reference_sequences = (0..30)
            .map(|index| {
                let first = alphabet[index % alphabet.len()] as char;
                let second = alphabet[(index / alphabet.len()) % alphabet.len()] as char;
                format!("CAAAA{first}{second}F")
            })
            .collect::<Vec<_>>();
        let reference = ReferenceData {
            sequences: reference_sequences,
            v_gene_frequencies: vec![("TRBV1".to_owned(), 1.0)],
            cdr3_length_frequencies: vec![(13, 1.0)],
            cluster_size_probabilities: Vec::new(),
            cluster_size_sample_sizes: Vec::new(),
        };

        let result = run_gliph(&records, &reference, Algorithm::Gliph2);

        assert!(!result.is_empty());
        assert!(result.iter().all(|row| row.number_subject == 3));
        assert!(result.iter().all(|row| row.number_unique_cdr3 >= 2));
        assert!(result.iter().all(|row| !row.hla_dpa1.is_empty()));
        assert_eq!(result, run_gliph(&records, &reference, Algorithm::Gliph2));

        let gliph1_result = run_gliph(&records, &reference, Algorithm::Gliph1);
        assert!(!gliph1_result.is_empty());
        assert!(gliph1_result.iter().all(|row| row.kind == "global"));
    }
}
