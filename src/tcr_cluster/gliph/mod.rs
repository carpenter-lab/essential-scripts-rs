//! GLIPH1 / GLIPH2 cluster discovery (motif enrichment, global similarity, merging).

mod gliph1;
mod global;
mod local;
mod merge;
mod singletons;

use crate::progress::progress_bar_with_message;
use crate::tcr_cluster::blosum::AminoAcidSubMatrix;
use crate::tcr_cluster::model::*;
use crate::tcr_cluster::pipeline::ClusterFinder;
use crate::tcr_cluster::stats::{ContingencyTable2x2, LogFactorialTable};
use crate::tcr_cluster::{blosum, stats};
use anyhow::{Result, bail};
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;
use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use self::{gliph1::*, global::*, local::*, merge::*, singletons::*};

const MOTIF_MINIMUM_FOLD_CHANGE: [f64; 3] = [1000.0, 100.0, 10.0];
const MIN_SHARED_MEMBERS: usize = 4;
const MIN_OVERLAP_FRACTION: f64 = 0.80;

/// The GLIPH1 / GLIPH2 cluster finder.
pub struct Gliph {
    pub algorithm: Algorithm,
    pub sim_depth: usize,
    pub sub_matrix: Option<AminoAcidSubMatrix>,
}

/// A GLIPH cluster with the records it covers: every record whose CDR3β is a member.
pub struct GliphCluster {
    pub cluster: Cluster,
    pub member_indices: Vec<usize>,
}

impl HasCluster for GliphCluster {
    fn cluster(&self) -> &Cluster {
        &self.cluster
    }

    fn member_indices(&self) -> &[usize] {
        &self.member_indices
    }
}

/// Attaches to each cluster the (ascending) indices of the records carrying one of its
/// CDR3βs, so scoring looks members up directly instead of scanning every record.
fn with_member_indices(clusters: Vec<Cluster>, records: &[ChainRecord]) -> Vec<GliphCluster> {
    let mut by_beta: HashMap<&str, Vec<usize>> = HashMap::new();
    for (index, record) in records.iter().enumerate() {
        by_beta.entry(record.beta.as_str()).or_default().push(index);
    }
    clusters
        .into_par_iter()
        .map(|cluster| {
            let mut member_indices: Vec<usize> = cluster
                .member_sequences
                .iter()
                .filter_map(|beta| by_beta.get(beta.as_str()))
                .flatten()
                .copied()
                .collect();
            member_indices.sort_unstable();
            member_indices.dedup();
            GliphCluster {
                cluster,
                member_indices,
            }
        })
        .collect()
}

impl ClusterFinder for Gliph {
    type Cluster = GliphCluster;

    fn prepare<'a>(&self, records: &'a [ChainRecord]) -> Result<Cow<'a, [ChainRecord]>> {
        let valid: Vec<ChainRecord> = records
            .iter()
            .filter(|record| is_valid_cdr3(&record.beta))
            .cloned()
            .collect();
        if valid.is_empty() {
            bail!("No valid CDR3b sequences remain after GLIPH sequence filtering");
        }
        Ok(Cow::Owned(valid))
    }

    fn find(
        &self,
        records: &[ChainRecord],
        reference: &ReferenceData,
    ) -> Result<Vec<GliphCluster>> {
        let filter = GlobalFilterConfig::new(0.05, 2, 0.1);
        let sample_sequences: Vec<String> = records
            .iter()
            .map(|record| record.beta.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let reference_sequences = reference.valid_sequences();
        if reference_sequences.len() <= sample_sequences.len() {
            bail!(
                "Reference repertoire has {} valid sequences, but the sample has {}; GLIPH requires a larger reference repertoire",
                reference_sequences.len(),
                sample_sequences.len()
            );
        }
        if reference.cdr3_length_frequencies.is_empty() || reference.v_gene_frequencies.is_empty() {
            bail!("Reference data is missing V-gene or CDR3-length frequency distributions");
        }

        let mut clusters = local_candidates(
            self.algorithm,
            &sample_sequences,
            &reference_sequences,
            self.sim_depth,
            &filter,
        );

        if self.algorithm == Algorithm::Gliph2 {
            let global = global_fisher_candidates(
                &sample_sequences,
                &reference_sequences,
                &filter,
                &self.sub_matrix,
            );
            clusters.extend(global);
            clusters = merge_gliph2_candidates(clusters, &sample_sequences);
        } else {
            let cutoff = if records.len() < 125 { 2 } else { 1 };
            clusters = gliph1_candidates(clusters, &sample_sequences, cutoff);
        }
        eprintln!(
            "Found {} {} groups after clustering",
            clusters.len(),
            self.algorithm
        );
        clusters = add_singleton_clusters(clusters, records, &reference_sequences);
        clusters.sort_unstable_by(|left, right| {
            left.fisher_score
                .total_cmp(&right.fisher_score)
                .then_with(|| left.pattern.cmp(&right.pattern))
                .then_with(|| left.kind.cmp(&right.kind))
        });
        Ok(with_member_indices(clusters, records))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tcr_cluster::model::{ChainRecord, Cluster, GlobalFilterConfig};
    use crate::tcr_cluster::scoring::expansion_score;
    use pretty_assertions::assert_eq;
    use rstest::{fixture, rstest};

    #[test]
    fn gliph2_global_filter_removes_insignificant_candidates() {
        let sample_sequences = vec!["CASSIRSSYEQYF".to_owned(), "CASSIRSSYEQAF".to_owned()];
        let reference_sequences = vec![
            "CASSIRSSYEQVF".to_owned(),
            "CASSIRSSYEQGF".to_owned(),
            "CASSIRSSYEQDF".to_owned(),
            "CASSIRSSYEQEF".to_owned(),
            "CASSIRSSYEQNF".to_owned(),
            "CASSIRSSYEQHF".to_owned(),
            "CASSIRSSYEQKF".to_owned(),
            "CASSIRSSYEQRF".to_owned(),
        ];
        let config = GlobalFilterConfig::new(0.01, 2, 1.0);

        let clusters = global_fisher_candidates(
            &sample_sequences,
            &reference_sequences,
            &config,
            &Some(AminoAcidSubMatrix::Blosum62),
        );

        assert!(
            clusters.is_empty(),
            "common reference-supported global tags should not pass the Fisher threshold"
        );
    }

    #[test]
    fn gliph2_global_filter_removes_overly_common_reference_tags() {
        let sample_sequences = vec!["CASSIRSSYEQYF".to_owned(), "CASSIRSSYEQAF".to_owned()];
        let reference_sequences = vec![
            "CASSIRSSYEQVF".to_owned(),
            "CASSIRSSYEQGF".to_owned(),
            "CAAAAAAAAEQYF".to_owned(),
            "CAAAAAAAAEQAF".to_owned(),
        ];
        let config = GlobalFilterConfig::new(1.0, 2, 0.0);

        let clusters = global_fisher_candidates(
            &sample_sequences,
            &reference_sequences,
            &config,
            &Some(AminoAcidSubMatrix::Blosum62),
        );

        assert!(
            clusters.is_empty(),
            "tags present in the reference should be filtered when max_reference_fraction is zero"
        );
    }

    #[fixture]
    fn reference_sequences() -> Vec<String> {
        (0..300)
            .map(|index| {
                let amino_acids = b"ACDEFGHIKLMNPQTVWY";
                let first = amino_acids[index % amino_acids.len()] as char;
                let second = amino_acids[(index / amino_acids.len()) % amino_acids.len()] as char;
                format!("CAAAAA{first}{second}F")
            })
            .collect()
    }

    #[rstest]
    fn gliph2_global_filter_keeps_significant_rare_candidates(reference_sequences: Vec<String>) {
        let sample_sequences = vec![
            "CASSIRSSYEQYF".to_owned(),
            "CASSIRSSYEQAF".to_owned(),
            "CASSIRSSYEQVF".to_owned(),
        ];

        let config = GlobalFilterConfig::new(0.01, 2, 0.01);

        let clusters = global_fisher_candidates(
            &sample_sequences,
            &reference_sequences,
            &config,
            &Some(AminoAcidSubMatrix::Blosum62),
        );

        assert!(
            clusters
                .iter()
                .any(|cluster| cluster.member_sequences.len() == 3),
            "rare enriched sample-only global tags should pass filtering"
        );
        assert!(
            clusters
                .iter()
                .all(|cluster| cluster.fisher_score <= config.fisher_threshold)
        );
    }
    #[test]
    fn gliph2_merge_uses_minimal_local_pattern_label() {
        let sample_sequences = vec![
            "CASSQFLAAQTQYF".to_owned(),
            "CASSQFLGGDTQYF".to_owned(),
            "CASSQFLSSYEQYF".to_owned(),
        ];
        let candidates = vec![Cluster {
            pattern: "SQFL".to_owned(),
            kind: ClusterKind::Local,
            fisher_score: 0.001,
            member_sequences: sample_sequences.clone(),
        }];

        let merged = merge_gliph2_candidates(candidates, &sample_sequences);

        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].pattern, "SQFL");
    }

    #[test]
    fn gliph2_merge_prefers_longest_local_pattern_label() {
        let sample_sequences = vec![
            "CASSQFLAAQTQYF".to_owned(),
            "CASSQFLGGDTQYF".to_owned(),
            "CASSQFLSSYEQYF".to_owned(),
        ];
        let candidates = vec![
            Cluster {
                pattern: "QF".to_owned(),
                kind: ClusterKind::Local,
                fisher_score: 0.01,
                member_sequences: sample_sequences.clone(),
            },
            Cluster {
                pattern: "SQFL".to_owned(),
                kind: ClusterKind::Local,
                fisher_score: 0.001,
                member_sequences: sample_sequences.clone(),
            },
        ];

        let merged = merge_gliph2_candidates(candidates, &sample_sequences);

        assert_eq!(merged.len(), 2);
        assert_eq!(merged[0].pattern, "QF");
        assert_eq!(merged[1].pattern, "SQFL");
        assert_eq!(merged[0].fisher_score, 0.01);
        assert_eq!(merged[1].fisher_score, 0.001);
    }

    #[test]
    fn gliph2_merge_uses_trimmed_global_pattern_when_no_local_pattern_exists() {
        let sample_sequences = vec!["CASSQFLAAQTQYF".to_owned(), "CASSQFLGAQTQYF".to_owned()];
        let candidates = vec![Cluster {
            pattern: "SQFL%AQT".to_owned(),
            kind: ClusterKind::Global,
            fisher_score: 0.001,
            member_sequences: sample_sequences.clone(),
        }];

        let merged = merge_gliph2_candidates(candidates, &sample_sequences);

        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].pattern, "SQFL%AQT");
    }

    #[rstest]
    fn gliph2_global_filter_honors_min_members(reference_sequences: Vec<String>) {
        let sample_sequences = vec![
            "CASSIRSSYEQYF".to_owned(),
            "CASSIRSSYEQAF".to_owned(),
            "CASSIRSSYEQVF".to_owned(),
        ];

        let config = GlobalFilterConfig::new(1.0, 4, 1.0);

        let clusters = global_fisher_candidates(
            &sample_sequences,
            &reference_sequences,
            &config,
            &Some(AminoAcidSubMatrix::Blosum62),
        );

        assert!(
            clusters.is_empty(),
            "global candidates smaller than min_members should be filtered"
        );
    }

    #[test]
    fn expansion_score_is_neutral_without_counts() {
        let records = vec![
            ChainRecord {
                beta: "CASSAAAAEQYF".to_owned(),
                ..ChainRecord::default()
            },
            ChainRecord {
                beta: "CASSBBBBEQYF".to_owned(),
                ..ChainRecord::default()
            },
        ];
        let members: HashSet<&str> = ["CASSAAAAEQYF"].into_iter().collect::<HashSet<_>>();
        let exp_ctx = ExpansionContext::build(&records);

        let score = expansion_score(&members, &exp_ctx, 100, "no-counts");

        assert_eq!(score, 1.0);
    }

    #[test]
    fn expansion_score_aggregates_counts_by_unique_cdr3() {
        let records = vec![
            ChainRecord {
                beta: "CASSAAAAEQYF".to_owned(),
                counts: Some(10.0),
                ..ChainRecord::default()
            },
            ChainRecord {
                beta: "CASSAAAAEQYF".to_owned(),
                counts: Some(15.0),
                ..ChainRecord::default()
            },
            ChainRecord {
                beta: "CASSBBBBEQYF".to_owned(),
                counts: Some(1.0),
                ..ChainRecord::default()
            },
            ChainRecord {
                beta: "CASSCCCCEQYF".to_owned(),
                counts: Some(1.0),
                ..ChainRecord::default()
            },
        ];
        let members: HashSet<&str> = ["CASSAAAAEQYF"].into_iter().collect::<HashSet<_>>();
        let exp_ctx = ExpansionContext::build(&records);

        let score = expansion_score(&members, &exp_ctx, 100, "aggregated-counts");

        assert!(
            score < 1.0,
            "expanded unique CDR3 should be more extreme than at least some background draws"
        );
    }

    #[test]
    fn expansion_score_is_deterministic() {
        let records = vec![
            ChainRecord {
                beta: "CASSAAAAEQYF".to_owned(),
                counts: Some(20.0),
                ..ChainRecord::default()
            },
            ChainRecord {
                beta: "CASSBBBBEQYF".to_owned(),
                counts: Some(1.0),
                ..ChainRecord::default()
            },
            ChainRecord {
                beta: "CASSCCCCEQYF".to_owned(),
                counts: Some(1.0),
                ..ChainRecord::default()
            },
            ChainRecord {
                beta: "CASSDDDDEQYF".to_owned(),
                counts: Some(1.0),
                ..ChainRecord::default()
            },
        ];
        let members: HashSet<&str> = ["CASSAAAAEQYF"].into_iter().collect::<HashSet<_>>();
        let exp_ctx = ExpansionContext::build(&records);

        assert_eq!(
            expansion_score(&members, &exp_ctx, 100, "deterministic"),
            expansion_score(&members, &exp_ctx, 100, "deterministic")
        );
    }
}
