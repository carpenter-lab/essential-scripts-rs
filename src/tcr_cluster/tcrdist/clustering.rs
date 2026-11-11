use crate::tcr_cluster::blosum::{AminoAcidSubMatrix, Chain};
use crate::tcr_cluster::model::{ChainRecord, Cluster, ClusterKind, HasCluster};
use crate::tcr_cluster::motif::{
    BlosumMotif, cluster_alpha_motif, select_longest_representative_motif,
};
use anyhow::Context;
use indicatif::ProgressBar;
use rayon::prelude::*;
use rust_igraph::{Graph, connected_components};
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::hash::{Hash, Hasher};

/// A cluster from the joint alpha/beta distance graph: the shared [`Cluster`] (built from the
/// beta chain) plus a motif for each chain.
#[derive(Debug, Clone)]
pub struct PairedCluster {
    /// The CDR3-edit threshold the cluster was found at (0 for identical-sequence clusters).
    pub threshold: u8,
    pub cluster_id: u32,
    pub member_indices: Vec<usize>,
    /// Written as the `pattern_alpha` column.
    pub alpha_motif: Option<BlosumMotif>,
    // `cluster.pattern` already carries the β consensus.
    #[allow(dead_code)]
    pub beta_motif: Option<BlosumMotif>,
    pub cluster: Cluster,
}

/// Two clusters are the same when they have the same members: everything else (motifs,
/// pattern, Fisher score) is derived from the members, and the threshold and component ID
/// only say where the cluster was found. Collecting thresholds in increasing order into a
/// `HashSet` therefore keeps each distinct cluster once, at the smallest threshold.
impl PartialEq<Self> for PairedCluster {
    fn eq(&self, other: &Self) -> bool {
        self.member_indices == other.member_indices
    }
}

impl Eq for PairedCluster {}

impl Hash for PairedCluster {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.member_indices.hash(state);
    }
}

impl HasCluster for PairedCluster {
    fn cluster(&self) -> &Cluster {
        &self.cluster
    }

    fn member_indices(&self) -> &[usize] {
        &self.member_indices
    }

    fn alpha_pattern(&self) -> Option<&str> {
        self.alpha_motif
            .as_ref()
            .map(|motif| motif.consensus.as_str())
    }
}

#[allow(clippy::too_many_arguments)]
pub fn build_clusters_for_threshold(
    threshold: u8,
    memberships: &[u32],
    records: &[ChainRecord],
    sub_matrix: &AminoAcidSubMatrix,
    trim_flanks: bool,
    fisher_for: &(dyn Fn(&str) -> f64 + Sync),
    cluster_kind: &ClusterKind,
    pb: &ProgressBar,
) -> anyhow::Result<Vec<PairedCluster>> {
    // Group record indices by their cluster component ID
    let mut groups: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for (record_idx, &cluster_id) in memberships.iter().enumerate() {
        groups.entry(cluster_id).or_default().push(record_idx);
    }

    // Skip singletons before any expensive deduplication or motif discovery.
    let mut groups: Vec<(u32, Vec<usize>)> = groups
        .into_iter()
        .filter(|(_, members)| members.len() >= 2)
        .collect();
    // Largest clusters first so the slowest motif searches start early instead of
    // becoming the tail of the parallel run.
    groups.sort_unstable_by_key(|(_, members)| std::cmp::Reverse(members.len()));

    // Each cluster is independent; `with_max_len(1)` lets rayon balance clusters of very
    // different sizes one at a time instead of in fixed chunks.
    // `set_length` alone keeps the previous threshold's position, so reset first.
    pb.reset();
    pb.set_length(groups.len() as u64);
    pb.set_message("Finding motifs");
    let discovered: Vec<PairedCluster> = groups
        .into_par_iter()
        .with_max_len(1)
        .filter_map(|(cluster_id, member_indices)| {
            // Deduplicate sequences to avoid PCR duplicate bias
            let unique_betas: Vec<String> = member_indices
                .iter()
                .map(|&idx| records[idx].beta.clone())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();

            // Components with a single CDR3β are covered by `build_exact_beta_clusters`,
            // which groups every record sharing that CDR3β regardless of the graph.
            if unique_betas.len() < 2 {
                pb.inc(1);
                return None;
            }

            let alpha_motif =
                cluster_alpha_motif(records, &member_indices, sub_matrix, trim_flanks);
            let beta_motif = select_longest_representative_motif(
                &unique_betas,
                sub_matrix,
                trim_flanks,
                &Chain::Beta,
            );
            let kind = *cluster_kind;

            let (pattern, fisher_score) = match beta_motif.clone() {
                Some(motif) => {
                    let p = motif.consensus();
                    let f = fisher_for(&p);
                    (p, f)
                }
                None => ("no motif found".to_string(), 1.0),
            };

            let cluster = Cluster {
                pattern,
                kind,
                fisher_score,
                member_sequences: unique_betas,
            };
            pb.inc(1);

            Some(PairedCluster {
                threshold,
                cluster_id,
                member_indices,
                alpha_motif,
                beta_motif,
                cluster,
            })
        })
        .collect();

    Ok(discovered)
}

/// [`ClusterKind::Single`] clusters of identical CDR3β (or CDR3α + CDR3β when `use_alpha`)
/// sequences, independent of the distance graph so V-gene differences do not split them.
/// Every sequence shared by two or more records gets one, plus each remaining record whose
/// index is not in `covered` (the records already in another cluster).
/// `fisher_for` receives the CDR3β and the number of records carrying it.
pub fn build_exact_beta_clusters(
    records: &[ChainRecord],
    sub_matrix: &AminoAcidSubMatrix,
    trim_flanks: bool,
    use_alpha: bool,
    covered: &HashSet<usize>,
    fisher_for: &(dyn Fn(&str, usize) -> f64 + Sync),
) -> Vec<PairedCluster> {
    let mut groups: BTreeMap<(&str, Option<&str>), Vec<usize>> = BTreeMap::new();
    for (idx, record) in records.iter().enumerate() {
        let alpha = record.alpha.as_deref().filter(|_| use_alpha);
        groups
            .entry((record.beta.as_str(), alpha))
            .or_default()
            .push(idx);
    }

    groups
        .into_iter()
        // Shared sequences always form a cluster; a sequence seen once only does when no other
        // cluster contains it, so every record appears in the output at least once.
        .filter(|(_, members)| members.len() >= 2 || !covered.contains(&members[0]))
        .enumerate()
        .map(|(cluster_id, ((beta, _), member_indices))| {
            let alpha_motif =
                cluster_alpha_motif(records, &member_indices, sub_matrix, trim_flanks);

            let members = vec![beta.to_string(); member_indices.len()];
            let beta_motif = BlosumMotif::new(
                beta.to_string(),
                beta.to_string(),
                0,
                BTreeMap::from([(beta.to_string(), member_indices.len())]),
                (0..member_indices.len()).collect(),
                &members,
            );
            let cluster = Cluster::new(
                beta.to_string(),
                ClusterKind::Single,
                fisher_for(beta, member_indices.len()),
                vec![beta.to_string()],
            );
            PairedCluster {
                threshold: 0,
                cluster_id: cluster_id as u32,
                member_indices,
                alpha_motif,
                beta_motif: Some(beta_motif),
                cluster,
            }
        })
        .collect()
}

/// Component label of every vertex of `graph`.
pub fn graph_components(graph: &Graph) -> anyhow::Result<Vec<u32>> {
    connected_components(graph)
        .context("Failed to compute connected components in the graph")
        .map(|cc| cc.membership)
}

pub fn graph_louvain(_graph: &Graph) -> anyhow::Result<Vec<u32>> {
    anyhow::bail!("NotImplemented: Louvain clustering")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(beta: &str, alpha: Option<&str>) -> ChainRecord {
        ChainRecord {
            beta: beta.to_owned(),
            alpha: alpha.map(str::to_owned),
            ..ChainRecord::default()
        }
    }

    #[test]
    fn same_members_at_a_larger_threshold_is_a_duplicate() {
        let cluster = |threshold: u8, cluster_id: u32, members: Vec<usize>| PairedCluster {
            threshold,
            cluster_id,
            member_indices: members,
            alpha_motif: None,
            beta_motif: None,
            cluster: Cluster {
                pattern: format!("p{threshold}"),
                kind: ClusterKind::Global,
                fisher_score: 1.0,
                member_sequences: Vec::new(),
            },
        };
        let found: HashSet<PairedCluster> = [
            cluster(1, 7, vec![0, 1]),
            cluster(2, 3, vec![0, 1]),
            cluster(2, 4, vec![0, 1, 2]),
        ]
        .into_iter()
        .collect();
        assert_eq!(found.len(), 2);
        let pair = found.iter().find(|c| c.member_indices == [0, 1]).unwrap();
        assert_eq!(pair.threshold, 1, "the smallest threshold is kept");
    }

    #[test]
    fn exact_beta_clusters_group_identical_cdr3b() {
        let records = [
            rec("CASSIRSSYEQYF", Some("CAVRDSNYQLIW")),
            rec("CASSLGGYGYTF", None),
            rec("CASSIRSSYEQYF", Some("CAVKKTTGGNKLTF")),
            rec("CASSIRSSYEQYF", Some("CAVRDSNYQLIW")),
            rec("CASRLGGYGETF", None),
        ];
        // Records 1 and 4 (seen once) are already in another cluster.
        let covered = HashSet::from([1, 4]);
        let exact = |use_alpha| {
            build_exact_beta_clusters(
                &records,
                &AminoAcidSubMatrix::Blosum62,
                true,
                use_alpha,
                &covered,
                &|_, _| 0.5,
            )
        };
        let paired = exact(true);
        let members: Vec<&[usize]> = paired.iter().map(|c| c.member_indices.as_slice()).collect();
        assert_eq!(
            members,
            [&[2][..], &[0, 3][..]],
            "differing CDR3α splits the group; the uncovered leftover is its own single"
        );
        assert_eq!(
            paired[1].alpha_pattern(),
            Some("CAVRDSNYQLIW"),
            "a shared CDR3α is the alpha pattern"
        );

        let clusters = exact(false);
        assert_eq!(
            clusters.len(),
            1,
            "covered seen-once records get no extra single"
        );
        let c = &clusters[0];
        assert_eq!(c.member_indices, vec![0, 2, 3]);
        assert_eq!(c.cluster.kind, ClusterKind::Single);
        assert_eq!(c.cluster.pattern, "CASSIRSSYEQYF");
        assert_eq!(
            c.cluster.member_sequences,
            vec!["CASSIRSSYEQYF".to_string()]
        );
        assert_eq!(c.beta_motif.as_ref().unwrap().support_count, 3);
    }
}
