use crate::tcr_cluster::blosum::{AminoAcidSubMatrix, blosum_nonnegative};
use crate::tcr_cluster::model::{ChainRecord, Cluster, ClusterKind, HasCluster};
use anyhow::Context;
use indicatif::ProgressBar;
use rayon::prelude::*;
use rust_igraph::{Graph, connected_components};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::hash::{Hash, Hasher};

const K_MIN: usize = 3;
const K_MAX: usize = 8;
const MIN_SUPPORT: f64 = 0.5;

#[derive(Debug, Clone)]
pub struct BlosumMotif {
    pub consensus: String, // e.g., "G%YG" (if degenerate) or "GGYG" (if exact)
    pub degenerate_motif: String, // e.g., "G[G/S]YG"
    pub variable_pos: usize,
    pub support_count: usize,
    pub support_fraction: f64,
    pub matched_kmers: BTreeMap<String, usize>,
}

impl Hash for BlosumMotif {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.consensus.hash(state);
        self.degenerate_motif.hash(state);
        self.variable_pos.hash(state);
        self.support_count.hash(state);
        self.support_fraction.to_bits().hash(state);
        self.matched_kmers.hash(state);
    }
}

impl BlosumMotif {
    pub fn new(
        consensus: String,
        degenerate_motif: String,
        variable_pos: usize,
        matched_kmers: BTreeMap<String, usize>,
        supporing_seq_ids: Vec<usize>,
        sequences: &[String],
    ) -> Self {
        let support_count = supporing_seq_ids.len();
        let support_fraction = support_count as f64 / sequences.len() as f64;

        BlosumMotif {
            consensus,
            degenerate_motif,
            variable_pos,
            support_count,
            support_fraction,
            matched_kmers,
        }
    }

    pub fn consensus(self) -> String {
        self.consensus
    }
}

impl PartialEq<Self> for BlosumMotif {
    fn eq(&self, other: &Self) -> bool {
        self.consensus == other.consensus
            && self.degenerate_motif == other.degenerate_motif
            && self.variable_pos == other.variable_pos
            && self.support_count == other.support_count
            && self.support_fraction == other.support_fraction
            && self.matched_kmers == other.matched_kmers
    }
}

impl Eq for BlosumMotif {}

struct Candidate {
    motif: BlosumMotif,
    center_freq: usize,
    seq_ids: Vec<usize>,
}

impl Candidate {
    pub fn new(
        consensus: String,
        degenerate_motif: String,
        variable_pos: usize,
        matched_kmers: BTreeMap<String, usize>,
        supporing_seq_ids: Vec<usize>,
        sequences: &[String],
        center_kmer: &String,
    ) -> Self {
        let center_freq = matched_kmers.get(center_kmer).copied().unwrap_or(0);
        let motif = BlosumMotif::new(
            consensus,
            degenerate_motif,
            variable_pos,
            matched_kmers,
            supporing_seq_ids.clone(),
            sequences,
        );
        Candidate {
            motif,
            center_freq,
            seq_ids: supporing_seq_ids,
        }
    }
}

pub fn trim_cdr3(seq: &str, trim_left: usize, trim_right: usize) -> &str {
    if seq.len() <= trim_left + trim_right {
        seq
    } else {
        &seq[trim_left..seq.len() - trim_right]
    }
}

pub fn discover_blosum_motifs(
    sequences: &[String],
    sub_matrix: &AminoAcidSubMatrix,
    k: usize,
    trim_flanks: bool,
    top_n: usize,
) -> Vec<BlosumMotif> {
    if sequences.is_empty() || k == 0 {
        return Vec::new();
    }

    let (kmer_to_seqs, wildcard_buckets) = wildcard_motif(sequences, k, trim_flanks);

    let mut candidates: Vec<Candidate> = Vec::new();
    let mut seen_clusters: HashSet<(String, Vec<usize>)> = HashSet::new();

    for center_kmer in kmer_to_seqs.keys() {
        let center_bytes = center_kmer.as_bytes();

        for pos in 0..k {
            let mut pattern = String::with_capacity(k);
            pattern.push_str(&center_kmer[..pos]);
            pattern.push('.');
            pattern.push_str(&center_kmer[pos + 1..]);

            let center_aa = center_bytes[pos];
            let Some(bucket) = wildcard_buckets.get(&pattern) else {
                continue;
            };

            let matched_kmers: Vec<&String> = bucket
                .iter()
                .filter(|cand| {
                    let cand_aa = cand.as_bytes()[pos];
                    cand_aa == center_aa || blosum_nonnegative(center_aa, cand_aa, sub_matrix)
                })
                .collect();

            let mut supporting_seqs = BTreeSet::new();
            let mut kmer_counts = BTreeMap::new();

            for &cand in &matched_kmers {
                if let Some(seqs) = kmer_to_seqs.get(cand) {
                    supporting_seqs.extend(seqs.iter().copied());
                    kmer_counts.insert(cand.clone(), seqs.len());
                }
            }

            let seq_ids: Vec<usize> = supporting_seqs.into_iter().collect();
            let cluster_key = (center_kmer.clone(), seq_ids.clone());
            if !seen_clusters.insert(cluster_key) {
                continue;
            }

            let mut observed_aas: Vec<u8> =
                matched_kmers.iter().map(|c| c.as_bytes()[pos]).collect();
            observed_aas.sort_unstable();
            observed_aas.dedup();

            observed_aas.sort_by(|&aa_a, &aa_b| {
                let mut kmer_a = center_kmer.clone().into_bytes();
                kmer_a[pos] = aa_a;
                let freq_a = kmer_to_seqs
                    .get(std::str::from_utf8(&kmer_a).unwrap_or(""))
                    .map_or(0, |s| s.len());

                let mut kmer_b = center_kmer.clone().into_bytes();
                kmer_b[pos] = aa_b;
                let freq_b = kmer_to_seqs
                    .get(std::str::from_utf8(&kmer_b).unwrap_or(""))
                    .map_or(0, |s| s.len());

                freq_b.cmp(&freq_a).then_with(|| aa_a.cmp(&aa_b))
            });

            // Format both `consensus` ("G%YG" if >1 AA, else "GGYG") and `degenerate_motif` ("G[G/S]YG")
            let (consensus, degenerate_motif) = if observed_aas.len() > 1 {
                let joined = observed_aas
                    .iter()
                    .map(|&b| (b as char).to_string())
                    .collect::<Vec<_>>()
                    .join("/");
                let deg = format!(
                    "{}\\[{}\\]
{}",
                    &center_kmer[..pos],
                    joined,
                    &center_kmer[pos + 1..]
                );
                let cons = format!("{}%{}", &center_kmer[..pos], &center_kmer[pos + 1..]);
                (cons, deg)
            } else {
                (center_kmer.clone(), center_kmer.clone())
            };

            candidates.push(Candidate::new(
                consensus,
                degenerate_motif,
                pos,
                kmer_counts,
                seq_ids,
                sequences,
                center_kmer,
            ));
        }
    }

    candidates.sort_by(|a, b| {
        b.motif
            .support_count
            .cmp(&a.motif.support_count)
            .then_with(|| b.center_freq.cmp(&a.center_freq))
            .then_with(|| {
                b.motif
                    .matched_kmers
                    .len()
                    .cmp(&a.motif.matched_kmers.len())
            })
            .then_with(|| a.motif.consensus.cmp(&b.motif.consensus))
    });

    let mut final_results = Vec::with_capacity(top_n);
    let mut seen_seq_sets: HashSet<Vec<usize>> = HashSet::new();

    for cand in candidates {
        if seen_seq_sets.insert(cand.seq_ids) {
            final_results.push(cand.motif);
            if final_results.len() >= top_n {
                break;
            }
        }
    }

    final_results
}

/// Finds the longest motif (k from `k_max` down to `k_min`) shared by at least
/// `min_support` of the sequences, so it stays representative of the cluster.
/// Falls back to the best-supported motif at `k_min` if none qualifies.
pub fn select_longest_representative_motif(
    sequences: &[String],
    sub_matrix: &AminoAcidSubMatrix,
    trim_flanks: bool,
) -> Option<BlosumMotif> {
    for k in (K_MIN..=K_MAX).rev() {
        let best = discover_blosum_motifs(sequences, sub_matrix, k, trim_flanks, 1)
            .into_iter()
            .next();
        if let Some(motif) = best
            && (motif.support_fraction >= MIN_SUPPORT || k == K_MIN)
        {
            return Some(motif);
        }
    }
    None
}

fn wildcard_motif(
    sequences: &[String],
    k: usize,
    trim_flanks: bool,
) -> (
    HashMap<String, BTreeSet<usize>>,
    HashMap<String, HashSet<String>>,
) {
    let mut kmer_to_seqs: HashMap<String, BTreeSet<usize>> = HashMap::new();
    let mut wildcard_buckets: HashMap<String, HashSet<String>> = HashMap::new();

    for (seq_idx, raw_seq) in sequences.iter().enumerate() {
        let seq = if trim_flanks {
            trim_cdr3(raw_seq, 3, 3)
        } else {
            raw_seq.as_str()
        };

        if seq.len() < k {
            continue;
        }

        for i in 0..=(seq.len() - k) {
            let kmer = &seq[i..i + k];
            kmer_to_seqs
                .entry(kmer.to_string())
                .or_default()
                .insert(seq_idx);

            for pos in 0..k {
                let mut pattern = String::with_capacity(k);
                pattern.push_str(&kmer[..pos]);
                pattern.push('.');
                pattern.push_str(&kmer[pos + 1..]);

                wildcard_buckets
                    .entry(pattern)
                    .or_default()
                    .insert(kmer.to_string());
            }
        }
    }
    (kmer_to_seqs, wildcard_buckets)
}

/// A cluster from the joint alpha/beta distance graph: the shared [`Cluster`] (built from the
/// beta chain) plus a motif for each chain.
#[derive(Debug, Clone)]
pub struct PairedCluster {
    pub threshold: u8,
    pub cluster_id: u32,
    pub member_indices: Vec<usize>,
    pub alpha_motif: Option<BlosumMotif>,
    pub beta_motif: Option<BlosumMotif>,
    pub cluster: Cluster,
}

impl PartialEq<Self> for PairedCluster {
    fn eq(&self, other: &Self) -> bool {
        self.threshold == other.threshold
            && self.cluster_id == other.cluster_id
            && self.member_indices == other.member_indices
            && self.alpha_motif == other.alpha_motif
            && self.beta_motif == other.beta_motif
            && self.cluster == other.cluster
    }
}

impl Eq for PairedCluster {}

impl Hash for PairedCluster {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.threshold.hash(state);
        self.cluster_id.hash(state);
        self.member_indices.hash(state);
        self.alpha_motif.hash(state);
        self.beta_motif.hash(state);
        self.cluster.hash(state);
    }
}

impl HasCluster for PairedCluster {
    fn cluster(&self) -> &Cluster {
        &self.cluster
    }

    fn member_indices(&self) -> &[usize] {
        &self.member_indices
    }
}

pub fn build_clusters_for_threshold(
    threshold: u8,
    memberships: &[u32],
    records: &[ChainRecord],
    sub_matrix: &AminoAcidSubMatrix,
    trim_flanks: bool,
    fisher_for: &(dyn Fn(&str) -> f64 + Sync),
    cluster_kind: &ClusterKind,
    max_members: usize,
    pb: &ProgressBar,
) -> anyhow::Result<Vec<PairedCluster>> {
    // Group record indices by their cluster component ID
    let mut groups: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for (record_idx, &cluster_id) in memberships.iter().enumerate() {
        groups.entry(cluster_id).or_default().push(record_idx);
    }

    // Skip singletons and uninformative giant components before any expensive
    // deduplication or motif discovery.
    let mut groups: Vec<(u32, Vec<usize>)> = groups
        .into_iter()
        .filter(|(_, members)| members.len() >= 2 && members.len() < max_members)
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
            let unique_alphas: Vec<String> = member_indices
                .iter()
                .map(|&idx| records[idx].alpha.clone())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .flatten()
                .collect();

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

            let alpha_motif = if unique_alphas.len() >= 2 {
                select_longest_representative_motif(&unique_alphas, sub_matrix, trim_flanks)
            } else {
                None
            };

            let beta_motif =
                select_longest_representative_motif(&unique_betas, sub_matrix, trim_flanks);
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
            let unique_alphas: Vec<String> = member_indices
                .iter()
                .filter_map(|&idx| records[idx].alpha.clone())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
            let alpha_motif = if unique_alphas.len() >= 2 {
                select_longest_representative_motif(&unique_alphas, sub_matrix, trim_flanks)
            } else {
                None
            };

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

    #[test]
    fn picks_longest_representative_motif() {
        let seqs: Vec<String> = ["CASSLGGYGYTF", "CASRLGGYGETF", "CASTLGGYGQYF"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let m = select_longest_representative_motif(&seqs, &AminoAcidSubMatrix::Blosum62, false)
            .unwrap();
        println!("Motif: {:?}", m);
        assert!(m.consensus.len() >= 5);
        assert!(m.support_fraction >= MIN_SUPPORT);
    }
}
