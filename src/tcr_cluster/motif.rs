//! BLOSUM-tolerant CDR3 motif discovery shared by the cluster finders.
//!
//! A motif is the longest k-mer (with at most one BLOSUM-similar degenerate position)
//! shared by enough of a cluster's unique CDR3s; GLIPH and tcrdist both use it to label
//! each cluster's CDR3α.

use crate::tcr_cluster::blosum::{AminoAcidSubMatrix, Chain, blosum_nonnegative};
use crate::tcr_cluster::model::ChainRecord;
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
    chain: &Chain,
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
                    cand_aa == center_aa
                        || blosum_nonnegative(center_aa, cand_aa, sub_matrix, chain)
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
    chain: &Chain,
) -> Option<BlosumMotif> {
    for k in (K_MIN..=K_MAX).rev() {
        let best = discover_blosum_motifs(sequences, sub_matrix, k, trim_flanks, 1, chain)
            .into_iter()
            .next();
        if let Some(motif) = best
            && (motif.support_fraction >= MIN_SUPPORT)
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

/// Motif of the CDR3α sequences of `records[i]` for each of `member_indices`.
///
/// Empty and missing CDR3αs are ignored. A single distinct CDR3α is its own motif (as a
/// single CDR3β is for single clusters); two or more go through
/// [`select_longest_representative_motif`].
pub fn cluster_alpha_motif(
    records: &[ChainRecord],
    member_indices: &[usize],
    sub_matrix: &AminoAcidSubMatrix,
    trim_flanks: bool,
) -> Option<BlosumMotif> {
    let unique_alphas: Vec<String> = member_indices
        .iter()
        .filter_map(|&idx| records[idx].alpha.as_deref())
        .filter(|alpha| !alpha.is_empty())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .map(str::to_owned)
        .collect();
    match unique_alphas.as_slice() {
        [] => None,
        [alpha] => Some(BlosumMotif::new(
            alpha.clone(),
            alpha.clone(),
            0,
            BTreeMap::from([(alpha.clone(), 1)]),
            vec![0],
            &unique_alphas,
        )),
        _ => select_longest_representative_motif(
            &unique_alphas,
            sub_matrix,
            trim_flanks,
            &Chain::Alpha,
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn picks_longest_representative_motif() {
        let seqs: Vec<String> = ["CASSLGGYGYTF", "CASRLGGYGETF", "CASTLGGYGQYF"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let m = select_longest_representative_motif(
            &seqs,
            &AminoAcidSubMatrix::Blosum62,
            false,
            &Chain::Beta,
        )
        .unwrap();
        println!("Motif: {:?}", m);
        assert!(m.consensus.len() >= 5);
        assert!(m.support_fraction >= MIN_SUPPORT);
    }

    fn alpha_record(alpha: Option<&str>) -> ChainRecord {
        ChainRecord {
            alpha: alpha.map(str::to_owned),
            ..ChainRecord::default()
        }
    }

    #[test]
    fn cluster_alpha_motif_handles_missing_single_and_shared_alphas() {
        let records = [
            alpha_record(None),
            alpha_record(Some("")),
            alpha_record(Some("CAVRDSNYQLIW")),
            alpha_record(Some("CAVRDSNYQLIW")),
            alpha_record(Some("CAVRDTNYQLIW")),
            alpha_record(Some("CAVKDSNYQLIW")),
        ];
        let motif = |indices: &[usize]| {
            cluster_alpha_motif(&records, indices, &AminoAcidSubMatrix::Blosum62, true)
                .map(|m| m.consensus)
        };

        assert_eq!(motif(&[0, 1]), None, "no CDR3α, no motif");
        assert_eq!(
            motif(&[0, 1, 2, 3]).as_deref(),
            Some("CAVRDSNYQLIW"),
            "a single distinct CDR3α is its own motif"
        );
        let shared = motif(&[2, 4, 5]).expect("similar CDR3αs share a motif");
        assert!(shared.len() >= K_MIN, "motif {shared}");
    }
}
