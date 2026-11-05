use crate::progress::progress_bar_with_message;
use crate::tcr_cluster::{Algorithm, MotifStats};
use indicatif::ParallelProgressIterator;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use rayon::slice::ParallelSliceMut;
use rustc_hash::FxHashMap;
use std::collections::BTreeSet;
use std::sync::OnceLock;

const MOTIF_LENGTHS: [usize; 4] = [2, 3, 4, 5];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ContingencyTable2x2 {
    pub a: usize, // Query, motif present
    pub b: usize, // Query, motif absent (N_query - a)
    pub c: usize, // Reference, motif present
    pub d: usize, // Reference, motif absent (N_ref - c)
}

impl ContingencyTable2x2 {
    /// Constructs a table from published parameters:
    /// - `a`: query motif hits
    /// - `n_query`: total unique CDR3s in query set
    /// - `c`: reference motif hits
    /// - `n_ref`: total unique CDR3s in reference set
    pub fn from_reported_counts(
        a: usize,
        n_query: usize,
        c: usize,
        n_ref: usize,
    ) -> Result<Self, &'static str> {
        if a > n_query {
            return Err("Query motif count cannot exceed total query CDR3s");
        }
        if c > n_ref {
            return Err("Reference motif count cannot exceed total reference CDR3s");
        }

        Ok(Self {
            a,
            b: n_query - a,
            c,
            d: n_ref - c,
        })
    }

    /// Sample Odds Ratio. Returns f64::INFINITY if denominator is zero.
    pub fn odds_ratio(&self) -> f64 {
        let num = (self.a as f64) * (self.d as f64);
        let den = (self.b as f64) * (self.c as f64);

        if den == 0.0 {
            if num == 0.0 { 1.0 } else { f64::INFINITY }
        } else {
            num / den
        }
    }
}

static GLOBAL_FACT_TABLE: OnceLock<LogFactorialTable> = OnceLock::new();

/// Precomputed log-factorial table: ln(n!)
/// Enables O(1) hypergeometric point probability evaluation.
pub struct LogFactorialTable {
    log_factorials: Vec<f64>,
}

impl LogFactorialTable {
    pub fn new(max_n: usize) -> Self {
        let mut log_factorials = Vec::with_capacity(max_n + 1);
        log_factorials.push(0.0); // ln(0!) = 0
        let mut acc = 0.0;
        for i in 1..=max_n {
            acc += (i as f64).ln();
            log_factorials.push(acc);
        }
        Self { log_factorials }
    }

    pub fn global(default_max_n: usize) -> &'static Self {
        GLOBAL_FACT_TABLE.get_or_init(|| LogFactorialTable::new(default_max_n))
    }

    #[inline]
    fn stirling_ln_fact(n: usize) -> f64 {
        let x = (n + 1) as f64; // Gamma(x) = n!
        let ln_sqrt_2pi = 0.918_938_533_204_672_8; // 0.5 * ln(2 * PI)
        let inv_x = 1.0 / x;
        let inv_x2 = inv_x * inv_x;

        // Ramanujan / Stirling asymptotic expansion for ln Gamma(x)
        (x - 0.5) * x.ln() - x
            + ln_sqrt_2pi
            + inv_x * (1.0 / 12.0 - inv_x2 * (1.0 / 360.0 - inv_x2 * (1.0 / 1260.0)))
    }

    #[inline(always)]
    pub fn ln_fact(&self, n: usize) -> f64 {
        if let Some(&val) = self.log_factorials.get(n) {
            val
        } else {
            Self::stirling_ln_fact(n)
        }
    }

    /// Computes ln(n choose k)
    #[inline(always)]
    pub fn ln_choose(&self, n: usize, k: usize) -> f64 {
        if k > n {
            f64::NEG_INFINITY
        } else {
            self.ln_fact(n) - self.ln_fact(k) - self.ln_fact(n - k)
        }
    }

    /// Calculates P(X = k) under the hypergeometric distribution:
    /// Total population = N, Total successes = M, Sample size = n.
    #[inline(always)]
    pub fn hypergeom_p(&self, k: usize, n: usize, m: usize, total: usize) -> f64 {
        if k > n || k > m || (n - k) > (total - m) {
            return 0.0;
        }
        let ln_p =
            self.ln_choose(m, k) + self.ln_choose(total - m, n - k) - self.ln_choose(total, n);
        ln_p.exp()
    }

    /// Calculates one-sided (enrichment) and two-sided p-values.
    pub fn fishers_exact(&self, table: &ContingencyTable2x2) -> (f64, f64) {
        let n_query = table.a + table.b;
        let n_ref = table.c + table.d;
        let total = n_query + n_ref;
        let m = table.a + table.c; // total motif-positive across both groups

        let min_k = (m + n_query).saturating_sub(total);
        let max_k = usize::min(n_query, m);

        let observed_a = table.a;
        let observed_prob = self.hypergeom_p(observed_a, n_query, m, total);

        let mut p_one_sided = 0.0;
        let mut p_two_sided = 0.0;
        let tolerance = 1e-12;

        for k in min_k..=max_k {
            let p_k = self.hypergeom_p(k, n_query, m, total);

            // One-sided: enrichment in query set (k >= observed_a)
            if k >= observed_a {
                p_one_sided += p_k;
            }

            // Two-sided: probabilities <= observed probability
            if p_k <= observed_prob * (1.0 + tolerance) {
                p_two_sided += p_k;
            }
        }

        (p_one_sided.clamp(0.0, 1.0), p_two_sided.clamp(0.0, 1.0))
    }
}

/// Sorted, de-duplicated borrowed views of `seqs`.
fn unique_sequences(seqs: &[String]) -> Vec<&str> {
    let mut unique: Vec<&str> = seqs.iter().map(String::as_str).collect();
    unique.par_sort_unstable();
    unique.dedup();
    unique
}

/// Fills `out` with the distinct motifs of `sequence`, borrowed from it (no allocation
/// beyond `out`'s reused capacity). Matches the set returned by [`motifs`].
fn distinct_motifs_into<'a>(sequence: &'a str, out: &mut Vec<&'a str>) {
    out.clear();
    let region = motif_region(sequence);
    for length in MOTIF_LENGTHS {
        if region.len() < length {
            continue;
        }
        out.extend((0..=region.len() - length).filter_map(|start| region.get(start..start + length)));
    }
    out.sort_unstable();
    out.dedup();
}

/// Motif -> unique sample sequences containing it, built in parallel.
fn sample_motif_members<'a>(unique_samples: &[&'a str]) -> FxHashMap<&'a str, Vec<&'a str>> {
    unique_samples
        .par_iter()
        .fold(
            || (FxHashMap::<&str, Vec<&str>>::default(), Vec::new()),
            |(mut map, mut scratch), &sequence| {
                distinct_motifs_into(sequence, &mut scratch);
                for &motif in &scratch {
                    map.entry(motif).or_default().push(sequence);
                }
                (map, scratch)
            },
        )
        .map(|(map, _)| map)
        .reduce(FxHashMap::default, |a, b| {
            let (mut big, small) = if a.len() >= b.len() { (a, b) } else { (b, a) };
            for (motif, members) in small {
                big.entry(motif).or_default().extend(members);
            }
            big
        })
}

/// Motif -> number of unique reference sequences containing it, restricted to motifs
/// present in `sample_members`. Built in parallel.
fn reference_motif_counts<'a>(
    unique_refs: &[&'a str],
    sample_members: &FxHashMap<&str, Vec<&str>>,
) -> FxHashMap<&'a str, usize> {
    unique_refs
        .par_iter()
        .fold(
            || (FxHashMap::<&str, usize>::default(), Vec::new()),
            |(mut map, mut scratch), &sequence| {
                distinct_motifs_into(sequence, &mut scratch);
                for &motif in &scratch {
                    if sample_members.contains_key(motif) {
                        *map.entry(motif).or_default() += 1;
                    }
                }
                (map, scratch)
            },
        )
        .map(|(map, _)| map)
        .reduce(FxHashMap::default, |a, b| {
            let (mut big, small) = if a.len() >= b.len() { (a, b) } else { (b, a) };
            for (motif, count) in small {
                *big.entry(motif).or_default() += count;
            }
            big
        })
}

/// Per-motif enrichment statistics of `sample_sequences` against `reference_sequences`.
///
/// Motif counting runs in parallel over borrowed `&str` slices; `String`s are only
/// allocated for the yielded motif and its members. Items are yielded in motif order.
pub fn motif_statistics_iter(
    sample_sequences: &[String],
    reference_sequences: &[String],
) -> impl ParallelIterator<Item = (String, MotifStats)> {
    let unique_samples = unique_sequences(sample_sequences);
    let unique_refs = unique_sequences(reference_sequences);
    let sample_members = sample_motif_members(&unique_samples);
    let reference_counts = reference_motif_counts(&unique_refs, &sample_members);

    let mut sample_counts: Vec<(&str, Vec<&str>)> = sample_members.into_iter().collect();
    sample_counts.par_sort_unstable_by_key(|&(motif, _)| motif);

    let total_samples = sample_sequences.len();
    let total_refs = reference_sequences.len();
    let log_fact_table = LogFactorialTable::new(total_refs + total_samples);
    sample_counts.into_par_iter().map(move |(motif, members)| {
        let reference_count = reference_counts.get(motif).copied().unwrap_or(0);
        let members: BTreeSet<String> = members.into_iter().map(str::to_owned).collect();
        let motif = motif.to_owned();
        let observed = members.len();
        // let fisher = hypergeometric_tail(
        //     total_samples + total_refs,
        //     observed + reference_count,
        //     total_samples,
        //     observed,
        // );
        let cont_table = ContingencyTable2x2::from_reported_counts(
            observed,
            total_samples,
            reference_count,
            total_refs,
        )
        .expect("Failed to create contingency table for global cluster");
        let (fisher, _) = log_fact_table.fishers_exact(&cont_table);
        let fold_change = (observed as f64 / (reference_count as f64 + 0.01))
            / total_samples as f64
            * total_refs as f64;
        (
            motif,
            MotifStats {
                members,
                fisher,
                fold_change,
            },
        )
    })
}

fn get_p_value_and_fc_sim(
    motif: &String,
    stat: &MotifStats,
    sim_depth: usize,
    sample_sequences: &[String],
    reference_sequences: &[String],
) -> (f64, f64) {
    let progress = progress_bar_with_message(
        reference_sequences.len(),
        "Preparing GLIPH1 reference motifs",
        true,
        true,
    )
    .expect("Failed to create progress bar for GLIPH1 reference motifs");
    let reference_motifs: Vec<BTreeSet<String>> = reference_sequences
        .par_iter()
        .map(|x| motifs(x))
        .progress_with(progress)
        .collect();
    // progress.finish_with_message("GLIPH1 reference motifs ready");

    let observed = stat.members.len();
    let simulations = sim_depth.max(1);
    let mut rng = crate::tcr_cluster::seeded_random(motif);
    let sample_size = sample_sequences.len();
    let mut total = 0_usize;
    let mut extreme = 0_usize;
    for _ in 0..simulations {
        let mut sampled_count = 0_usize;
        for _ in 0..sample_size {
            let index =
                crate::tcr_cluster::random_uniform_index(&mut rng, reference_sequences.len());
            if reference_motifs[index].contains(motif) {
                sampled_count += 1;
            }
        }
        total += sampled_count;
        if sampled_count >= observed {
            extreme += 1;
        }
    }
    let average = total as f64 / simulations as f64;
    let ratio = if average > 0.0 {
        observed as f64 / average
    } else {
        1.0 / (simulations as f64 * sample_size.max(1) as f64)
    };
    (
        (extreme as f64 / simulations as f64).max(1.0 / simulations as f64),
        ratio,
    )
}

#[inline]
pub fn get_p_value_and_fc(
    method: Algorithm,
    motif: &String,
    stat: &MotifStats,
    sim_depth: usize,
    sample_sequences: &[String],
    reference_sequences: &[String],
) -> (f64, f64) {
    match method {
        Algorithm::Gliph1 => get_p_value_and_fc_sim(
            motif,
            stat,
            sim_depth,
            sample_sequences,
            reference_sequences,
        ),
        Algorithm::Gliph2 => (stat.fisher, stat.fold_change),
    }
}

pub fn motif_region(sequence: &str) -> &str {
    if sequence.len() > 6 {
        &sequence[3..sequence.len() - 3]
    } else {
        ""
    }
}

/// True if `pattern` (with `%` matching any single residue) occurs in the sequence's motif region.
fn pattern_matches(region: &[u8], pattern: &[u8]) -> bool {
    pattern.len() <= region.len()
        && region
            .windows(pattern.len())
            .any(|w| w.iter().zip(pattern).all(|(&r, &p)| p == b'%' || p == r))
}

/// One-sided Fisher p-value for an arbitrary-length pattern (`%` = wildcard),
/// comparing unique sample sequences against the reference.
pub fn pattern_fisher(
    pattern: &str,
    sample_sequences: &[String],
    reference_sequences: &[String],
    table: &LogFactorialTable,
) -> f64 {
    let pat = pattern.as_bytes();
    let count = |seqs: &[String]| {
        seqs.par_iter()
            .filter(|s| pattern_matches(motif_region(s).as_bytes(), pat))
            .count()
    };
    let observed = count(sample_sequences);
    let reference = count(reference_sequences);
    match ContingencyTable2x2::from_reported_counts(
        observed,
        sample_sequences.len(),
        reference,
        reference_sequences.len(),
    ) {
        Ok(t) => table.fishers_exact(&t).0,
        Err(_) => 1.0,
    }
}

pub fn motifs(sequence: &str) -> BTreeSet<String> {
    let region = motif_region(sequence).as_bytes();
    let mut found = BTreeSet::new();
    for length in MOTIF_LENGTHS {
        if region.len() < length {
            continue;
        }
        for start in 0..=region.len() - length {
            found.insert(String::from_utf8_lossy(&region[start..start + length]).to_string());
        }
    }
    found
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    #[test]
    fn test_once_lock_single_initialization() {
        // First call initializes the table
        let ptr1 = LogFactorialTable::global(100_000) as *const _;
        // Second call with different max_n reuses the already initialized table
        let ptr2 = LogFactorialTable::global(500_000) as *const _;

        assert_eq!(
            ptr1, ptr2,
            "Pointers must be identical (same memory instance)"
        );
    }

    #[test]
    fn trims_structural_boundaries_before_motif_search() {
        assert_eq!(motif_region("CASSIRSSYEQYF"), "SIRSSYE");
        assert!(motifs("CASSIRSSYEQYF").contains("RS"));
        assert!(!motifs("CASSIRSSYEQYF").contains("CA"));
    }

    #[test]
    fn test_table_construction() {
        // Query: 15 hits out of 100 unique CDR3s
        // Reference: 20 hits out of 1000 unique CDR3s
        let table = ContingencyTable2x2::from_reported_counts(15, 100, 20, 1000).unwrap();

        assert_eq!(table.a, 15);
        assert_eq!(table.b, 85);
        assert_eq!(table.c, 20);
        assert_eq!(table.d, 980);

        // Odds ratio: (15 * 980) / (85 * 20) = 14700 / 1700 ≈ 8.647
        let or = table.odds_ratio();
        assert!((or - 8.64705).abs() < 1e-4);
    }

    #[test]
    fn test_fishers_exact_classic_case() {
        // Standard Fisher tea-tasting or clinical example:
        // Table: [[12, 5], [3, 10]] -> N_total = 30
        let table = ContingencyTable2x2 {
            a: 12,
            b: 5,
            c: 3,
            d: 10,
        };

        let lookup = LogFactorialTable::new(50);
        let (p_one, p_two) = lookup.fishers_exact(&table);

        // Ground truth from R's fisher.test(matrix(c(12, 3, 5, 10), nrow=2)):
        // One-sided alternative="less": p ≈ 0.01266
        // Two-sided: p ≈ 0.02533
        assert!((p_one - 0.01266).abs() < 1e-5);
        assert!((p_two - 0.02533).abs() < 1e-5);
    }

    #[test]
    fn test_cdr3_motif_enrichment() {
        // Query repertoire: 200 CDR3s, 30 contain motif "YNEQF"
        // Ref repertoire: 10,000 CDR3s, 250 contain motif "YNEQF"
        let n_query = 200;
        let n_ref = 10000;
        let a = 30;
        let c = 250;

        let table = ContingencyTable2x2::from_reported_counts(a, n_query, c, n_ref).unwrap();
        let lookup = LogFactorialTable::new(n_query + n_ref);

        let (p_enrich, p_two) = lookup.fishers_exact(&table);

        // Highly enriched: OR ≈ 7.0
        assert!(table.odds_ratio() > 6.5);
        assert!(p_enrich < 1e-10);
        assert!(p_two < 1e-10);
    }

    #[test]
    fn test_stirling_accuracy_and_out_of_bounds_recovery() {
        // Table only sized to 10,000
        let lookup = LogFactorialTable::new(10_000);

        // Access index 960,194 (as in the panic)
        let val_stirling = lookup.ln_fact(960_194);
        assert!(val_stirling > 0.0);

        // Compare against full table computation for a known boundary
        let full_lookup = LogFactorialTable::new(20_000);
        let exact_val = full_lookup.ln_fact(15_000);
        let fallback_val = lookup.ln_fact(15_000);

        // Relative error must be < 1e-12
        let rel_error = ((exact_val - fallback_val) / exact_val).abs();
        assert!(rel_error < 1e-12, "Rel error {} too high", rel_error);
    }

    #[test]
    fn test_large_repertoire_without_panic() {
        // Query set: 383,632 CDR3s, 50 motif hits
        // Reference set: 576,562 CDR3s, 20 motif hits
        // Total N = 960,194 (matching the error dimensions)
        let n_query = 383_632;
        let n_ref = 576_562;
        let a = 50;
        let c = 20;

        let table = ContingencyTable2x2::from_reported_counts(a, n_query, c, n_ref).unwrap();

        // Even if lookup table was under-allocated to 500,000, it must not panic:
        let lookup = LogFactorialTable::new(500_000);
        let (p_enrich, p_two) = lookup.fishers_exact(&table);

        assert!((0.0..=1.0).contains(&p_enrich));
        assert!((0.0..=1.0).contains(&p_two));
        assert!(table.odds_ratio() > 1.0);
    }

    #[test]
    fn borrowed_motifs_match_owned_motifs() {
        let mut scratch = Vec::new();
        for seq in ["CASSIRSSYEQYF", "CASSLGQAYEQYF", "CASSSSSSF", "CAS", "CASSF"] {
            distinct_motifs_into(seq, &mut scratch);
            let borrowed: BTreeSet<String> = scratch.iter().map(|m| m.to_string()).collect();
            assert_eq!(borrowed, motifs(seq), "motif sets differ for {seq}");
            assert_eq!(borrowed.len(), scratch.len(), "scratch not deduplicated for {seq}");
        }
    }

    #[test]
    fn motif_statistics_counts_unique_sequences() {
        let to_vec = |v: &[&str]| v.iter().map(|s| s.to_string()).collect::<Vec<_>>();
        // Duplicate sample sequence must only count once as a member.
        let sample = to_vec(&["CASSIRSSYEQYF", "CASSIRSSYEQYF", "CASSLGQAYEQYF"]);
        let reference = to_vec(&["CASSIRSTYEQYF", "CASSPGQAYEQYF", "CASSPGQAYEQYF"]);
        let stats: Vec<(String, MotifStats)> =
            motif_statistics_iter(&sample, &reference).collect();

        let motifs_in_order: Vec<&String> = stats.iter().map(|(m, _)| m).collect();
        let mut sorted = motifs_in_order.clone();
        sorted.sort();
        assert_eq!(motifs_in_order, sorted, "motifs must be yielded in sorted order");

        let get = |m: &str| &stats.iter().find(|(k, _)| k == m).unwrap().1;
        assert_eq!(get("SIRS").members.len(), 1);
        assert_eq!(get("SS").members.len(), 1); // only CASSIRSSYEQYF has SS in its region
        assert_eq!(get("GQA").members.len(), 1);
        assert_eq!(get("YE").members.len(), 2);
        // GQA appears in one unique reference sequence (the duplicate counts once):
        // fold_change = (1 / (1 + 0.01)) / 3 * 3
        assert!((get("GQA").fold_change - 1.0 / 1.01).abs() < 1e-12);
    }

    #[test]
    fn test_invalid_counts() {
        // Hits cannot exceed total sequences
        assert!(ContingencyTable2x2::from_reported_counts(105, 100, 10, 50).is_err());
    }
}
