use crate::gliph::{Algorithm, MotifStats};
use crate::progress::progress_bar_with_message;
use indicatif::ParallelProgressIterator;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
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

fn get_motif_vec(seqs: &[String]) -> Vec<(&String, BTreeSet<String>)> {
    seqs.par_iter()
        .map(|sequence| {
            let motifs = motifs(sequence);
            (sequence, motifs)
        })
        .collect()
}

pub fn motif_statistics_par_iter(
    sample_sequences: &[String],
    reference_sequences: &[String],
) -> impl ParallelIterator<Item = (String, MotifStats)> {
    let mut sample_counts: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    let mut reference_counts: HashMap<String, BTreeSet<String>> = HashMap::new();
    let sample_motifs = get_motif_vec(sample_sequences);
    for (sequence, sequence_motifs) in &sample_motifs {
        for motif in sequence_motifs {
            sample_counts
                .entry(motif.clone())
                .or_default()
                .insert((*sequence).clone());
        }
    }
    let sample_motifs_set: HashSet<String> = sample_counts.keys().cloned().collect();
    let reference_motifs = get_motif_vec(reference_sequences);

    for (sequence, sequence_motifs) in reference_motifs {
        for motif in sequence_motifs {
            if sample_motifs_set.contains(&motif) {
                reference_counts
                    .entry(motif)
                    .or_default()
                    .insert(sequence.clone());
            }
        }
    }

    let total_samples = sample_sequences.len();
    let total_refs = reference_sequences.len();
    let log_fact_table = LogFactorialTable::new(total_refs + total_samples);
    sample_counts.into_par_iter().map(move |(motif, members)| {
        let reference_count = reference_counts.get(&motif).map_or(0, |s| s.len());
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
    )
    .expect("Failed to create progress bar for GLIPH1 reference motifs");
    let reference_motifs: Vec<BTreeSet<String>> = reference_sequences
        .par_iter()
        .progress_with(progress)
        .map(|x| motifs(x))
        .collect();
    // progress.finish_with_message("GLIPH1 reference motifs ready");

    let observed = stat.members.len();
    let simulations = sim_depth.max(1);
    let mut rng = crate::gliph::analysis::seeded_random(motif);
    let sample_size = sample_sequences.len();
    let mut total = 0_usize;
    let mut extreme = 0_usize;
    for _ in 0..simulations {
        let mut sampled_count = 0_usize;
        for _ in 0..sample_size {
            let index =
                crate::gliph::analysis::random_uniform_index(&mut rng, reference_sequences.len());
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
    fn test_invalid_counts() {
        // Hits cannot exceed total sequences
        assert!(ContingencyTable2x2::from_reported_counts(105, 100, 10, 50).is_err());
    }
}
