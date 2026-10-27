use crate::gliph::RoundDecimal;
use crate::gliph::analysis::{random_uniform_index, seeded_random};
use crate::gliph::stats::{ContingencyTable2x2, LogFactorialTable};
use crate::gliph::types::*;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;

pub fn randomized_frequency_score<T, K>(
    observed: &[T],
    reference_frequencies: &[(K, f64)],
    simulations: usize,
    seed: &str,
    gliph1: bool,
) -> f64
where
    T: Hash + Eq,
{
    if observed.is_empty() || reference_frequencies.is_empty() {
        return 1.0;
    }

    let n = observed.len();
    let n_f64 = n as f64;

    // K requires no bounds because we only read the float weight
    let weights: Vec<f64> = reference_frequencies
        .iter()
        .map(|(_, frequency)| *frequency)
        .collect();

    // HashMap borrows &T directly — works for &usize, &&str, &String with 0 allocations
    let mut observed_counts: HashMap<&T, usize> = HashMap::with_capacity(n);
    for value in observed {
        *observed_counts.entry(value).or_default() += 1;
    }

    let observed_score = observed_counts
        .values()
        .map(|&count| count as f64 / n_f64)
        .product::<f64>()
        .round_dp(3);

    let mut rng = seeded_random(seed);
    let mut at_least_as_extreme = 0;

    // Reuse the buffer across all simulations
    let mut counts = vec![0_usize; weights.len()];
    let sampler = CumulativeSampler::new(&weights);
    for _ in 0..simulations {
        counts.fill(0);
        for _ in 0..observed.len() {
            counts[sampler.sample(&mut rng)] += 1;
        }

        let mut raw_product = 1_u64;
        let mut non_zero_categories = 0i32;

        for &count in &counts {
            if count > 0 {
                raw_product *= count as u64;
                non_zero_categories += 1;
            }
        }

        let simulated_score = (raw_product as f64 / n_f64.powi(non_zero_categories)).round_dp(3);
        let extreme = if gliph1 {
            simulated_score >= observed_score
        } else {
            simulated_score > observed_score
        };
        if extreme {
            at_least_as_extreme += 1;
        }
    }

    let total_sims = simulations.max(1) as f64;
    (at_least_as_extreme as f64 / total_sims).max(1.0 / total_sims)
}

pub fn expansion_score(
    member_sequences: &HashSet<&str>,
    exp_ctx: &ExpansionContext,
    simulations: usize,
    seed: &str,
) -> f64 {
    if !exp_ctx.has_counts || exp_ctx.log_background.is_empty() {
        return 1.0;
    }

    let bg_len = exp_ctx.log_background.len();

    // 1. Direct O(K) lookup of cluster members (No full repertoire scanning)
    let mut observed_log_sum = 0.0;
    let mut cluster_size = 0;

    for &seq in member_sequences {
        if let Some(&count) = exp_ctx.repertoire_counts.get(seq) {
            observed_log_sum += count.ln_1p();
            cluster_size += 1;
        }
    }

    if cluster_size == 0 || bg_len < cluster_size {
        return 1.0;
    }

    let simulations = simulations.max(1);
    let mut rng = seeded_random(seed);
    let mut at_least_as_extreme = 0;

    // 2. CRITICAL: Allocate the scratch buffer ONCE outside the loop
    let mut sampled: Vec<usize> = Vec::with_capacity(cluster_size);

    for _ in 0..simulations {
        sampled.clear();

        // Fast sampling without replacement:
        // For cluster sizes typical in TCR/BCR data (<= 64), a contiguous slice check
        // runs in L1 cache and is orders of magnitude faster than a HashSet.
        while sampled.len() < cluster_size {
            let idx = random_uniform_index(&mut rng, bg_len);
            if !sampled.contains(&idx) {
                sampled.push(idx);
            }
        }

        // 3. Simple array lookups of pre-calculated logarithms (0 ln_1p calls!)
        let simulated_log_sum: f64 = sampled.iter().map(|&idx| exp_ctx.log_background[idx]).sum();

        // 4. Compare raw sums directly (avoids float division)
        if simulated_log_sum >= observed_log_sum {
            at_least_as_extreme += 1;
        }
    }

    let sims_f = simulations as f64;
    (at_least_as_extreme as f64 / sims_f).max(1.0 / sims_f)
}

pub fn sample_scores(
    cluster_record: &ClusterRecord,
    reference: &ReferenceData,
    algorithm: Algorithm,
    sim_depth: usize,
    precomputed_size_col: Option<usize>,
    hla_context: &HlaContext,
    exp_ctx: &ExpansionContext,
) -> SampleScore {
    // Cluster size score using precalculated column index
    let size_index = cluster_record
        .cluster
        .member_sequences
        .len()
        .saturating_sub(1);

    let cluster_size = precomputed_size_col
        .and_then(|col| {
            reference
                .cluster_size_probabilities
                .get(size_index)?
                .get(col)
        })
        .copied()
        .unwrap_or(1.0);

    // CDR3 Lengths: Borrow directly from member_set (no BTreeSet, no heap strings)
    // If randomized_frequency_score is generic over T: AsRef<str> or Eq:
    let length_values: Vec<String> = cluster_record
        .member_sequences
        .iter()
        .map(|seq| seq.len().to_string())
        .collect();

    let length_score = randomized_frequency_score(
        &length_values,
        &reference.cdr3_length_frequencies, // Pre-convert once in ReferenceData!
        sim_depth,
        &format!("length:{}", cluster_record.cluster.pattern),
        algorithm == Algorithm::Gliph1,
    );

    // V-gene values: Filter and collect references or clone only when strictly required
    let v_values: Vec<String> = cluster_record
        .records
        .iter()
        .map(|record| record.v_gene.clone())
        .filter(|gene| !gene.is_empty())
        .collect();

    let v_score = if v_values.is_empty() {
        1.0
    } else {
        randomized_frequency_score(
            &v_values,
            &reference.v_gene_frequencies,
            sim_depth,
            &format!("v:{}", cluster_record.cluster.pattern),
            algorithm == Algorithm::Gliph1,
        )
    };

    // Expansion and HLA scoring
    let expansion_score = expansion_score(
        &cluster_record.member_sequences,
        exp_ctx,
        sim_depth,
        &format!("expansion:{}", cluster_record.cluster.pattern),
    );

    // HLA scoring
    let (hla, hla_hit) = hla_score(cluster_record, hla_context);

    let scale = if algorithm == Algorithm::Gliph1 {
        0.064
    } else {
        1.0
    };

    let final_score = cluster_size * length_score * v_score * expansion_score * hla * scale;

    SampleScore {
        final_score,
        hla_score: hla,
        vb_score: v_score,
        expansion_score,
        length_score,
        cluster_size_score: cluster_size,
        hla_hit,
    }
}

fn hla_score(records: &ClusterRecord, ctx: &HlaContext) -> (f64, Vec<(usize, String)>) {
    if ctx.population_size == 0 {
        return (1.0, Vec::new());
    }
    let log_fact_table = LogFactorialTable::new(ctx.population_size);

    let member_sequences = &records.member_sequences;
    let records = &records.records;
    let mut cluster_subjects = HashSet::with_capacity(member_sequences.len());
    for record in records {
        if member_sequences.contains(record.beta.as_str()) {
            cluster_subjects.insert(record.subject.as_str());
        }
    }

    let draws = cluster_subjects.len();
    if draws == 0 {
        return (1.0, Vec::new());
    }

    let mut best = 1.0_f64;
    let mut hits = Vec::new();

    // Iterate over pre-computed unique alleles per locus
    for locus_idx in 0..8 {
        for allele in &ctx.locus_unique_alleles[locus_idx] {
            let mut total_with = 0;
            let mut cluster_with = 0;

            for (subject, loci) in &ctx.subject_loci {
                if loci[locus_idx].contains(allele) {
                    total_with += 1;
                    if cluster_subjects.contains(subject.as_str()) {
                        cluster_with += 1;
                    }
                }
            }

            if cluster_with > 1 {
                let cont_table = ContingencyTable2x2::from_reported_counts(
                    cluster_with,
                    draws,
                    total_with,
                    ctx.population_size,
                )
                .expect("Failed to create contingency table HLA");
                // let p_value = stats::hypergeometric_tail(
                //     ctx.population_size,
                //     total_with,
                //     draws,
                //     cluster_with,
                // );
                let (p_value, _) = log_fact_table.fishers_exact(&cont_table);

                best = best.min(p_value);
                if p_value <= ctx.threshold {
                    hits.push((locus_idx, allele.clone()));
                }
            }
        }
    }

    (best, hits)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(subject: &str, beta: &str, dpa1: &str) -> ChainRecord {
        let mut hla: [String; 8] = Default::default();
        hla[0] = dpa1.to_owned();
        ChainRecord {
            subject: subject.to_owned(),
            beta: beta.to_owned(),
            hla,
            ..Default::default()
        }
    }

    #[test]
    fn hla_score_reports_full_resolution_hits() {
        let mut records = vec![
            record("S1", "CASSAAF", "DPA1*01:03/DPA1*02:01"),
            record("S2", "CASSAAF", "DPA1*01:03"),
            record("S3", "CASSAAF", "DPA1*01:03/DPA1*01:04"),
        ];
        for i in 4..=20 {
            records.push(record(&format!("S{i}"), "CASSXXF", "DPA1*02:01"));
        }
        let cluster = Cluster {
            pattern: "SSAA".to_owned(),
            kind: ClusterKind::Local,
            fisher_score: 0.01,
            member_sequences: vec!["CASSAAF".to_owned()],
        };
        let ctx = HlaContext::build(&records, 0.05);
        let cluster_record = ClusterRecord::new(&records, &cluster);

        let (score, hits) = hla_score(&cluster_record, &ctx);

        assert!(score <= 0.05, "score {score}");
        assert!(hits.contains(&(0, "DPA1*01:03".to_owned())), "{hits:?}");
    }
}
