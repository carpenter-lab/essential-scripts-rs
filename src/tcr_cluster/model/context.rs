use super::{Algorithm, ChainRecord, ReferenceData};
use std::collections::{HashMap, HashSet};

/// Pre-parsed HLA lookup data built once per run
pub struct HlaContext {
    /// Maps each subject -> 8 loci, each locus holding a set of canonical alleles
    pub subject_loci: HashMap<String, [HashSet<String>; 8]>,
    /// Pre-extracted unique alleles per locus across the entire population
    pub locus_unique_alleles: [Vec<String>; 8],
    /// Total population count (number of distinct subjects)
    pub population_size: usize,
    pub threshold: f64,
}

impl HlaContext {
    pub fn build(records: &[ChainRecord], threshold: f64) -> Self {
        let mut subject_loci: HashMap<String, [HashSet<String>; 8]> = HashMap::new();

        // 1. Pre-parse every subject's alleles once
        for record in records {
            let entry = subject_loci.entry(record.subject.clone()).or_default();

            for (locus_idx, allele_str) in record.hla.iter().enumerate() {
                if locus_idx < 8 {
                    for allele in canonical_hla_borrowed(allele_str) {
                        entry[locus_idx].insert(allele.to_string());
                    }
                }
            }
        }

        // 2. Pre-extract all unique alleles per locus across the entire population
        let mut locus_unique_alleles: [Vec<String>; 8] = Default::default();
        for locus_idx in 0..8 {
            let mut unique = HashSet::new();
            for loci in subject_loci.values() {
                for allele in &loci[locus_idx] {
                    unique.insert(allele.clone());
                }
            }
            locus_unique_alleles[locus_idx] = unique.into_iter().collect();
        }

        let population_size = subject_loci.len();

        Self {
            subject_loci,
            locus_unique_alleles,
            population_size,
            threshold,
        }
    }
}

/// Zero-allocation helper to split a cell into full-resolution alleles
#[inline]
fn canonical_hla_borrowed(value: &str) -> impl Iterator<Item = &str> {
    value
        .split([',', ';', '|', '/'])
        .map(str::trim)
        .filter(|allele| !allele.is_empty())
}

pub struct ExpansionContext {
    /// Maps beta sequence -> total count across the repertoire
    pub repertoire_counts: HashMap<String, f64>,
    /// Pre-computed ln_1p() values for each sequence in the background
    pub log_background: Vec<f64>,
    /// Whether any record had non-empty counts
    pub has_counts: bool,
}

impl ExpansionContext {
    pub fn build(records: &[ChainRecord]) -> Self {
        let has_counts = records.iter().any(|r| r.counts.is_some());
        if !has_counts {
            return Self {
                repertoire_counts: HashMap::new(),
                log_background: Vec::new(),
                has_counts: false,
            };
        }

        let mut repertoire_counts: HashMap<String, f64> = HashMap::with_capacity(records.len());
        for record in records {
            let count = record.counts.unwrap_or(1.0).max(1.0);
            *repertoire_counts.entry(record.beta.clone()).or_default() += count;
        }

        // Precompute log(1 + count) once for all background items
        // Sorted so that simulated draws do not depend on HashMap iteration order.
        let mut log_background: Vec<f64> = repertoire_counts
            .values()
            .map(|&count| count.ln_1p())
            .collect();
        log_background.sort_by(f64::total_cmp);

        Self {
            repertoire_counts,
            log_background,
            has_counts,
        }
    }
}

/// Everything needed to score any cluster, regardless of which algorithm discovered it.
pub struct ScoringContext<'a> {
    pub reference: &'a ReferenceData,
    pub algorithm: Algorithm,
    pub sim_depth: usize,
    pub hla: &'a HlaContext,
    pub expansion: &'a ExpansionContext,
    pub size_column: Option<usize>,
}

impl<'a> ScoringContext<'a> {
    pub fn new(
        reference: &'a ReferenceData,
        algorithm: Algorithm,
        sim_depth: usize,
        hla: &'a HlaContext,
        expansion: &'a ExpansionContext,
        n_records: usize,
    ) -> Self {
        Self {
            reference,
            algorithm,
            sim_depth,
            hla,
            expansion,
            size_column: reference.closest_size_column(n_records),
        }
    }
}
