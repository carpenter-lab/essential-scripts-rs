use clap::ValueEnum;
use serde::Serialize;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::fmt::{Display, Formatter};
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct GlobalFilterConfig {
    pub fisher_threshold: f64,
    pub min_members: usize,
    pub max_reference_fraction: f64,
}

impl GlobalFilterConfig {
    pub fn new(fisher_threshold: f64, min_members: usize, max_reference_fraction: f64) -> Self {
        Self {
            fisher_threshold,
            min_members,
            max_reference_fraction,
        }
    }
}

impl Default for GlobalFilterConfig {
    fn default() -> Self {
        Self::new(0.01, 2, 0.01)
    }
}

pub const HLA_COLUMNS: [&str; 8] = [
    "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
];

pub const OUTPUT_COLUMNS: [&str; 28] = [
    "index",
    "pattern",
    "number_subject",
    "number_unique_cdr3",
    "TcRb",
    "Fisher_score",
    "Freq",
    "subject",
    "condition",
    "final_score",
    "hla_score",
    "vb_score",
    "expansion_score",
    "length_score",
    "cluster_size_score",
    "type",
    "ulTcRb",
    "V",
    "J",
    "TcRa",
    "HLA-DPA1",
    "HLA-DPB1",
    "HLA-DQA1",
    "HLA-DQB1",
    "HLA-DRB1",
    "HLA-DRB3",
    "HLA-DRB4",
    "HLA-DRB5",
];

#[derive(Clone, Debug, Default)]
pub struct ChainRecord {
    pub subject: String,
    pub condition: String,
    pub beta: String,
    pub alpha: String,
    pub v_gene: String,
    pub j_gene: String,
    pub counts: Option<f64>,
    pub hla: [String; 8],
}

#[derive(Clone, Debug, Default)]
pub struct HlaRecord {
    pub alleles: [String; 8],
}

#[derive(Clone, Debug)]
pub struct ReferenceData {
    pub sequences: Vec<String>,
    pub v_gene_frequencies: Vec<(String, f64)>,
    pub cdr3_length_frequencies: Vec<(usize, f64)>,
    pub cluster_size_probabilities: Vec<Vec<f64>>,
    pub cluster_size_sample_sizes: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, ValueEnum, Default)]
pub enum References {
    HumanV1CD4,
    HumanV1CD8,
    HumanV1CD4CD8,
    HumanV2CD4,
    HumanV2CD8,
    #[default]
    HumanV2CD4CD8,
    MouseV1CD4,
    MouseV1CD8,
    MouseV1CD4CD8,
    MouseV2CD4,
    MouseV2CD8,
    MouseV2CD4CD8,
}

impl Display for References {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            References::HumanV1CD4 => write!(f, "human_v1.0_CD4"),
            References::HumanV1CD8 => write!(f, "human_v1.0_CD8"),
            References::HumanV1CD4CD8 => write!(f, "human_v1.0_CD48"),
            References::HumanV2CD4 => write!(f, "human_v2.0_CD4"),
            References::HumanV2CD8 => write!(f, "human_v2.0_CD8"),
            References::HumanV2CD4CD8 => write!(f, "human_v2.0_CD48"),
            References::MouseV1CD4 => write!(f, "mouse_v1.0_CD4"),
            References::MouseV1CD8 => write!(f, "mouse_v1.0_CD8"),
            References::MouseV1CD4CD8 => write!(f, "mouse_v1.0_CD48"),
            References::MouseV2CD4 => write!(f, "mouse_v2.0_CD4"),
            References::MouseV2CD8 => write!(f, "mouse_v2.0_CD8"),
            References::MouseV2CD4CD8 => write!(f, "mouse_v2.0_CD48"),
        }
    }
}

#[derive(Clone, Debug, PartialEq, PartialOrd, Ord, Eq, Copy)]
pub enum ClusterKind {
    Local,
    Global,
    Single,
}

impl Display for ClusterKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ClusterKind::Local => write!(f, "local"),
            ClusterKind::Global => write!(f, "global"),
            ClusterKind::Single => write!(f, "single"),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Cluster {
    pub pattern: String,
    pub kind: ClusterKind,
    pub fisher_score: f64,
    pub member_sequences: Vec<String>,
}

impl Cluster {
    pub fn new(
        pattern: String,
        kind: ClusterKind,
        fisher_score: f64,
        member_sequences: Vec<String>,
    ) -> Self {
        Self {
            pattern,
            kind,
            fisher_score,
            member_sequences: member_sequences.into_iter().collect(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct ClusterRecord<'a> {
    pub records: Vec<&'a ChainRecord>,
    pub subject_count: usize,
    pub cluster: &'a Cluster,
    pub member_sequences: Arc<HashSet<&'a str>>,
    pub cluster_subjects: Arc<HashSet<&'a str>>,
}

impl<'a> ClusterRecord<'a> {
    pub fn new(records: &'a [ChainRecord], cluster: &'a Cluster) -> Self {
        let member_set: HashSet<&str> = cluster
            .member_sequences
            .iter()
            .map(|s| s.as_str())
            .collect();
        let cluster_records: Vec<&ChainRecord> = records
            .iter()
            .filter(|record| member_set.contains(&record.beta.as_str()))
            .collect();

        let cluster_subjects: HashSet<&'a str> = cluster_records
            .iter()
            .map(|record| record.subject.as_str())
            .collect();
        let subject_count = cluster_subjects.len();
        Self {
            records: cluster_records,
            subject_count,
            cluster,
            member_sequences: Arc::from(member_set),
            cluster_subjects: Arc::from(cluster_subjects),
        }
    }
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
}

#[derive(Clone, Debug)]
pub struct SampleScore {
    pub final_score: f64,
    pub hla_score: f64,
    pub vb_score: f64,
    pub expansion_score: f64,
    pub length_score: f64,
    pub cluster_size_score: f64,
    pub hla_hit: Vec<(usize, String)>,
}

#[derive(Debug, PartialEq, Serialize)]
pub struct OutputRow {
    pub index: usize,
    pub pattern: String,
    pub number_subject: usize,
    pub number_unique_cdr3: usize,
    #[serde(rename = "TcRb")]
    pub tc_rb: String,
    #[serde(rename = "Fisher_score")]
    pub fisher_score: f64,
    #[serde(rename = "Freq")]
    pub frequency: usize,
    pub subject: String,
    pub condition: String,
    pub final_score: f64,
    pub hla_score: f64,
    pub vb_score: f64,
    pub expansion_score: f64,
    pub length_score: f64,
    pub cluster_size_score: f64,
    #[serde(rename = "type")]
    pub kind: String,
    #[serde(rename = "ulTcRb")]
    pub unique_tc_rb: String,
    #[serde(rename = "V")]
    pub v_gene: String,
    #[serde(rename = "J")]
    pub j_gene: String,
    #[serde(rename = "TcRa")]
    pub tc_ra: String,
    #[serde(rename = "HLA-DPA1")]
    pub hla_dpa1: String,
    #[serde(rename = "HLA-DPB1")]
    pub hla_dpb1: String,
    #[serde(rename = "HLA-DQA1")]
    pub hla_dqa1: String,
    #[serde(rename = "HLA-DQB1")]
    pub hla_dqb1: String,
    #[serde(rename = "HLA-DRB1")]
    pub hla_drb1: String,
    #[serde(rename = "HLA-DRB3")]
    pub hla_drb3: String,
    #[serde(rename = "HLA-DRB4")]
    pub hla_drb4: String,
    #[serde(rename = "HLA-DRB5")]
    pub hla_drb5: String,
}

impl OutputRow {
    #[allow(clippy::too_many_arguments)]
    fn new(
        index: usize,
        pattern: String,
        number_subject: usize,
        number_unique_cdr3: usize,
        tc_rb: String,
        fisher_score: f64,
        frequency: usize,
        subject: String,
        condition: String,
        final_score: f64,
        hla_score: f64,
        vb_score: f64,
        expansion_score: f64,
        length_score: f64,
        cluster_size_score: f64,
        kind: String,
        unique_tc_rb: String,
        v_gene: String,
        j_gene: String,
        tc_ra: String,
        hla_dpa1: String,
        hla_dpb1: String,
        hla_dqa1: String,
        hla_dqb1: String,
        hla_drb1: String,
        hla_drb3: String,
        hla_drb4: String,
        hla_drb5: String,
    ) -> Self {
        OutputRow {
            index,
            pattern,
            number_subject,
            number_unique_cdr3,
            tc_rb,
            fisher_score,
            frequency,
            subject,
            condition,
            final_score,
            hla_score,
            vb_score,
            expansion_score,
            length_score,
            cluster_size_score,
            kind,
            unique_tc_rb,
            v_gene,
            j_gene,
            tc_ra,
            hla_dpa1,
            hla_dpb1,
            hla_dqa1,
            hla_dqb1,
            hla_drb1,
            hla_drb3,
            hla_drb4,
            hla_drb5,
        }
    }

    fn record_frequency(record: &ChainRecord) -> usize {
        record
            .counts
            .filter(|count| count.is_finite() && *count > 0.0)
            .map(|count| count.round().max(1.0) as usize)
            .unwrap_or(1)
    }

    pub fn new_from_records(
        index: usize,
        cluster_record: &ClusterRecord,
        sample_score: &SampleScore,
        record: &ChainRecord,
    ) -> Self {
        let mut hla = record.hla.clone();
        for (locus, allele) in &sample_score.hla_hit {
            hla[*locus] = crate::gliph::analysis::mark_hla(&hla[*locus], allele);
        }
        Self::new(
            index + 1,
            cluster_record.cluster.pattern.clone(),
            cluster_record.subject_count,
            cluster_record.member_sequences.len(),
            record.beta.clone(),
            cluster_record.cluster.fisher_score,
            Self::record_frequency(record),
            record.subject.clone(),
            record.condition.clone(),
            sample_score.final_score,
            sample_score.hla_score,
            sample_score.vb_score,
            sample_score.expansion_score,
            sample_score.length_score,
            sample_score.cluster_size_score,
            cluster_record.cluster.kind.to_string(),
            cluster_record.cluster.member_sequences.join(" "),
            record.v_gene.clone(),
            record.j_gene.clone(),
            record.alpha.clone(),
            hla[0].clone(),
            hla[1].clone(),
            hla[2].clone(),
            hla[3].clone(),
            hla[4].clone(),
            hla[5].clone(),
            hla[6].clone(),
            hla[7].clone(),
        )
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum, Default)]
pub enum Algorithm {
    Gliph1,
    #[default]
    Gliph2,
}

impl Display for Algorithm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Algorithm::Gliph1 => write!(f, "GLIPH1"),
            Algorithm::Gliph2 => write!(f, "GLIPH2"),
        }
    }
}

pub struct CumulativeSampler {
    cumulative_weights: Vec<f64>,
    total: f64,
}

impl CumulativeSampler {
    pub(crate) fn new(weights: &[f64]) -> Self {
        let mut cumulative_weights = Vec::with_capacity(weights.len());
        let mut running_sum = 0.0;

        for &weight in weights {
            if weight > 0.0 && weight.is_finite() {
                running_sum += weight;
            }
            cumulative_weights.push(running_sum);
        }

        Self {
            cumulative_weights,
            total: running_sum,
        }
    }

    #[inline]
    pub(crate) fn sample(&self, state: &mut u64) -> usize {
        let len = self.cumulative_weights.len();
        if len == 0 {
            return 0;
        }

        *state = state.wrapping_mul(6364136223846793005).wrapping_add(1);

        if self.total <= 0.0 {
            return (*state as usize) % len;
        }

        let target = ((*state >> 11) as f64 / ((1_u64 << 53) as f64)) * self.total;

        // Binary Search O(log W) instead of linear scan O(W)
        match self.cumulative_weights.binary_search_by(|&val| {
            val.partial_cmp(&target)
                .unwrap_or(std::cmp::Ordering::Greater)
        }) {
            Ok(idx) => idx.min(len - 1),
            Err(idx) => idx.min(len - 1),
        }
    }
}

pub trait RoundDecimal {
    fn round_dp(self, places: u32) -> f64;
}

impl RoundDecimal for f64 {
    fn round_dp(self, places: u32) -> f64 {
        let scale = 10_f64.powi(places as i32);
        (self * scale).round() / scale
    }
}

#[derive(Clone, Debug)]
pub struct MotifStats {
    pub members: BTreeSet<String>,
    pub fisher: f64,
    pub fold_change: f64,
}

#[derive(Clone, Debug)]
pub struct TaggedMember {
    pub sequence: String,
    pub residue: u8,
    pub position: usize,
}

#[derive(Clone, Debug)]
pub struct UnionFind {
    pub parent: Vec<usize>,
    pub rank: Vec<u8>,
}

impl UnionFind {
    pub fn new(size: usize) -> Self {
        Self {
            parent: (0..size).collect(),
            rank: vec![0; size],
        }
    }

    pub fn root(&mut self, node: usize) -> usize {
        if self.parent[node] != node {
            self.parent[node] = self.root(self.parent[node]);
        }
        self.parent[node]
    }

    pub fn join(&mut self, left: usize, right: usize) {
        let mut left = self.root(left);
        let mut right = self.root(right);
        if left == right {
            return;
        }
        if self.rank[left] < self.rank[right] {
            std::mem::swap(&mut left, &mut right);
        }
        self.parent[right] = left;
        if self.rank[left] == self.rank[right] {
            self.rank[left] += 1;
        }
    }
}

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
        let log_background: Vec<f64> = repertoire_counts
            .values()
            .map(|&count| count.ln_1p())
            .collect();

        Self {
            repertoire_counts,
            log_background,
            has_counts,
        }
    }
}
