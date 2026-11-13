use super::chain::is_valid_cdr3;
use clap::ValueEnum;
use std::fmt::{Display, Formatter};

#[derive(Clone, Debug)]
pub struct ReferenceData {
    pub sequences: Vec<String>,
    pub v_gene_frequencies: Vec<(String, f64)>,
    pub cdr3_length_frequencies: Vec<(usize, f64)>,
    pub cluster_size_probabilities: Vec<Vec<f64>>,
    pub cluster_size_sample_sizes: Vec<usize>,
}

impl ReferenceData {
    /// Reference sequences that pass the CDR3 validity filter.
    pub fn valid_sequences(&self) -> Vec<String> {
        self.sequences
            .iter()
            .filter(|sequence| is_valid_cdr3(sequence))
            .cloned()
            .collect()
    }

    /// Index of the precomputed cluster-size column whose sample size is closest to `n_records`.
    pub fn closest_size_column(&self, n_records: usize) -> Option<usize> {
        let total_records = n_records.max(1) as f64;
        self.cluster_size_sample_sizes
            .iter()
            .enumerate()
            .min_by_key(|&(_, sample_size)| {
                let ratio = *sample_size as f64 / total_records;
                ((1.0 - ratio).abs() * 1_000_000.0) as u64
            })
            .map(|(idx, _)| idx)
    }
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

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum, Default)]
pub enum Algorithm {
    Gliph1,
    #[default]
    Gliph2,
}

impl Display for Algorithm {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Algorithm::Gliph1 => write!(f, "GLIPH1"),
            Algorithm::Gliph2 => write!(f, "GLIPH2"),
        }
    }
}
