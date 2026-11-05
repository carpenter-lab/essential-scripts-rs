use crate::tcr_cluster::blosum::is_amino_acid_sequence;
use std::borrow::Cow;

pub const HLA_LOCI: usize = 8;

pub const HLA_COLUMNS: [&str; HLA_LOCI] = [
    "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
];

/// Allele used when a gene name is given without one (e.g. `TRBV19` -> `TRBV19*01`).
pub const DEFAULT_ALLELE: &str = "01";

/// A single input TCR chain (beta, optionally paired with alpha) from one subject.
///
/// Gene names are stored exactly as the user wrote them (so they are reported back
/// unchanged); use the `*_allele` methods when an allele-qualified name is needed.
#[derive(Clone, Debug, Default)]
pub struct ChainRecord {
    pub subject: String,
    pub condition: String,
    pub beta: String,
    pub alpha: Option<String>,
    /// Beta V gene (`TRBV` input column).
    pub v_gene: Option<String>,
    /// Beta J gene (`TRBJ` input column).
    pub j_gene: Option<String>,
    /// Alpha V gene (optional `TRAV` input column).
    pub va_gene: Option<String>,
    /// Alpha J gene (optional `TRAJ` input column).
    pub ja_gene: Option<String>,
    pub counts: Option<f64>,
    pub hla: [String; HLA_LOCI],
}

/// Appends `*01` to a gene name that has no allele; names with an allele are returned as is.
pub fn with_default_allele(gene: &str) -> Cow<'_, str> {
    if gene.contains('*') {
        Cow::Borrowed(gene)
    } else {
        Cow::Owned(format!("{gene}*{DEFAULT_ALLELE}"))
    }
}

impl ChainRecord {
    /// Beta V gene with an allele, e.g. `TRBV19*01`.
    pub fn vb_allele(&self) -> Option<Cow<'_, str>> {
        self.v_gene.as_deref().map(with_default_allele)
    }

    /// Beta J gene with an allele, e.g. `TRBJ2-1*01`.
    pub fn jb_allele(&self) -> Option<Cow<'_, str>> {
        self.j_gene.as_deref().map(with_default_allele)
    }

    /// Alpha V gene with an allele, e.g. `TRAV12-2*01`.
    pub fn va_allele(&self) -> Option<Cow<'_, str>> {
        self.va_gene.as_deref().map(with_default_allele)
    }

    /// Alpha J gene with an allele, e.g. `TRAJ33*01`.
    pub fn ja_allele(&self) -> Option<Cow<'_, str>> {
        self.ja_gene.as_deref().map(with_default_allele)
    }

    /// Whether this record has a non-empty CDR3α.
    pub fn has_alpha_cdr3(&self) -> bool {
        self.alpha.as_deref().is_some_and(|alpha| !alpha.is_empty())
    }

    /// Whether this record carries any alpha-chain gene annotation.
    pub fn has_alpha_genes(&self) -> bool {
        self.va_gene.is_some() || self.ja_gene.is_some()
    }

    /// Clonal frequency reported in the output (counts rounded, defaulting to 1).
    pub fn frequency(&self) -> usize {
        self.counts
            .filter(|count| count.is_finite() && *count > 0.0)
            .map(|count| count.round().max(1.0) as usize)
            .unwrap_or(1)
    }
}

#[derive(Clone, Debug, Default)]
pub struct HlaRecord {
    pub alleles: [String; HLA_LOCI],
}

/// Whether a CDR3 looks like a complete, canonical amino-acid CDR3 (`C...F`, at least 8 residues).
pub fn is_valid_cdr3(sequence: &str) -> bool {
    sequence.len() >= 8
        && sequence.starts_with('C')
        && sequence.ends_with('F')
        && is_amino_acid_sequence(sequence)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adds_default_allele_only_when_missing() {
        assert_eq!(with_default_allele("TRBV19"), "TRBV19*01");
        assert_eq!(with_default_allele("TRAV12-2*02"), "TRAV12-2*02");
        let record = ChainRecord {
            v_gene: Some("TRBV9".to_owned()),
            ja_gene: Some("TRAJ33".to_owned()),
            ..ChainRecord::default()
        };
        assert_eq!(record.vb_allele().as_deref(), Some("TRBV9*01"));
        assert_eq!(record.ja_allele().as_deref(), Some("TRAJ33*01"));
        assert_eq!(record.va_allele(), None);
        assert!(record.has_alpha_genes());
        assert_eq!(record.v_gene.as_deref(), Some("TRBV9"), "input is not rewritten");
    }
}
