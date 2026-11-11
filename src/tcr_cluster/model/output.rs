use super::{ChainRecord, ClusterMembers, ClusterScore};
use serde::Serialize;

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
    "VB",
    "JB",
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

/// Appended after [`OUTPUT_COLUMNS`] only when the input provided CDR3α sequences.
pub const ALPHA_PATTERN_COLUMN: &str = "pattern_alpha";

/// Appended after [`OUTPUT_COLUMNS`] (and [`ALPHA_PATTERN_COLUMN`]) only when the input
/// provided alpha V/J genes.
pub const ALPHA_GENE_COLUMNS: [&str; 2] = ["VA", "JA"];

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
    #[serde(rename = "VB")]
    pub v_gene: Option<String>,
    #[serde(rename = "JB")]
    pub j_gene: Option<String>,
    #[serde(rename = "TcRa")]
    pub tc_ra: Option<String>,
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
    /// Alpha V gene as given in the input; written as `VA` only when alpha genes are present.
    #[serde(skip)]
    pub va_gene: Option<String>,
    /// Alpha J gene as given in the input; written as `JA` only when alpha genes are present.
    #[serde(skip)]
    pub ja_gene: Option<String>,
    /// Motif of the cluster's CDR3α sequences; written as `pattern_alpha` only when the
    /// input has CDR3α sequences.
    #[serde(skip)]
    pub alpha_pattern: Option<String>,
}

impl OutputRow {
    pub fn has_alpha_genes(&self) -> bool {
        self.va_gene.is_some() || self.ja_gene.is_some()
    }

    pub fn has_alpha_cdr3(&self) -> bool {
        self.tc_ra.as_deref().is_some_and(|alpha| !alpha.is_empty())
    }

    /// Sets the cluster's CDR3α motif.
    pub fn with_alpha_pattern(mut self, pattern: Option<&str>) -> Self {
        self.alpha_pattern = pattern.map(str::to_owned);
        self
    }
}

impl OutputRow {
    /// Builds the row for one `record` that belongs to `members`' cluster.
    /// `index` is the zero-based cluster index; the output is one-based.
    pub fn new(
        index: usize,
        members: &ClusterMembers,
        score: &ClusterScore,
        record: &ChainRecord,
    ) -> Self {
        let mut hla = record.hla.clone();
        for (locus, allele) in &score.hla_hit {
            hla[*locus] = mark_hla(&hla[*locus], allele);
        }
        let [
            hla_dpa1,
            hla_dpb1,
            hla_dqa1,
            hla_dqb1,
            hla_drb1,
            hla_drb3,
            hla_drb4,
            hla_drb5,
        ] = hla;
        let cluster = members.cluster;
        OutputRow {
            index: index + 1,
            pattern: cluster.pattern.clone(),
            number_subject: members.subject_count,
            number_unique_cdr3: members.member_sequences.len(),
            tc_rb: record.beta.clone(),
            fisher_score: cluster.fisher_score,
            frequency: record.frequency(),
            subject: record.subject.clone(),
            condition: record.condition.clone(),
            final_score: score.final_score,
            hla_score: score.hla_score,
            vb_score: score.vb_score,
            expansion_score: score.expansion_score,
            length_score: score.length_score,
            cluster_size_score: score.cluster_size_score,
            kind: cluster.kind.to_string(),
            // Historically always emitted empty; kept for output compatibility.
            unique_tc_rb: String::new(),
            v_gene: record.v_gene.clone(),
            j_gene: record.j_gene.clone(),
            tc_ra: record.alpha.clone(),
            hla_dpa1,
            hla_dpb1,
            hla_dqa1,
            hla_dqb1,
            hla_drb1,
            hla_drb3,
            hla_drb4,
            hla_drb5,
            va_gene: record.va_gene.clone(),
            ja_gene: record.ja_gene.clone(),
            alpha_pattern: None,
        }
    }
}

/// Suffixes `allele` with `!` wherever it appears as a whole token in `value`.
fn mark_hla(value: &str, allele: &str) -> String {
    let mut out = String::with_capacity(value.len() + 1);
    let mut token = String::new();
    let flush = |token: &mut String, out: &mut String| {
        let trimmed = token.trim();
        if !trimmed.is_empty() && trimmed == allele && !trimmed.ends_with('!') {
            let end = token.trim_end().len();
            out.push_str(&token[..end]);
            out.push('!');
            out.push_str(&token[end..]);
            token.clear();
        } else {
            out.push_str(token);
        }
        token.clear();
    };
    for ch in value.chars() {
        if matches!(ch, ',' | ';' | '|' | '/') {
            flush(&mut token, &mut out);
            out.push(ch);
        } else {
            token.push(ch);
        }
    }
    flush(&mut token, &mut out);
    out
}

#[cfg(test)]
mod tests {
    use super::mark_hla;

    #[test]
    fn marks_identified_hla_allele() {
        assert_eq!(
            mark_hla("DRB1*15:01, DRB1*03:01", "DRB1*15:01"),
            "DRB1*15:01!, DRB1*03:01"
        );
        assert_eq!(
            mark_hla("DPA1*01:03/DPA1*02:01", "DPA1*02:01"),
            "DPA1*01:03/DPA1*02:01!"
        );
        assert_eq!(mark_hla("DPA1*01:03", "DPA1*01:04"), "DPA1*01:03");
        assert_eq!(mark_hla("01:03", "01:03"), "01:03!");
    }
}
