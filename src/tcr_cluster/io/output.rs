use crate::tcr_cluster::model::*;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

fn write_rows<W: Write>(writer: W, rows: &[OutputRow]) -> Result<()> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);
    // The pattern_alpha and VA/JA columns only appear when the input supplied CDR3α
    // sequences or alpha genes, so beta-only runs keep the original GLIPH column layout.
    let include_alpha_pattern = rows.iter().any(OutputRow::has_alpha_cdr3);
    let include_alpha_genes = rows.iter().any(OutputRow::has_alpha_genes);
    let header = OUTPUT_COLUMNS
        .iter()
        .chain(include_alpha_pattern.then_some(&ALPHA_PATTERN_COLUMN))
        .chain(
            include_alpha_genes
                .then_some(&ALPHA_GENE_COLUMNS)
                .into_iter()
                .flatten(),
        );
    writer
        .write_record(header)
        .context("Failed to write GLIPH TSV header")?;
    for row in rows {
        let (pattern, va, ja) = (&row.alpha_pattern, &row.va_gene, &row.ja_gene);
        let result = match (include_alpha_pattern, include_alpha_genes) {
            (false, false) => writer.serialize(row),
            (true, false) => writer.serialize((row, pattern)),
            (false, true) => writer.serialize((row, va, ja)),
            (true, true) => writer.serialize((row, pattern, va, ja)),
        };
        result.context("Failed to serialize a GLIPH result row")?;
    }
    writer.flush().context("Failed to flush GLIPH TSV output")
}

pub fn write_output(path: &Path, rows: &[OutputRow]) -> Result<()> {
    if path == Path::new("-") {
        let stdout = io::stdout();
        return write_rows(stdout.lock(), rows);
    }
    let file = File::create(path)
        .with_context(|| format!("Failed to create GLIPH output file {}", path.display()))?;
    write_rows(file, rows)
}

#[cfg(test)]
mod tests {
    use super::write_output;
    use crate::tcr_cluster::model::{
        ALPHA_GENE_COLUMNS, ALPHA_PATTERN_COLUMN, ChainRecord, OUTPUT_COLUMNS, OutputRow,
    };
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn writes_exact_tsv_column_order() {
        let temp = tempdir().unwrap();
        let path = temp.path().join("result.tsv");
        write_output(
            &path,
            &[OutputRow {
                index: 1,
                pattern: "RS".to_owned(),
                number_subject: 1,
                number_unique_cdr3: 2,
                tc_rb: "CASSIRSSYEQYF".to_owned(),
                fisher_score: 0.01,
                frequency: 2,
                subject: "S1".to_owned(),
                condition: "case".to_owned(),
                final_score: 0.1,
                hla_score: 0.2,
                vb_score: 0.3,
                expansion_score: 1.0,
                length_score: 0.4,
                cluster_size_score: 0.5,
                kind: "local".to_owned(),
                unique_tc_rb: "CASSIRSSYEQYF CASSPRSSYEQYF".to_owned(),
                v_gene: Some("TRBV9".to_owned()),
                j_gene: Some("TRBJ2-1".to_owned()),
                tc_ra: Some(String::new()),
                hla_dpa1: "01:03".to_owned(),
                hla_dpb1: String::new(),
                hla_dqa1: String::new(),
                hla_dqb1: String::new(),
                hla_drb1: String::new(),
                hla_drb3: String::new(),
                hla_drb4: String::new(),
                hla_drb5: String::new(),
                va_gene: None,
                ja_gene: None,
                alpha_pattern: None,
            }],
        )
        .unwrap();

        let text = fs::read_to_string(path).unwrap();
        assert_eq!(text.lines().next().unwrap(), OUTPUT_COLUMNS.join("\t"));
        assert_eq!(text.lines().count(), 2);
    }

    fn row_for(record: &ChainRecord) -> OutputRow {
        OutputRow {
            index: 1,
            pattern: "RS".to_owned(),
            number_subject: 1,
            number_unique_cdr3: 1,
            tc_rb: record.beta.clone(),
            fisher_score: 0.01,
            frequency: 1,
            subject: "S1".to_owned(),
            condition: "case".to_owned(),
            final_score: 0.1,
            hla_score: 0.2,
            vb_score: 0.3,
            expansion_score: 1.0,
            length_score: 0.4,
            cluster_size_score: 0.5,
            kind: "local".to_owned(),
            unique_tc_rb: String::new(),
            v_gene: record.v_gene.clone(),
            j_gene: record.j_gene.clone(),
            tc_ra: record.alpha.clone(),
            hla_dpa1: String::new(),
            hla_dpb1: String::new(),
            hla_dqa1: String::new(),
            hla_dqb1: String::new(),
            hla_drb1: String::new(),
            hla_drb3: String::new(),
            hla_drb4: String::new(),
            hla_drb5: String::new(),
            va_gene: record.va_gene.clone(),
            ja_gene: record.ja_gene.clone(),
            alpha_pattern: None,
        }
    }

    #[test]
    fn appends_alpha_gene_columns_as_given_when_present() {
        let temp = tempdir().unwrap();
        let path = temp.path().join("result.tsv");
        let with_alpha = ChainRecord {
            beta: "CASSIRSSYEQYF".to_owned(),
            v_gene: Some("TRBV9".to_owned()),
            va_gene: Some("TRAV12-2".to_owned()),
            ja_gene: Some("TRAJ33".to_owned()),
            ..ChainRecord::default()
        };
        let beta_only = ChainRecord {
            beta: "CASSPRSSYEQYF".to_owned(),
            ..ChainRecord::default()
        };
        write_output(&path, &[row_for(&with_alpha), row_for(&beta_only)]).unwrap();

        let text = fs::read_to_string(path).unwrap();
        let lines: Vec<Vec<&str>> = text.lines().map(|l| l.split('\t').collect()).collect();
        let expected_header: Vec<&str> = OUTPUT_COLUMNS
            .iter()
            .chain(&ALPHA_GENE_COLUMNS)
            .copied()
            .collect();
        assert_eq!(lines[0], expected_header);
        let n = expected_header.len();
        assert!(lines.iter().all(|l| l.len() == n));
        assert_eq!(
            lines[1][n - 2..],
            ["TRAV12-2", "TRAJ33"],
            "no allele added in output"
        );
        assert_eq!(
            lines[1][OUTPUT_COLUMNS.iter().position(|c| *c == "VB").unwrap()],
            "TRBV9"
        );
        assert_eq!(lines[2][n - 2..], ["", ""]);
        assert!(!text.contains("*01"));
    }

    #[test]
    fn appends_alpha_pattern_column_when_input_has_cdr3a() {
        let temp = tempdir().unwrap();
        let path = temp.path().join("result.tsv");
        let paired = ChainRecord {
            beta: "CASSIRSSYEQYF".to_owned(),
            alpha: Some("CAVRDSNYQLIW".to_owned()),
            va_gene: Some("TRAV12-2".to_owned()),
            ..ChainRecord::default()
        };
        let beta_only = ChainRecord {
            beta: "CASSPRSSYEQYF".to_owned(),
            ..ChainRecord::default()
        };
        write_output(
            &path,
            &[
                row_for(&paired).with_alpha_pattern(Some("RDSNYQ")),
                row_for(&beta_only),
            ],
        )
        .unwrap();

        let text = fs::read_to_string(path).unwrap();
        let lines: Vec<Vec<&str>> = text.lines().map(|l| l.split('\t').collect()).collect();
        let expected_header: Vec<&str> = OUTPUT_COLUMNS
            .iter()
            .chain(&[ALPHA_PATTERN_COLUMN])
            .chain(&ALPHA_GENE_COLUMNS)
            .copied()
            .collect();
        assert_eq!(lines[0], expected_header);
        let n = expected_header.len();
        assert!(lines.iter().all(|l| l.len() == n));
        assert_eq!(lines[1][n - 3..], ["RDSNYQ", "TRAV12-2", ""]);
        assert_eq!(lines[2][n - 3..], ["", "", ""]);
    }

    #[test]
    fn beta_only_output_keeps_original_columns() {
        let temp = tempdir().unwrap();
        let path = temp.path().join("result.tsv");
        let record = ChainRecord {
            beta: "CASSIRSSYEQYF".to_owned(),
            ..ChainRecord::default()
        };
        write_output(&path, &[row_for(&record)]).unwrap();
        let text = fs::read_to_string(path).unwrap();
        assert_eq!(text.lines().next().unwrap(), OUTPUT_COLUMNS.join("\t"));
        assert!(
            text.lines()
                .all(|l| l.split('\t').count() == OUTPUT_COLUMNS.len())
        );
    }
}
