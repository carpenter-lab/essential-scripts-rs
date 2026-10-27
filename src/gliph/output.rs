use crate::gliph::types::*;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

fn write_rows<W: Write>(writer: W, rows: &[OutputRow]) -> Result<()> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);
    writer
        .write_record(OUTPUT_COLUMNS)
        .context("Failed to write GLIPH TSV header")?;
    for row in rows {
        writer
            .serialize(row)
            .context("Failed to serialize a GLIPH result row")?;
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
    use crate::gliph::types::{OUTPUT_COLUMNS, OutputRow};
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
                v_gene: "TRBV9".to_owned(),
                j_gene: "TRBJ2-1".to_owned(),
                tc_ra: String::new(),
                hla_dpa1: "01:03".to_owned(),
                hla_dpb1: String::new(),
                hla_dqa1: String::new(),
                hla_dqb1: String::new(),
                hla_drb1: String::new(),
                hla_drb3: String::new(),
                hla_drb4: String::new(),
                hla_drb5: String::new(),
            }],
        )
        .unwrap();

        let text = fs::read_to_string(path).unwrap();
        assert_eq!(text.lines().next().unwrap(), OUTPUT_COLUMNS.join("\t"));
        assert_eq!(text.lines().count(), 2);
    }
}
