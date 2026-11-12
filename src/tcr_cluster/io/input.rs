use crate::tcr_cluster::model::*;
use anyhow::{Context, Result, bail};
use csv::{ReaderBuilder, StringRecord};
use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

fn header_index(headers: &StringRecord, name: &str) -> Option<usize> {
    headers.iter().position(|header| header == name)
}

fn required_index(headers: &StringRecord, name: &str, path: &Path) -> Result<usize> {
    header_index(headers, name)
        .with_context(|| format!("Missing required column '{name}' in {}", path.display()))
}

fn cell(record: &StringRecord, index: Option<usize>) -> Option<&str> {
    index
        .and_then(|i| record.get(i))
        .map(str::trim)
        .filter(|&s| !s.is_empty() && !s.eq_ignore_ascii_case("NA"))
}

fn read_tsv(path: &Path) -> Result<csv::Reader<File>> {
    ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(false)
        .from_path(path)
        .with_context(|| format!("Failed to read TSV file {}", path.display()))
}

fn read_hla(path: &Path) -> Result<HashMap<String, HlaRecord>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(true)
        .from_path(path)
        .with_context(|| format!("Failed to read HLA TSV file {}", path.display()))?;
    let mut by_subject = HashMap::new();

    for (row_index, row) in reader.records().enumerate() {
        let row = row.with_context(|| {
            format!(
                "Failed to parse HLA row {} in {}",
                row_index + 2,
                path.display()
            )
        })?;
        let subject = match cell(&row, Some(0)) {
            Some(s) => s,
            None => bail!(
                "Empty subject on HLA row {} in {}",
                row_index + 2,
                path.display()
            ),
        };
        if row.len() > 13 {
            bail!(
                "HLA row {} in {} has {} columns, expected at most 13 (subject plus 12 allele columns)",
                row_index + 2,
                path.display(),
                row.len()
            );
        }
        let alleles = group_alleles(&row);
        if by_subject
            .insert(subject.to_string(), HlaRecord { alleles })
            .is_some()
        {
            bail!(
                "Duplicate HLA record for subject '{subject}' in {}",
                path.display()
            );
        }
    }

    if by_subject.is_empty() {
        bail!("No HLA records found in {}", path.display());
    }
    Ok(by_subject)
}

/// Collapses the 12 L/R allele columns into the 8 output loci, joining distinct alleles with '/'.
fn group_alleles(row: &StringRecord) -> [String; 8] {
    let mut groups: [Vec<&str>; 8] = Default::default();
    for column in 1..13 {
        let value = cell(row, Some(column));
        if value.is_none() {
            continue;
        }
        if let Some(v) = value {
            let locus = match column {
                1 | 2 => 0,
                3 | 4 => 1,
                5 | 6 => 2,
                7 | 8 => 3,
                9 | 10 => 4,
                _ if v.starts_with("DRB3") => 5,
                _ if v.starts_with("DRB4") => 6,
                _ if v.starts_with("DRB5") => 7,
                _ => continue,
            };
            if !groups[locus].contains(&v) {
                groups[locus].push(v);
            }
        }
    }
    groups.map(|alleles| alleles.join("/"))
}

pub fn read_chain_inputs(chains_path: &Path, hla_path: Option<&Path>) -> Result<Vec<ChainRecord>> {
    let hla_by_subject = hla_path.map(read_hla).transpose()?;
    let mut reader = read_tsv(chains_path)?;
    let headers = reader
        .headers()
        .with_context(|| format!("Failed to read TSV header from {}", chains_path.display()))?
        .clone();
    let subject_condition_index = required_index(&headers, "subject:condition", chains_path)?;
    let beta_index = required_index(&headers, "CDR3b", chains_path)?;
    let alpha_index = header_index(&headers, "CDR3a");
    let v_index = header_index(&headers, "TRBV");
    let j_index = header_index(&headers, "TRBJ");
    // Alpha genes are optional; without these columns only the beta genes are used.
    let va_index = header_index(&headers, "TRAV");
    let ja_index = header_index(&headers, "TRAJ");
    let gene = |row: &StringRecord, index: Option<usize>| {
        cell(row, index)
            .filter(|s| !s.is_empty())
            .map(str::to_owned)
    };
    let counts_index = header_index(&headers, "count");
    let mut chains = Vec::new();

    for (row_index, row) in reader.records().enumerate() {
        let row = row.with_context(|| {
            format!(
                "Failed to parse chain row {} in {}",
                row_index + 2,
                chains_path.display()
            )
        })?;
        let subject_condition = cell(&row, Some(subject_condition_index));
        let (subject, condition) = subject_condition
            .and_then(|sc| {
                sc.split_once(':')
                    .map(|(s, c)| (s.trim().to_owned(), c.trim().to_owned()))
            })
            .with_context(|| {
                format!(
                    "Expected 'subject:condition' value on chain row {} in {}",
                    row_index + 2,
                    chains_path.display()
                )
            })?;
        let subject = subject.trim().to_owned();
        let condition = condition.trim().to_owned();
        let beta = match cell(&row, Some(beta_index)) {
            Some(b) => b.to_ascii_uppercase(),
            None => bail!(
                "Empty CDR3β on row {} in {}",
                row_index + 2,
                chains_path.display()
            ),
        };
        if subject.is_empty() || condition.is_empty() {
            bail!(
                "Empty subject or condition in 'subject:condition' on row {} in {}",
                row_index + 2,
                chains_path.display()
            );
        }
        let hla = hla_by_subject
            .as_ref()
            .and_then(|hla_by_subject| hla_by_subject.get(&subject));
        let counts = counts_index
            .map(|index| {
                let value = cell(&row, Some(index));
                match value {
                    Some(value) => {
                        let parsed = value.parse::<f64>().with_context(|| {
                            format!(
                                "Invalid counts value '{value}' on chain row {} in {}",
                                row_index + 2,
                                chains_path.display()
                            )
                        })?;
                        if !parsed.is_finite() || parsed < 0.0 {
                            bail!(
                                "Counts must be a finite non-negative number on chain row {} in {}",
                                row_index + 2,
                                chains_path.display()
                            );
                        }
                        Ok(parsed)
                    }
                    None => Ok(1.0),
                }
            })
            .transpose()?;

        chains.push(ChainRecord {
            subject,
            condition,
            beta,
            alpha: {
                let alpha = cell(&row, alpha_index).map(str::to_ascii_uppercase);
                if let Some(value) = alpha.as_deref()
                    && !value.is_empty()
                    && !crate::tcr_cluster::blosum::is_amino_acid_sequence(value)
                {
                    bail!(
                        "Invalid CDR3a amino-acid sequence on chain row {} in {}",
                        row_index + 2,
                        chains_path.display()
                    );
                }
                alpha
            },
            v_gene: gene(&row, v_index),
            j_gene: gene(&row, j_index),
            va_gene: gene(&row, va_index),
            ja_gene: gene(&row, ja_index),
            counts,
            hla: hla.map_or_else(Default::default, |hla| hla.alleles.clone()),
        });
    }

    if chains.is_empty() {
        bail!("No CDR3 chain records found in {}", chains_path.display());
    }
    Ok(chains)
}

#[cfg(test)]
mod tests {
    use super::read_chain_inputs;
    use pretty_assertions::assert_eq;
    use std::fs;
    use tempfile::tempdir;

    const HLA_HEADER: &str = "PID\tHLA-DPA1(L)\tHLA-DPA1(R)\tHLA-DPB1(L)\tHLA-DPB1(R)\tHLA-DQA1(L)\tHLA-DQA1(R)\tHLA-DQB1(L)\tHLA-DQB1(R)\tHLA-DRB1(L)\tHLA-DRB1(R)\tHLA-DRB3/4/5(L)\tHLA-DRB3/4/5(R)\n";
    const HLA_ROW: &str = "S1\tDPA1*01:03\tDPA1*02:01\tDPB1*04:01\tDPB1*04:01\tDQA1*05:01\tNA\tDQB1*02:01\tDQB1*03:02\tDRB1*15:01\tDRB1*04:03\tDRB3*02:02\tDRB4*01:03\n";

    #[test]
    fn reads_chain_and_hla_tsv_and_joins_by_subject() {
        let temp = tempdir().unwrap();
        let chains = temp.path().join("chains.tsv");
        let hla = temp.path().join("hla.tsv");
        fs::write(
            &chains,
            "subject:condition\tCDR3b\tCDR3a\tTRBV\tTRBJ\tcount\nS1:case\tCASSIRSSYEQYF\tCAVR\tTRBV9\tTRBJ2-1\t4\n",
        )
            .unwrap();
        fs::write(&hla, format!("{HLA_HEADER}{HLA_ROW}")).unwrap();

        let records = read_chain_inputs(&chains, Some(&hla)).unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].subject, "S1");
        assert_eq!(records[0].condition, "case");
        assert_eq!(records[0].beta, "CASSIRSSYEQYF");
        assert_eq!(records[0].hla[0], "DPA1*01:03/DPA1*02:01");
        assert_eq!(records[0].hla[1], "DPB1*04:01");
        assert_eq!(records[0].hla[2], "DQA1*05:01");
        assert_eq!(records[0].hla[3], "DQB1*02:01/DQB1*03:02");
        assert_eq!(records[0].hla[4], "DRB1*15:01/DRB1*04:03");
        assert_eq!(records[0].hla[5], "DRB3*02:02");
        assert_eq!(records[0].hla[6], "DRB4*01:03");
        assert_eq!(records[0].hla[7], "");
        assert_eq!(records[0].counts, Some(4.0));
    }

    #[test]
    fn accepts_hla_rows_with_missing_trailing_columns() {
        let temp = tempdir().unwrap();
        let chains = temp.path().join("chains.tsv");
        let hla = temp.path().join("hla.tsv");
        fs::write(
            &chains,
            "subject:condition\tCDR3b\nS1:case\tCASSIRSSYEQYF\n",
        )
        .unwrap();
        fs::write(
            &hla,
            format!("{HLA_HEADER}S1\tDPA1*01:03\tDPA1*02:01\tDPB1*04:01\n"),
        )
        .unwrap();

        let records = read_chain_inputs(&chains, Some(&hla)).unwrap();

        assert_eq!(records[0].hla[0], "DPA1*01:03/DPA1*02:01");
        assert_eq!(records[0].hla[1], "DPB1*04:01");
        assert_eq!(records[0].hla[7], "");
    }

    #[test]
    fn accepts_chain_subject_missing_from_hla_file() {
        let temp = tempdir().unwrap();
        let chains = temp.path().join("chains.tsv");
        let hla = temp.path().join("hla.tsv");
        fs::write(
            &chains,
            "subject:condition\tCDR3b\nS1:case\tCASSIRSSYEQYF\n",
        )
        .unwrap();
        fs::write(&hla, format!("{HLA_HEADER}{}", HLA_ROW.replace("S1", "S2"))).unwrap();

        let records = read_chain_inputs(&chains, Some(&hla)).unwrap();

        assert_eq!(records.len(), 1);
        assert!(records[0].hla.iter().all(String::is_empty));
    }

    #[test]
    fn reads_chain_inputs_without_hla_file() {
        let temp = tempdir().unwrap();
        let chains = temp.path().join("chains.tsv");
        fs::write(
            &chains,
            "subject:condition\tCDR3b\nS1:case\tCASSIRSSYEQYF\n",
        )
        .unwrap();

        let records = read_chain_inputs(&chains, None).unwrap();

        assert_eq!(records.len(), 1);
        assert!(records[0].hla.iter().all(String::is_empty));
    }

    #[test]
    fn rejects_malformed_subject_condition_value() {
        let temp = tempdir().unwrap();
        let chains = temp.path().join("chains.tsv");
        fs::write(&chains, "subject:condition\tCDR3b\nS1\tCASSIRSSYEQYF\n").unwrap();

        assert!(read_chain_inputs(&chains, None).is_err());
    }

    #[test]
    fn reads_optional_alpha_genes_as_given() {
        let temp = tempdir().unwrap();
        let chains = temp.path().join("chains.tsv");
        fs::write(
            &chains,
            "subject:condition\tCDR3b\tTRBV\tTRBJ\tCDR3a\tTRAV\tTRAJ\n\
             S1:case\tCASSIRSSYEQYF\tTRBV9\tTRBJ2-1\tCAVRDSNYQLIW\tTRAV12-2\tTRAJ33*01\n\
             S2:case\tCASSPRSSYEQYF\tTRBV9\tTRBJ2-1\t\tNA\t\n",
        )
        .unwrap();

        let records = read_chain_inputs(&chains, None).unwrap();

        assert_eq!(records[0].va_gene.as_deref(), Some("TRAV12-2"));
        assert_eq!(records[0].ja_gene.as_deref(), Some("TRAJ33*01"));
        assert_eq!(records[1].va_gene, None);
        assert_eq!(records[1].ja_gene, None);
    }

    #[test]
    fn missing_alpha_gene_columns_default_to_beta_only() {
        let temp = tempdir().unwrap();
        let chains = temp.path().join("chains.tsv");
        fs::write(
            &chains,
            "subject:condition\tCDR3b\tTRBV\tTRBJ\nS1:case\tCASSIRSSYEQYF\tTRBV9\tTRBJ2-1\n",
        )
        .unwrap();

        let records = read_chain_inputs(&chains, None).unwrap();

        assert_eq!(records[0].v_gene.as_deref(), Some("TRBV9"));
        assert!(!records[0].has_alpha_genes());
    }
}
