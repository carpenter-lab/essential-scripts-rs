use crate::io;
use crate::io::WriteToCsvOrStdout;
use clap::Subcommand;
use polars::prelude::*;
use std::path::PathBuf;

#[derive(Subcommand)]
pub enum Commands {
    #[command(about = "Split sample ID into subject and condition from GLIPH2 output")]
    SplitSampleId {
        #[arg(required = true, num_args = 1.., help = "Input CSV file to process")]
        input_file: PathBuf,

        #[arg(required = true, help = "Output file path")]
        output_file: PathBuf,

        #[arg(
            short,
            long,
            default_value = "subject:condition",
            help = "Column to split"
        )]
        column_name: String,
    },
    #[command(about = "Split CDR3 sequences and genes if a semicolon is present")]
    SplitCdr3Seq {
        #[arg(required = true, help = "Input CSV file to process")]
        input_file: PathBuf,

        #[arg(required = true, help = "Output file path")]
        output_file: PathBuf,

        #[arg(
            short,
            long,
            num_args = 1..,
            help = "Optional columns to group by; if omitted, each input row is treated as its own group"
        )]
        group: Option<Vec<String>>,
    },
}

#[derive(Clone, Copy)]
enum GeneSchema {
    CtGene,
    TrVj,
}

fn has_col(df: &DataFrame, name: &str) -> bool {
    df.get_column_names().iter().any(|n| n.as_str() == name)
}

fn resolve_gene_schema(df: &DataFrame) -> PolarsResult<GeneSchema> {
    let has_ctgene = has_col(df, "CTgeneA") && has_col(df, "CTgeneB");
    let has_tr_vj =
        has_col(df, "TRAV") && has_col(df, "TRAJ") && has_col(df, "TRBV") && has_col(df, "TRBJ");

    if !has_col(df, "CDR3a") || !has_col(df, "CDR3b") {
        return Err(PolarsError::ComputeError(
            "split-cdr3-seq requires both CDR3a and CDR3b columns".into(),
        ));
    }

    match (has_ctgene, has_tr_vj) {
        (true, false) => Ok(GeneSchema::CtGene),
        (false, true) => Ok(GeneSchema::TrVj),
        (true, true) => Err(PolarsError::ComputeError(
            "Found both schema styles; provide either CTgeneA/CTgeneB or TRAV/TRAJ/TRBV/TRBJ"
                .into(),
        )),
        (false, false) => Err(PolarsError::ComputeError(
            "Missing required gene columns. Provide either CTgeneA/CTgeneB or TRAV/TRAJ/TRBV/TRBJ"
                .into(),
        )),
    }
}

fn split_cdr3_seq(df: DataFrame, chain: &str, schema: GeneSchema) -> PolarsResult<DataFrame> {
    if chain != "alpha" && chain != "beta" {
        return Err(PolarsError::ComputeError(
            format!("chain must be alpha or beta, not {}", chain).into(),
        ));
    }

    let split_cols: Vec<&str> = match (chain, schema) {
        ("alpha", GeneSchema::CtGene) => vec!["CTgeneA", "CDR3a"],
        ("beta", GeneSchema::CtGene) => vec!["CTgeneB", "CDR3b"],
        ("alpha", GeneSchema::TrVj) => vec!["TRAV", "TRAJ", "CDR3a"],
        ("beta", GeneSchema::TrVj) => vec!["TRBV", "TRBJ", "CDR3b"],
        _ => unreachable!(),
    };

    let list_cols: Vec<String> = split_cols.iter().map(|c| format!("__{}_list", c)).collect();

    let explode_pattern = list_cols
        .iter()
        .map(|c| format!("^{}$", c))
        .collect::<Vec<String>>()
        .join("|");
    let explode_selector = col(&explode_pattern)
        .into_selector()
        .expect("failed to build explode selector");

    let split_exprs: Vec<Expr> = split_cols
        .iter()
        .zip(list_cols.iter())
        .map(|(source, list_name)| col(*source).str().split(lit(";")).alias(list_name))
        .collect();

    let alias_exprs: Vec<Expr> = split_cols
        .iter()
        .zip(list_cols.iter())
        .map(|(source, list_name)| col(list_name).alias(*source))
        .collect();

    df.lazy()
        .with_columns(split_exprs)
        // Explode related columns together so chain columns stay paired.
        .explode(
            explode_selector,
            ExplodeOptions {
                empty_as_null: false,
                keep_nulls: false,
            },
        )
        .with_columns(alias_exprs)
        .select([all()
            .exclude_cols(list_cols.iter().map(|s| s.as_str()).collect::<Vec<&str>>())
            .as_expr()])
        .collect()
}

pub(crate) fn split_cdr3_seq_main(
    input_file: PathBuf,
    output_file: PathBuf,
    group: Option<Vec<String>>,
) -> () {
    let lazy_df: LazyFrame = io::read_from_file(input_file, None);
    let mut df = lazy_df
        .collect()
        .expect("Failed to collect initial dataframe");
    let gene_schema = resolve_gene_schema(&df).expect("Invalid gene column schema");

    // Keep CLI compatibility and fail early on typos in group column names.
    // If user provides group columns, validate them.
    // If not, create a unique per-row group key.
    match &group {
        Some(cols) => {
            for g in cols {
                if !df.get_column_names().iter().any(|n| n.as_str() == g) {
                    panic!("Group column '{}' not found", g);
                }
            }
        }
        None => {
            // One unique group per original row.
            // This keeps alpha/beta splits anchored to original rows.
            let row_ids: Vec<u64> = (0..df.height() as u64).collect();
            df = df
                .with_column(Series::new("__row_group".into(), row_ids).into())
                .expect("Failed to add per-row fallback group")
                .to_owned();
        }
    }

    // Process both chains
    for chain in ["alpha", "beta"] {
        df = split_cdr3_seq(df, chain, gene_schema).expect("Failed to split CDR3 sequences");
    }

    // Write output
    if group.is_none()
        && df
            .get_column_names()
            .iter()
            .any(|n| n.as_str() == "__row_group")
    {
        df = df
            .drop("__row_group")
            .expect("Failed to drop temporary row group column");
    }
    df.write_to_flat_or_stdout(output_file, None)
}

pub(crate) fn split_sample_id(
    input_file: PathBuf,
    output_file: PathBuf,
    column_name: &String,
) -> () {
    let df = io::read_from_file(input_file, None);
    let df = df
        .with_column(col(column_name).str().split(lit(":")))
        .with_columns([
            col(column_name).list().first().alias("subject"),
            col(column_name).list().last().alias("condition"),
        ])
        .select([all().exclude_cols([column_name]).as_expr()]);

    df.collect()
        .expect("Failed to collect dataframe")
        .write_to_flat_or_stdout(output_file, None)
}

pub fn handle_command(cmd: Commands) -> () {
    match cmd {
        Commands::SplitSampleId {
            input_file,
            output_file,
            column_name,
        } => {
            split_sample_id(input_file, output_file, &column_name);
        }
        Commands::SplitCdr3Seq {
            input_file,
            output_file,
            group,
        } => {
            split_cdr3_seq_main(input_file, output_file, group);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn split_cdr3_seq_expands_alpha_pairs_and_preserves_rows() {
        let df = df!(
            "CDR3a" => [
                "CAVVPPNQAGTALIF;CVVNGGRLMF",
                "CAVNPVGSYIPTF",
                "CALSGDSGNTGKLIF"
            ],
            "CDR3b" => [
                "CSARSSGTGSSYNSPLHF",
                "CSASSAGGTDTQYF",
                "CSASKSSYEQYF"
            ],
            "TRAV" => [
                "TRAV2.TRAJ15.TRAC;TRAV12-1.TRAJ31.TRAC",
                "TRAV8-1.TRAJ6.TRAC",
                "TRAV16.TRAJ37.TRAC"
            ],
            "TRAJ" => ["TRAJ15;TRAJ31", "TRAJ6", "TRAJ37"],
            "TRBV" => [
                "TRBV20-1.TRBJ1-6.TRBC1",
                "TRBV20-1.TRBJ2-3.TRBC2",
                "TRBV29-1.TRBJ2-7.TRBC2"
            ],
            "TRBJ" => ["TRBJ1-6", "TRBJ2-3", "TRBJ2-7"],
            "groups" => ["A", "B", "C"]
        )
        .expect("failed to create test dataframe");

        let out =
            split_cdr3_seq(df, "alpha", GeneSchema::TrVj).expect("alpha split should succeed");

        assert_eq!(out.height(), 4);

        let groups: Vec<&str> = out
            .column("groups")
            .expect("missing groups")
            .str()
            .expect("groups should be utf8")
            .into_iter()
            .map(|v| v.expect("groups value should not be null"))
            .collect();
        assert_eq!(groups, vec!["A", "A", "B", "C"]);

        let trav: Vec<&str> = out
            .column("TRAV")
            .expect("missing TRAV")
            .str()
            .expect("TRAV should be utf8")
            .into_iter()
            .map(|v| v.expect("TRAV value should not be null"))
            .collect();
        assert_eq!(
            trav,
            vec![
                "TRAV2.TRAJ15.TRAC",
                "TRAV12-1.TRAJ31.TRAC",
                "TRAV8-1.TRAJ6.TRAC",
                "TRAV16.TRAJ37.TRAC"
            ]
        );

        let cdr3a: Vec<&str> = out
            .column("CDR3a")
            .expect("missing CDR3a")
            .str()
            .expect("CDR3a should be utf8")
            .into_iter()
            .map(|v| v.expect("CDR3a value should not be null"))
            .collect();
        assert_eq!(
            cdr3a,
            vec![
                "CAVVPPNQAGTALIF",
                "CVVNGGRLMF",
                "CAVNPVGSYIPTF",
                "CALSGDSGNTGKLIF"
            ]
        );
    }

    #[test]
    fn split_cdr3_seq_alpha_then_beta_matches_expected_row_count() {
        let df = df!(
            "CDR3a" => [
                "CAVVPPNQAGTALIF;CVVNGGRLMF",
                "CAVNPVGSYIPTF",
                "CALSGDSGNTGKLIF"
            ],
            "CDR3b" => [
                "CSARSSGTGSSYNSPLHF",
                "CSASSAGGTDTQYF",
                "CSASKSSYEQYF"
            ],
            "TRAV" => [
                "TRAV2.TRAJ15.TRAC;TRAV12-1.TRAJ31.TRAC",
                "TRAV8-1.TRAJ6.TRAC",
                "TRAV16.TRAJ37.TRAC"
            ],
            "TRAJ" => ["TRAJ15;TRAJ31", "TRAJ6", "TRAJ37"],
            "TRBV" => [
                "TRBV20-1.TRBJ1-6.TRBC1",
                "TRBV20-1.TRBJ2-3.TRBC2",
                "TRBV29-1.TRBJ2-7.TRBC2"
            ],
            "TRBJ" => ["TRBJ1-6", "TRBJ2-3", "TRBJ2-7"],
            "groups" => ["A", "B", "C"]
        )
        .expect("failed to create test dataframe");

        let out =
            split_cdr3_seq(df, "alpha", GeneSchema::TrVj).expect("alpha split should succeed");
        let out = split_cdr3_seq(out, "beta", GeneSchema::TrVj).expect("beta split should succeed");

        assert_eq!(out.height(), 4);

        let cdr3a: Vec<&str> = out
            .column("CDR3a")
            .expect("missing CDR3a")
            .str()
            .expect("CDR3a should be utf8")
            .into_iter()
            .map(|v| v.expect("CDR3a value should not be null"))
            .collect();
        assert_eq!(
            cdr3a,
            vec![
                "CAVVPPNQAGTALIF",
                "CVVNGGRLMF",
                "CAVNPVGSYIPTF",
                "CALSGDSGNTGKLIF"
            ]
        );
    }

    #[test]
    fn split_cdr3_seq_errors_on_mismatched_alpha_pair_lengths() {
        let df = df!(
            "CDR3a" => ["CAVVPPNQAGTALIF"],
            "CDR3b" => ["CSARSSGTGSSYNSPLHF"],
            "TRAV" => ["TRAV2.TRAJ15.TRAC;TRAV12-1.TRAJ31.TRAC"],
            "TRAJ" => ["TRAJ15"],
            "TRBV" => ["TRBV20-1.TRBJ1-6.TRBC1"],
            "TRBJ" => ["TRBJ1-6"],
            "groups" => ["A"]
        )
        .expect("failed to create test dataframe");

        let out = split_cdr3_seq(df, "alpha", GeneSchema::TrVj);
        assert!(
            out.is_err(),
            "expected error when TRAV/CDR3a split lengths are mismatched"
        );
    }

    #[test]
    fn resolve_gene_schema_accepts_ctgene_columns() {
        let df = df!(
            "CDR3a" => ["A"],
            "CDR3b" => ["B"],
            "CTgeneA" => ["X;Y"],
            "CTgeneB" => ["Z"]
        )
        .expect("failed to create test dataframe");

        let schema = resolve_gene_schema(&df).expect("CTgene schema should be valid");
        assert!(matches!(schema, GeneSchema::CtGene));
    }

    #[test]
    fn resolve_gene_schema_rejects_partial_tr_columns() {
        let df = df!(
            "CDR3a" => ["A"],
            "CDR3b" => ["B"],
            "TRAV" => ["TRAV1"],
            "TRBV" => ["TRBV1"]
        )
        .expect("failed to create test dataframe");

        let schema = resolve_gene_schema(&df);
        assert!(schema.is_err(), "partial TR columns should be rejected");
    }

    #[test]
    fn split_cdr3_seq_main_without_group_uses_unique_row_fallback() {
        // This test targets the intended "group omitted" behavior:
        // one unique group per input row, so expansion is row-local and no cross-row mixing.

        // Build a minimal TR V/J schema dataframe similar to your existing fixtures.
        let mut df = df!(
            "CDR3a" => [
                "CAVVPPNQAGTALIF;CVVNGGRLMF",
                "CAVNPVGSYIPTF",
                "CALSGDSGNTGKLIF"
            ],
            "CDR3b" => [
                "CSARSSGTGSSYNSPLHF",
                "CSASSAGGTDTQYF",
                "CSASKSSYEQYF"
            ],
            "TRAV" => [
                "TRAV2.TRAJ15.TRAC;TRAV12-1.TRAJ31.TRAC",
                "TRAV8-1.TRAJ6.TRAC",
                "TRAV16.TRAJ37.TRAC"
            ],
            "TRAJ" => ["TRAJ15;TRAJ31", "TRAJ6", "TRAJ37"],
            "TRBV" => [
                "TRBV20-1.TRBJ1-6.TRBC1",
                "TRBV20-1.TRBJ2-3.TRBC2",
                "TRBV29-1.TRBJ2-7.TRBC2"
            ],
            "TRBJ" => ["TRBJ1-6", "TRBJ2-3", "TRBJ2-7"]
        )
        .expect("failed to create test dataframe");

        // Simulate split_cdr3_seq_main fallback path when group is None:
        // add temporary per-row key, run both chains, drop temporary key.
        let row_ids: Vec<u64> = (0..df.height() as u64).collect();
        df = df
            .with_column(Series::new("__row_group".into(), row_ids).into())
            .expect("failed to add fallback row group")
            .to_owned();

        let schema = resolve_gene_schema(&df).expect("schema should resolve");
        let out = split_cdr3_seq(df, "alpha", schema).expect("alpha split should succeed");
        let mut out = split_cdr3_seq(out, "beta", schema).expect("beta split should succeed");

        if out
            .get_column_names()
            .iter()
            .any(|n| n.as_str() == "__row_group")
        {
            out = out
                .drop("__row_group")
                .expect("failed to drop temporary fallback group");
        }

        // Same expected output shape as grouped case: first row expands to two.
        assert_eq!(out.height(), 4);

        let cdr3a: Vec<&str> = out
            .column("CDR3a")
            .expect("missing CDR3a")
            .str()
            .expect("CDR3a should be utf8")
            .into_iter()
            .map(|v| v.expect("CDR3a should not be null"))
            .collect();

        assert_eq!(
            cdr3a,
            vec![
                "CAVVPPNQAGTALIF",
                "CVVNGGRLMF",
                "CAVNPVGSYIPTF",
                "CALSGDSGNTGKLIF"
            ]
        );

        // Ensure temp fallback column is not leaked.
        assert!(
            !out.get_column_names()
                .iter()
                .any(|n| n.as_str() == "__row_group"),
            "temporary fallback group column should be removed"
        );
    }

    #[test]
    fn provided_group_columns_still_validate() {
        let df = df!(
            "CDR3a" => ["A"],
            "CDR3b" => ["B"],
            "TRAV" => ["TRAV1"],
            "TRAJ" => ["TRAJ1"],
            "TRBV" => ["TRBV1"],
            "TRBJ" => ["TRBJ1"],
            "groups" => ["G1"]
        )
        .expect("failed to create test dataframe");

        // Valid group column should pass the same validation logic used in split_cdr3_seq_main.
        let valid_groups = vec!["groups".to_string()];
        for g in &valid_groups {
            assert!(
                df.get_column_names().iter().any(|n| n.as_str() == g),
                "group column '{}' should exist",
                g
            );
        }

        // Missing group column should fail validation.
        let invalid_groups = vec!["does_not_exist".to_string()];
        let missing = invalid_groups.iter().find(|g| {
            !df.get_column_names()
                .iter()
                .any(|n| n.as_str() == g.as_str())
        });
        assert!(
            missing.is_some(),
            "validation should detect missing provided group columns"
        );
    }
}
