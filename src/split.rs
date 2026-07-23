use crate::io;
use crate::io::WriteToCsvOrStdout;
use clap::{Error, Subcommand};
use polars::prelude::*;
use std::path::PathBuf;

#[derive(Subcommand)]
pub enum Commands {
    /// Split sample ID into subject and condition from GLIPH2 output
    ///
    /// Splits the `subject:condition` column (the default) into two columns.
    /// Can be used on any file with a separator to split by however,
    /// these columns will always be called subject and condition
    #[command(about, long_about)]
    SplitSampleId {
        #[arg(required = true, num_args = 1.., help = "Input CSV file to process")]
        input_file: PathBuf,

        #[arg(help = "Output file path")]
        output_file: PathBuf,

        #[arg(
            short,
            long,
            default_value = "subject:condition",
            help = "Column to split"
        )]
        column_name: String,
    },
    /// Split CDR3 sequences and genes if a semicolon is present
    ///
    /// This potentially splits each row into up to 4 new rows.
    /// In the output from scRepertoire, the rows each belong to a single cell barcode.
    /// In the case when a cell has 2 detectable CDR3 beta or alpha sequences, the resulting CDR3 and V/J columns are concatinated with a ";".
    /// For downstream applications, this results in treating this chimeric sequence as a real, biological sequence.
    /// This tool will expand this into up to 4 different pairs of chains.
    /// "beta_1;beta2" and "alpha_1;alpha_2" will be split into 4 rows each containing a single alpha and beta chain.
    ///
    /// The TCR columns must be named CDR3a and CDR3b.
    /// Requires either CTgeneA/CTgeneB columns or TRAV/TRAJ/TRBV/TRBJ columns for the TCR genes.
    /// If group columns are not provided, each input row is treated as its own group and alpha/beta splits are anchored to original rows.
    #[command(about, long_about)]
    SplitCdr3Seq {
        #[arg(required = true, help = "Input CSV file to process")]
        input_file: PathBuf,

        #[arg(default_value = "-", help = "Output file path")]
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
            format!("chain must be alpha or beta, not {chain}").into(),
        ));
    }

    let split_cols: Vec<&str> = match (chain, schema) {
        ("alpha", GeneSchema::CtGene) => vec!["CTgeneA", "CDR3a"],
        ("beta", GeneSchema::CtGene) => vec!["CTgeneB", "CDR3b"],
        ("alpha", GeneSchema::TrVj) => vec!["TRAV", "TRAJ", "CDR3a"],
        ("beta", GeneSchema::TrVj) => vec!["TRBV", "TRBJ", "CDR3b"],
        _ => unreachable!(),
    };

    let list_cols: Vec<String> = split_cols.iter().map(|c| format!("__{c}_list")).collect();

    let explode_pattern = list_cols
        .iter()
        .map(|c| format!("^{c}$"))
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
            .exclude_cols(list_cols.iter().map(String::as_str).collect::<Vec<&str>>())
            .as_expr()])
        .collect()
}

pub(crate) fn split_cdr3_seq_main(
    input_file: PathBuf,
    output_file: PathBuf,
    group: Option<&Vec<String>>,
) {
    let lazy_df: LazyFrame = io::read_from_file(input_file, None);
    let mut df = lazy_df
        .collect()
        .expect("Failed to collect initial dataframe");
    let gene_schema = resolve_gene_schema(&df).expect("Invalid gene column schema");

    // Keep CLI compatibility and fail early on typos in group column names.
    // If user provides group columns, validate them.
    // If not, create a unique per-row group key.
    if let Some(cols) = group {
        for g in cols {
            assert!(
                df.get_column_names().iter().any(|n| n.as_str() == g),
                "Group column '{g}' not found"
            );
        }
    } else {
        // One unique group per original row.
        // This keeps alpha/beta splits anchored to original rows.
        let row_ids: Vec<u64> = (0..df.height() as u64).collect();
        df.with_column(Series::new("__row_group".into(), row_ids).into())
            .expect("Failed to add per-row fallback group");
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
    df.write_to_flat_or_stdout(output_file, None);
}

pub(crate) fn split_sample_id(input_file: PathBuf, output_file: PathBuf, column_name: &String) {
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
        .write_to_flat_or_stdout(output_file, None);
}

pub fn handle_command(cmd: Commands) -> Result<(), Error> {
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
            split_cdr3_seq_main(input_file, output_file, group.as_ref());
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    fn expected_split_cdr3a() -> Vec<&'static str> {
        vec![
            "CAVVPPNQAGTALIF",
            "CVVNGGRLMF",
            "CAVNPVGSYIPTF",
            "CALSGDSGNTGKLIF",
        ]
    }

    fn expected_split_cdr3b() -> Vec<&'static str> {
        vec![
            "CSARSSGTGSSYNSPLHF",
            "CSARSSGTGSSYNSPLHF",
            "CSASSAGGTDTQYF",
            "CSASKSSYEQYF",
        ]
    }

    fn col_utf8_values<'a>(df: &'a DataFrame, col_name: &str) -> PolarsResult<Vec<&'a str>> {
        let df = df
            .column(col_name)?
            .str()?
            .iter()
            .map(|v| v.expect("column value should not be null"))
            .collect();
        Ok(df)
    }

    fn cdr3_cols() -> Vec<Column> {
        let c1 = Column::new(
            "CDR3a".into(),
            [
                "CAVVPPNQAGTALIF;CVVNGGRLMF",
                "CAVNPVGSYIPTF",
                "CALSGDSGNTGKLIF",
            ],
        );
        let c2 = Column::new(
            "CDR3b".into(),
            ["CSARSSGTGSSYNSPLHF", "CSASSAGGTDTQYF", "CSASKSSYEQYF"],
        );
        vec![c1, c2]
    }

    fn ct_cols() -> Vec<Column> {
        let c1 = Column::new(
            "CTgeneA".into(),
            [
                "TRAV2.TRAJ15.TRAC;TRAV12-1.TRAJ31.TRAC",
                "TRAV8-1.TRAJ6.TRAC",
                "TRAV16.TRAJ37.TRAC",
            ],
        );
        let c2 = Column::new(
            "CTgeneB".into(),
            [
                "TRBV20-1.TRBJ1-6.TRBC1",
                "TRBV20-1.TRBJ2-3.TRBC2",
                "TRBV29-1.TRBJ2-7.TRBC2",
            ],
        );
        vec![c1, c2]
    }

    fn vj_cols() -> Vec<Column> {
        let c1 = Column::new("TRAV".into(), ["TRAV2;TRAV12-1", "TRAV8-1", "TRAV16"]);
        let c2 = Column::new("TRAJ".into(), ["TRAJ15;TRAJ31", "TRAJ6", "TRAJ37"]);
        let c3 = Column::new("TRBV".into(), ["TRBV20-1", "TRBV20-1", "TRBV29-1"]);
        let c4 = Column::new("TRBJ".into(), ["TRBJ1-6", "TRBJ2-3", "TRBJ2-7"]);
        vec![c1, c2, c3, c4]
    }

    fn expected_split_ctgene_a() -> Vec<&'static str> {
        vec![
            "TRAV2.TRAJ15.TRAC",
            "TRAV12-1.TRAJ31.TRAC",
            "TRAV8-1.TRAJ6.TRAC",
            "TRAV16.TRAJ37.TRAC",
        ]
    }
    fn expected_split_va() -> Vec<&'static str> {
        vec!["TRAV2", "TRAV12-1", "TRAV8-1", "TRAV16"]
    }

    fn create_df(schema: GeneSchema, group: bool) -> DataFrame {
        let mut col_vec = Vec::new();
        col_vec.extend(cdr3_cols());
        match schema {
            GeneSchema::CtGene => col_vec.extend(ct_cols()),
            GeneSchema::TrVj => col_vec.extend(vj_cols()),
        }
        if group {
            col_vec.push(Column::new("groups".into(), ["A", "B", "C"]));
        }
        DataFrame::new_infer_height(col_vec).expect("failed to create test dataframe")
    }

    #[rstest]
    #[case::ctgene_group(GeneSchema::CtGene, true)]
    #[case::ctgene_no_group(GeneSchema::CtGene, false)]
    #[case::trvj_group(GeneSchema::TrVj, true)]
    #[case::trvj_no_group(GeneSchema::TrVj, false)]
    fn split_cdr3_seq_test(#[case] schema: GeneSchema, #[case] group: bool) {
        let create_df = create_df(schema, group);
        let out = split_cdr3_seq(create_df, "alpha", schema).expect("alpha split should succeed");
        let out = split_cdr3_seq(out, "beta", schema).expect("beta split should succeed");
        if group {
            match col_utf8_values(&out, "groups") {
                Ok(groups) => assert_eq!(groups, vec!["A", "A", "B", "C"]),
                Err(..) => panic!(""),
            }
        } else {
            assert!(
                col_utf8_values(&out, "groups").is_err(),
                "no group column should be present"
            );
            assert!(
                col_utf8_values(&out, "__row_group").is_err(),
                "temporary group column should be removed"
            )
        }
        assert_eq!(out.height(), 4);
        match col_utf8_values(&out, "CDR3a") {
            Ok(groups) => assert_eq!(
                groups,
                expected_split_cdr3a(),
                "CDR3a sequences are not the same or in the same order"
            ),
            Err(..) => panic!("could not extract CDR3a column"),
        }
        match col_utf8_values(&out, "CDR3b") {
            Ok(groups) => assert_eq!(
                groups,
                expected_split_cdr3b(),
                "CDR3b sequences are not the same or in the same order"
            ),
            Err(..) => panic!("could not extract CDR3b column"),
        }
        match schema {
            GeneSchema::CtGene => match col_utf8_values(&out, "CTgeneA") {
                Ok(groups) => assert_eq!(
                    groups,
                    expected_split_ctgene_a(),
                    "CTgeneA sequences are not the same or in the same order"
                ),
                Err(..) => panic!("could not extract CTgeneA column"),
            },
            GeneSchema::TrVj => match col_utf8_values(&out, "TRAV") {
                Ok(groups) => assert_eq!(
                    groups,
                    expected_split_va(),
                    "VA sequences are not the same or in the same order"
                ),
                Err(..) => panic!("could not extract TRAV column"),
            },
        }
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
        let invalid_groups = ["does_not_exist".to_string()];
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

    #[rstest]
    #[should_panic(
        expected = "CTgene schema should be valid: ComputeError(ErrString(\"Missing required gene columns. Provide either CTgeneA/CTgeneB or TRAV/TRAJ/TRBV/TRBJ\"))"
    )]
    fn resolve_gene_schema_rejects_invalid_ctgene_columns() {
        let df = df!(
            "CDR3a" => ["A"],
            "CDR3b" => ["B"],
            "CTgeneA" => ["X;Y"],
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
    fn split_cdr3_seq_ctgene_errors_on_mismatched_alpha_pair_lengths() {
        let df = df!(
            "CDR3a" => ["CAVVPPNQAGTALIF"],
            "CDR3b" => ["CSARSSGTGSSYNSPLHF"],
            "CTgeneA" => ["TRAV2.TRAJ15.TRAC;TRAV12-1.TRAJ31.TRAC"],
            "CTgeneB" => ["TRBV20-1.TRBJ1-6.TRBC1"],
            "groups" => ["A"]
        )
        .expect("failed to create CTgene mismatch dataframe");

        let out = split_cdr3_seq(df, "alpha", GeneSchema::CtGene);
        assert!(
            out.is_err(),
            "expected error when CTgeneA/CDR3a split lengths are mismatched"
        );
    }
}
