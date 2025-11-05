use crate::io::WriteToCsvOrStdout;
use clap::Subcommand;
use polars::prelude::pivot::pivot_stable;
use polars::prelude::*;
use std::error::Error;
use std::path::PathBuf;

#[derive(Subcommand)]
pub enum Commands {
    #[command(about = "Aggregate CellRanger TCR output from multiple samples")]
    AggregateCellRangerTCR {
        #[arg(required = true, num_args = 1.., help = "Input CSV files to process")]
        input_files: Vec<PathBuf>,

        #[arg(required = true, help = "Output file path")]
        output_file: PathBuf,
    },
}

fn process_input_lazy(df: LazyFrame) -> Result<LazyFrame, PolarsError> {
    // Collect and pivot
    let df2 = df.group_by(["sample", "barcode", "chain"]).agg([
        col("v_gene").alias("v_gene"),
        col("d_gene").alias("d_gene"),
        col("j_gene").alias("j_gene"),
        col("cdr3").alias("cdr3"),
    ]);

    let pivoted = pivot_stable(
        &df2.collect()?,
        ["chain"],
        Some(["sample", "barcode"]),
        Some(["v_gene", "d_gene", "j_gene", "cdr3"]),
        false,
        None,
        None,
    )?;

    let result = pivoted
        .lazy()
        .explode(Selector::Matches(PlSmallStr::from_str("_TRB$")))
        .explode(Selector::Matches(PlSmallStr::from_str("_TRA$")))
        .group_by([
            col("sample"),
            col("barcode"),
            col("cdr3_TRB"),
            col("v_gene_TRB"),
            col("j_gene_TRB"),
            col("cdr3_TRA"),
            col("v_gene_TRA"),
            col("j_gene_TRA"),
        ])
        .agg([len().alias("count")])
        .with_column(col("sample").str().replace_all(lit("-"), lit(":"), false))
        .select([
            col("sample"),
            col("cdr3_TRB").alias("CDR3b"),
            col("v_gene_TRB").alias("TRBV"),
            col("j_gene_TRB").alias("TRBJ"),
            col("cdr3_TRA").alias("CDR3a"),
            col("count"),
        ]);

    Ok(result)
}

pub(crate) fn aggregate_cellranger_tcr_output(
    input_files: &Vec<PathBuf>,
    output_file: &PathBuf,
) -> Result<(), Box<dyn Error>> {
    println!("Processing {:?} files", input_files);
    let mut dfs = Vec::new();

    for input_file in input_files {
        let df = LazyCsvReader::new(PlPath::from_str(input_file.to_str().expect(""))).finish()?;
        dfs.push(df);
    }

    let concatenated = concat(
        &dfs,
        UnionArgs {
            parallel: true,
            rechunk: false,
            to_supertypes: false,
            diagonal: false,
            from_partitioned_ds: false,
            maintain_order: false,
        },
    )?;
    let processed = process_input_lazy(concatenated)?;

    processed.write_to_csv_or_stdout(output_file)
}

pub fn handle_command(cmd: Commands) -> () {
    match cmd {
        Commands::AggregateCellRangerTCR {
            input_files,
            output_file,
        } => {
            aggregate_cellranger_tcr_output(&input_files, &output_file).expect("");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let paths: Vec<PathBuf> = Vec::from([
            "filtered_contig_annotations.csv",
            "filtered_contig_annotations1.csv",
        ])
        .into_iter()
        .map(PathBuf::from)
        .collect();
        let result = aggregate_cellranger_tcr_output(&paths, &PathBuf::from("test_out.csv"));
        assert!(result.is_ok());
        let result_data = LazyCsvReader::new(PlPath::from_str("test_out.csv"))
            .finish()
            .unwrap()
            .sort(["sample", "CDR3b"], Default::default());
        let expected_data = LazyCsvReader::new(PlPath::from_str("test_out_expected.csv"))
            .finish()
            .unwrap()
            .sort(["sample", "CDR3b"], Default::default());
        assert_eq!(
            result_data.collect().unwrap(),
            expected_data.collect().unwrap()
        );
    }
}
