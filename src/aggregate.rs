use crate::io::WriteToCsvOrStdout;
use clap::Subcommand;
use polars::prelude::*;
use std::path::PathBuf;

#[derive(Subcommand)]
pub enum Commands {
    #[command(about = "Aggregate CellRanger TCR output from multiple samples")]
    AggregateCellRangerTCR {
        #[arg(required = true, num_args = 1.., help = "Input CSV files to process")]
        input_files: Vec<String>,

        #[arg(required = true, help = "Output file path")]
        output_file: PathBuf,

        #[arg(
            short,
            long,
            action = clap::ArgAction::SetTrue,
            default_value_t = false,
            help = "Keep alpha chain in output [default: false]"
        )]
        keep_alpha: bool,
    },
}

fn process_input_lazy(df: LazyFrame, keep_alpha: bool) -> LazyFrame {
    let df2 = df
        .group_by(["sample", "barcode", "chain"])
        .agg([
            col("v_gene").alias("v_gene"),
            col("d_gene").alias("d_gene"),
            col("j_gene").alias("j_gene"),
            col("cdr3").alias("cdr3"),
        ])
        .group_by(["sample", "barcode"])
        .agg([
            col("cdr3")
                .filter(col("chain").eq(lit("TRB")))
                .alias("cdr3_TRB")
                .flatten(),
            col("v_gene")
                .filter(col("chain").eq(lit("TRB")))
                .alias("v_gene_TRB")
                .flatten(),
            col("j_gene")
                .filter(col("chain").eq(lit("TRB")))
                .alias("j_gene_TRB")
                .flatten(),
            col("cdr3")
                .filter(col("chain").eq(lit("TRA")))
                .alias("cdr3_TRA")
                .flatten(),
            col("v_gene")
                .filter(col("chain").eq(lit("TRA")))
                .alias("v_gene_TRA")
                .flatten(),
            col("j_gene")
                .filter(col("chain").eq(lit("TRA")))
                .alias("j_gene_TRA")
                .flatten(),
        ])
        .explode(col("^.*_TRA$").into_selector().unwrap())
        .explode(col("^.*_TRB$").into_selector().unwrap())
        .group_by([col("sample"), col("^.*_TRA$"), col("^.*_TRB$")])
        .agg([len().alias("count")])
        .with_column(col("sample").str().replace_all(lit("-"), lit(":"), false));
    let result = match keep_alpha {
        false => df2.select([
            col("sample"),
            col("cdr3_TRB").alias("CDR3b"),
            col("v_gene_TRB").alias("TRBV"),
            col("j_gene_TRB").alias("TRBJ"),
            col("cdr3_TRA").alias("CDR3a"),
            col("count"),
        ]),
        true => df2.select([
            col("sample"),
            col("cdr3_TRB").alias("CDR3b"),
            col("v_gene_TRB").alias("TRBV"),
            col("j_gene_TRB").alias("TRBJ"),
            col("cdr3_TRA").alias("CDR3a"),
            col("v_gene_TRA").alias("TRAV"),
            col("j_gene_TRA").alias("TRAJ"),
            col("count"),
        ]),
    };
    result
}

pub(crate) fn aggregate_cellranger_tcr_output(
    input_files: &Vec<String>,
    output_file: &PathBuf,
    keep_alpha: bool,
) -> () {
    let concat_args = UnionArgs {
        parallel: true,
        rechunk: true,
        to_supertypes: false,
        diagonal: false,
        from_partitioned_ds: false,
        maintain_order: true,
    };

    println!("Processing {} files", input_files.len());
    let mut dfs = Vec::new();

    for input_file in input_files {
        let df = LazyCsvReader::new(PlPath::new(input_file)).finish();

        match df {
            Ok(df) => dfs.push(df.with_new_streaming(true)),
            Err(error) => {
                eprintln!("Failed to read file: {}", input_file);
                eprintln!("Error: {}", error);
                std::process::exit(1);
            }
        }
    }

    let concatenated = concat(&dfs, concat_args)
        .expect("Failed to concatenate input files")
        .with_new_streaming(true);
    let processed = process_input_lazy(concatenated, keep_alpha)
        .filter(col("CDR3b").is_not_null())
        .filter(col("CDR3b").neq(lit("")));

    processed.write_to_csv_or_stdout(output_file)
}

pub fn handle_command(cmd: Commands) -> () {
    match cmd {
        Commands::AggregateCellRangerTCR {
            input_files,
            output_file,
            keep_alpha,
        } => {
            aggregate_cellranger_tcr_output(&input_files, &output_file, keep_alpha);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use polars_testing::assert_dataframe_equal;

    #[test]
    fn test_process_input_lazy() {
        let df = df!(
            "sample" => [
                "U91846-Inf",
                "U91846-Inf",
                "U91846-Inf",
                "U91846-Inf",
                "U91846-Inf",
                "U90979-Peptide",
                "U90979-Peptide",
                "U90979-Peptide",
                "U90979-Peptide",
                "U90979-Peptide",
                "U90979-Peptide"
            ],
            "barcode" => [
                "ACACTGATCCTTGCCA-1",
                "ACACTGATCCTTGCCA-1",
                "AGATCTGGTGCAGACA-1",
                "AGCATACCACACAGAG-1",
                "AGCATACCACACAGAG-1",
                "AAACCTGTCATCGCTC-1",
                "AAACCTGTCATCGCTC-1",
                "AAGGCAGTCGGTGTCG-1",
                "AAGGCAGTCGGTGTCG-1",
                "AAGGCAGTCGGTGTCG-1",
                "AAGGCAGTCGGTGTCG-1"
            ],
            "chain" => ["TRB", "TRA", "TRB", "TRA", "TRB", "TRB", "TRA", "TRB", "TRB", "TRA", "TRA"],
            "v_gene" => [
                "TRBV7-2",
                "TRAV1-2",
                "TRBV6-5",
                "TRAV10",
                "TRBV15",
                "TRBV27",
                "TRAV12-3",
                "TRBV5-4",
                "TRBV20-1",
                "TRAV41",
                "TRAV8-6"
            ],
            "d_gene" => [None, None, None, None, None, Some("TRBD1"), None, None, None, None, None],
            "j_gene" => [
                "TRBJ1-4",
                "TRAJ47",
                "TRBJ1-3",
                "TRAJ36",
                "TRBJ2-1",
                "TRBJ2-7",
                "TRAJ36",
                "TRBJ2-3",
                "TRBJ1-1",
                "TRAJ47",
                "TRAJ56"
            ],
            "cdr3" => [
                "CASSLGLWEKLFF",
                "CAVRDGNKLVF",
                "CASRRGNTIYF",
                "CVVSSQTGANNLFF",
                "CATSRDRGRNNEQFF",
                "CASSLLGQGWDEQYF",
                "CAMSHQTGANNLFF",
                "CASSLMGGASDTQYF",
                "CSALILNAEAFF",
                "CAVNRNKLVF",
                "CAVSGDTGANSKLTF"
            ]
        )
            .unwrap()
            .lazy();

        let df = process_input_lazy(df, false)
            .collect()
            .unwrap()
            .sort(["sample", "CDR3b", "CDR3a"], Default::default())
            .unwrap();

        let df_expected = df!(
            "sample" => [
                "U90979:Peptide",
                "U90979:Peptide",
                "U91846:Inf",
                "U90979:Peptide",
                "U90979:Peptide",
                "U90979:Peptide",
                "U91846:Inf",
                "U91846:Inf"
            ],
            "CDR3b" => [
                "CASSLMGGASDTQYF",
                "CSALILNAEAFF",
                "CASRRGNTIYF",
                "CASSLMGGASDTQYF",
                "CSALILNAEAFF",
                "CASSLLGQGWDEQYF",
                "CASSLGLWEKLFF",
                "CATSRDRGRNNEQFF"
            ],
            "TRBV" => [
                "TRBV5-4",
                "TRBV20-1",
                "TRBV6-5",
                "TRBV5-4",
                "TRBV20-1",
                "TRBV27",
                "TRBV7-2",
                "TRBV15"
            ],
            "TRBJ" => [
                "TRBJ2-3",
                "TRBJ1-1",
                "TRBJ1-3",
                "TRBJ2-3",
                "TRBJ1-1",
                "TRBJ2-7",
                "TRBJ1-4",
                "TRBJ2-1"
            ],
            "CDR3a" => [
                Some("CAVNRNKLVF"),
                Some("CAVSGDTGANSKLTF"),
                None,
                Some("CAVSGDTGANSKLTF"),
                Some("CAVNRNKLVF"),
                Some("CAMSHQTGANNLFF"),
                Some("CAVRDGNKLVF"),
                Some("CVVSSQTGANNLFF")
            ],
            "count" => [1, 1, 1, 1, 1, 1, 1, 1]
        )
        .unwrap()
        .lazy()
        .with_column(col("count").cast(DataType::UInt32))
        .collect()
        .unwrap()
        .sort(["sample", "CDR3b", "CDR3a"], Default::default())
        .unwrap();

        assert_dataframe_equal!(&df, &df_expected);
    }
}
