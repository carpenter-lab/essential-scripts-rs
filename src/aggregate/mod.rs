use clap::{Error, Subcommand};
use std::path::PathBuf;

#[cfg(feature = "base_cmd")]
mod core;

#[derive(Subcommand)]
pub enum Commands {
    /// Aggregate CellRanger TCR output from multiple samples
    ///
    /// Parse a set of input files with the Cell Ranger TCR format (`filtered_contig_annotations.csv`) and aggregate them into a single output file.
    /// The output will contain one row per unique combination of sample, barcode, and TCR chain, with the corresponding CDR3 sequences and gene segments.
    /// Optionally, the alpha chain can be retained in the output.
    /// The internal `sample_id` column is used, so there is no need for unique filenames if running on direct outputs from Cell Ranger
    #[command(about, long_about)]
    AggregateCellRangerTCR {
        #[arg(required = true, num_args = 1.., help = "Input CSV files to process")]
        input_files: Vec<PathBuf>,

        #[arg(help = "Output file path")]
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

#[cfg(feature = "base_cmd")]
pub fn handle_command(cmd: Commands) -> Result<(), Error> {
    match cmd {
        Commands::AggregateCellRangerTCR {
            input_files,
            output_file,
            keep_alpha,
        } => {
            #[cfg(feature = "base_cmd")]
            core::aggregate_cellranger_tcr_output(input_files, output_file, keep_alpha);
        }
    }
    Ok(())
}

#[cfg(test)]
#[cfg(feature = "base_cmd")]
mod tests {
    use super::*;
    use crate::test_helpers::{get_timestamp, write_text_file};
    use temp_testdir::TempDir;

    #[test]
    fn test_handle_aggregate_cellranger_tcr_command() {
        let temp_dir = TempDir::default();
        let input_file = temp_dir
            .with_file_name(get_timestamp())
            .with_extension("csv");
        let output_file = temp_dir
            .with_file_name(get_timestamp())
            .with_extension("tsv");

        write_text_file(
            &input_file,
            "sample,barcode,chain,v_gene,d_gene,j_gene,cdr3\n\
             sample1,barcode1,TRB,TRBV1,,TRBJ1,CASSIRSSYEQYF\n",
        );

        let cmd = Commands::AggregateCellRangerTCR {
            input_files: vec![input_file],
            output_file: output_file.clone(),
            keep_alpha: false,
        };

        let result = handle_command(cmd);

        assert!(result.is_ok());
        assert!(output_file.exists());
    }
}
