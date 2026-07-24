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
