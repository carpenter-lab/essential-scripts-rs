use clap::{Error, Subcommand};
use std::path::PathBuf;

#[cfg(feature = "base_cmd")]
mod core;

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

#[cfg(feature = "base_cmd")]
pub fn handle_command(cmd: Commands) -> Result<(), Error> {
    match cmd {
        Commands::SplitSampleId {
            input_file,
            output_file,
            column_name,
        } => {
            #[cfg(feature = "base_cmd")]
            core::split_sample_id(input_file, output_file, &column_name);
        }
        Commands::SplitCdr3Seq {
            input_file,
            output_file,
            group,
        } => {
            #[cfg(feature = "base_cmd")]
            core::split_cdr3_seq_main(input_file, output_file, group.as_ref());
        }
    }
    Ok(())
}
