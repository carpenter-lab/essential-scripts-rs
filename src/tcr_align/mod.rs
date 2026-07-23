#[cfg(feature = "tcr")]
mod align;
#[cfg(feature = "tcr")]
mod core;
#[cfg(feature = "tcr")]
mod dataframe;

use clap::{Error, Subcommand};
use std::path::PathBuf;

#[cfg(not(feature = "tcr"))]
const CMD_ABOUT: &str = "Score a GLIPH2 output file with a TCR alignment pipeline. Requires re-installing with --features tcr";
#[cfg(feature = "tcr")]
const CMD_ABOUT: &str = "Score a GLIPH2 output file with a TCR alignment pipeline";

#[derive(Subcommand)]
pub enum Commands {
    #[command(
        about = CMD_ABOUT,
    )]
    ScoreTCRAlignments {
        #[arg(help = "Input CSV file to process")]
        input_file: Option<PathBuf>,

        #[arg(help = "Output file path")]
        output_file: Option<PathBuf>,

        #[arg(short, long, default_value_t = 1000)]
        replicates: usize,
    },
}

pub fn handle_command(cmd: Commands) -> Result<(), Error> {
    #[cfg(feature = "tcr")]
    match cmd {
        Commands::ScoreTCRAlignments {
            input_file,
            output_file,
            replicates,
        } => {
            core::tcr_score(input_file, output_file, replicates);
            Ok(())
        }
    }
    #[cfg(not(feature = "tcr"))]
    match cmd {
        Commands::ScoreTCRAlignments { .. } => Err(Error::raw(
            clap::error::ErrorKind::MissingSubcommand,
            "This command requires the `tcr` feature. Rebuild with `--features tcr`",
        )),
    }
}
