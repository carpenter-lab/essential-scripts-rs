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
            let input_file = input_file.ok_or_else(|| {
                Error::raw(
                    clap::error::ErrorKind::MissingRequiredArgument,
                    "Missing <input_file>",
                )
            })?;
            let output_file = output_file.ok_or_else(|| {
                Error::raw(
                    clap::error::ErrorKind::MissingRequiredArgument,
                    "Missing <output_file>",
                )
            })?;

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

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(feature = "tcr")]
    #[test]
    fn handle_command_requires_input_file() {
        let cmd = Commands::ScoreTCRAlignments {
            input_file: None,
            output_file: Some(PathBuf::from("out.csv")),
            replicates: 1,
        };

        let err = handle_command(cmd).expect_err("expected missing input_file error");
        assert_eq!(err.kind(), clap::error::ErrorKind::MissingRequiredArgument);
    }

    #[cfg(feature = "tcr")]
    #[test]
    fn handle_command_requires_output_file() {
        let cmd = Commands::ScoreTCRAlignments {
            input_file: Some(PathBuf::from("in.csv")),
            output_file: None,
            replicates: 1,
        };

        let err = handle_command(cmd).expect_err("expected missing output_file error");
        assert_eq!(err.kind(), clap::error::ErrorKind::MissingRequiredArgument);
    }

    #[cfg(not(feature = "tcr"))]
    #[test]
    fn handle_command_requires_tcr_feature_when_disabled() {
        let cmd = Commands::ScoreTCRAlignments {
            input_file: None,
            output_file: None,
            replicates: 1,
        };

        let err = handle_command(cmd).expect_err("expected feature-gated error");
        assert_eq!(err.kind(), clap::error::ErrorKind::MissingSubcommand);
    }
}
