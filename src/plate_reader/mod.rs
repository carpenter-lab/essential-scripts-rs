use clap::{Error, Subcommand};
use std::path::PathBuf;

#[cfg(feature = "base_cmd")]
mod core;

#[derive(Subcommand)]
#[command(about = "Reformat plate reader data into useful format")]
pub enum Commands {
    ReformatPlateReaderData {
        #[arg(required = true, help = "Input Excel file to process")]
        input_file: PathBuf,

        #[arg(
            required = true,
            help = "Output directory path. Will create one CSV per sheet."
        )]
        output_path: PathBuf,
    },
}
#[cfg(feature = "base_cmd")]
pub fn handle_command(cmd: Commands) -> Result<(), Error> {
    match cmd {
        Commands::ReformatPlateReaderData {
            input_file,
            output_path,
        } => {
            #[cfg(feature = "base_cmd")]
            core::reformat_plate_reader_data(&input_file, &output_path).map_err(|e| {
                Error::raw(
                    clap::error::ErrorKind::Io,
                    format!("Failed to reformat plate reader data: {e}"),
                )
            })?;
        }
    }
    Ok(())
}
