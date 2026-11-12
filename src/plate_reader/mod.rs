use anyhow::{Context, Result};
use clap::Subcommand;
use std::path::PathBuf;

mod core;
mod dataframe;

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
pub fn handle_command(cmd: Commands) -> Result<()> {
    match cmd {
        Commands::ReformatPlateReaderData {
            input_file,
            output_path,
        } => {
            core::reformat_plate_reader_data(&input_file, &output_path)
                .context("Failed to reformat plate reader data")?;
        }
    }
    Ok(())
}
