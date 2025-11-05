use clap::Subcommand;
use std::path::PathBuf;

#[derive(Subcommand)]
#[command(about = "Reformat plate reader data into useful format")]
pub enum Commands {
    ReformatPlateReaderData {
        #[arg(required = true, help = "Input file to process")]
        input_file: PathBuf,

        #[arg(required = true, help = "Output file path")]
        output_file: PathBuf,
    },
}

const LETTERS: [char; 8] = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'];

fn generate_well_ids(rows: &[char], cols: &[u32]) -> Vec<String> {
    rows.iter()
        .flat_map(|r| cols.iter().map(move |c| format!("{}{}", r, c)))
        .collect()
}

pub fn handle_command(cmd: Commands) -> () {
    match cmd {
        Commands::ReformatPlateReaderData {
            input_file: _input_file,
            output_file: _output_file,
        } => {
            todo!()
        }
    }
}
