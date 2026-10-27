#[cfg(feature = "base_cmd")]
pub mod helper;
use crate::cores::{thread_default, thread_parser};
use clap::Subcommand;
use std::path::PathBuf;

#[derive(Subcommand)]
pub enum Commands {
    /// A reimplementation of TCRMatch
    TcrMatch {
        #[arg(help = "Input file to process")]
        input_file: PathBuf,

        #[arg(short = 'd', long, default_value_t = helper::Database::Iedb)]
        database: helper::Database,

        #[arg(
            long,
            short = 't',
            value_parser = thread_parser(),
            default_value_t = thread_default(),
            help = "Number of cores to use. 0 = all available; negative values subtract from available cores"
        )]
        threads: i32,

        #[arg(short = 's', long, default_value_t = 0.97)]
        threshold: f32,

        #[arg(short = 'a', long)]
        airr: bool,

        #[arg(short = 'k', long)]
        keep_flanking: bool,

        #[arg(short = 'r', long)]
        trust4: bool,

        #[arg(short = 'm', long)]
        max_ram_gb: Option<usize>,
    },
}

#[cfg(feature = "base_cmd")]
pub fn handle_command(cmd: Commands) -> anyhow::Result<()> {
    match cmd {
        Commands::TcrMatch {
            input_file,
            database,
            threads,
            threshold,
            airr,
            keep_flanking,
            trust4,
            max_ram_gb,
        } => {
            let threads = crate::cores::process_cores(Some(threads))?;
            let stdout = std::io::stdout();
            let mut handle = stdout.lock();

            #[cfg(feature = "base_cmd")]
            helper::run_tcrmatch(
                &input_file,
                &database,
                threads,
                threshold,
                airr,
                keep_flanking,
                trust4,
                max_ram_gb,
                &mut handle,
            )
        }
    }
}
