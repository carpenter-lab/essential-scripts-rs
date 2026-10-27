#[cfg(feature = "base_cmd")]
mod analysis;
mod blosum;
#[cfg(feature = "base_cmd")]
mod input;
#[cfg(feature = "base_cmd")]
mod output;
#[cfg(feature = "base_cmd")]
mod reference;
pub mod score;
pub mod stats;
mod types;

pub use types::*;

use crate::cores::{thread_default, thread_parser};
use clap::Subcommand;
use std::path::PathBuf;

use crate::gliph::blosum::AminoAcidSubMatrix;
#[cfg(feature = "base_cmd")]
use anyhow::{Context, Result};

#[derive(Subcommand)]
pub enum Commands {
    /// A reimplementation of the GLIPH algorithms
    RunGliph {
        #[arg(
            help = "TSV containing subject:condition, TcRb, and optional TcRa/V/J/counts columns"
        )]
        chains_file: PathBuf,

        #[arg(help = "Output TSV file path; use '-' for stdout")]
        output_file: PathBuf,

        #[arg(
            long,
            help = "TSV with a header row, keyed by subject, followed by left/right allele columns for DPA1, DPB1, DQA1, DQB1, DRB1 and DRB3/4/5"
        )]
        hla_file: Option<PathBuf>,

        #[arg(long, value_enum, default_value_t = Default::default(), help = "GLIPH algorithm to use"
        )]
        method: Algorithm,

        #[arg(
            long,
            default_value_t = 1000,
            help = "Number of simulations to run for enrichment scores. Minimum score is 1/sim-depth."
        )]
        sim_depth: usize,

        #[arg(
            long,
            default_value_t = 0.05,
            help = "Mark (with '!') HLA alleles whose cluster enrichment p-value is at or below this threshold"
        )]
        hla_threshold: f64,

        #[arg(long, value_enum, default_value_t = Default::default(), help = "Reference dataset to use for global similarity scoring"
        )]
        reference: References,

        #[arg(long)]
        reference_cache_dir: Option<PathBuf>,

        #[arg(long, action = clap::ArgAction::SetTrue)]
        refresh_reference_cache: bool,

        #[arg(
            long,
            value_parser = thread_parser(),
            default_value_t = thread_default(),
            help = "Number of cores to use. 0 = all available; negative values subtract from available cores"
        )]
        cores: i32,

        #[arg(
            long,
            value_enum,
            help = "Substitution matrix to use for creating groups"
        )]
        sub_matrix: Option<AminoAcidSubMatrix>,
    },
}

#[cfg(feature = "base_cmd")]
pub fn handle_command(command: Commands) -> Result<()> {
    match command {
        Commands::RunGliph {
            chains_file,
            output_file,
            hla_file,
            method,
            sim_depth,
            hla_threshold,
            reference,
            reference_cache_dir,
            refresh_reference_cache,
            cores,
            sub_matrix,
        } => {
            let cores = crate::cores::process_cores(Some(cores))?;
            if sim_depth == 0 {
                return Err(clap::Error::raw(
                    clap::error::ErrorKind::ValueValidation,
                    "--sim-depth must be greater than zero",
                )
                .into());
            }
            if !(0.0..=1.0).contains(&hla_threshold) {
                return Err(clap::Error::raw(
                    clap::error::ErrorKind::ValueValidation,
                    "--hla-threshold must be between 0 and 1",
                )
                .into());
            }
            let records = input::read_chain_inputs(&chains_file, hla_file.as_deref())?;
            let reference = reference::load_reference(
                &reference,
                reference_cache_dir.as_deref(),
                refresh_reference_cache,
            )?;
            let hla_context = HlaContext::build(&records, hla_threshold);
            let expansion_ctx = ExpansionContext::build(&records);
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(cores)
                .build()
                .context("Failed to build thread pool")?;
            let rows = pool
                .install(|| {
                    analysis::analyze(
                        &records,
                        &reference,
                        method,
                        sim_depth,
                        &hla_context,
                        &expansion_ctx,
                        sub_matrix,
                    )
                })
                .context("GLIPH analysis failed")?;
            output::write_output(&output_file, &rows)?;
            Ok(())
        }
    }
}

#[cfg(not(feature = "base_cmd"))]
pub fn handle_command(_command: Commands) -> anyhow::Result<()> {
    Err(anyhow::Error::new(clap::Error::raw(
        clap::error::ErrorKind::MissingSubcommand,
        "This command requires the `base_cmd` feature. Rebuild with `--features base_cmd`",
    )))
}
