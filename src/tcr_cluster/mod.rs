//! Unified TCR clustering: GLIPH (1/2) and tcrdist share the same input, scoring and output model.
//!
//! Layout:
//! - `model`: shared data types (records, clusters, scores, output rows, scoring context)
//! - `io`: chain/reference readers and the TSV writer
//! - `gliph` / `tcrdist`: algorithm-specific cluster finders implementing [`ClusterFinder`]
//! - `scoring`, `stats`, `blosum`: scoring and statistics building blocks
//! - `pipeline`: find → score → format, shared by all finders

mod blosum;
mod cli;
mod model;
pub mod scoring;
pub mod stats;

mod gliph;
mod io;
mod motif;
mod pipeline;
mod tcrdist;

pub use cli::Commands;
pub use model::*;
pub use pipeline::ClusterFinder;

use anyhow::{Context, Result};
use blosum::AminoAcidSubMatrix;
use cli::CommonArgs;

pub fn handle_command(command: Commands) -> Result<()> {
    match command {
        Commands::RunGliph { args } => {
            let finder = gliph::Gliph {
                algorithm: args.method,
                sim_depth: args.sim_depth,
                sub_matrix: args.sub_matrix.clone(),
            };
            run_command(&finder, &args).context("GLIPH analysis failed")
        }
        Commands::RunTcrDist { args, tcrdist } => {
            let mode = tcrdist.distance_mode;
            // Germline CDR modes are only ever chosen explicitly, so the db is
            // downloaded only when asked for.
            let gene_db = mode
                .filter(|mode| mode.uses_germline_cdrs())
                .map(|_| {
                    tcrdist::TcrDistDb::load(
                        tcrdist.organism,
                        args.reference_cache_dir.as_deref(),
                        args.refresh_reference_cache,
                    )
                })
                .transpose()
                .context("Failed to load the tcrdist3 gene database")?;
            let finder = tcrdist::TcrDist {
                sub_matrix: args
                    .sub_matrix
                    .clone()
                    .unwrap_or(AminoAcidSubMatrix::Blosum62),
                params: DistanceParams::new(mode.unwrap_or(DistanceMode::Beta)),
                infer_mode: mode.is_none(),
                gene_db,
                cluster_method: tcrdist.cluster_method,
                min_neighbor_overlap: tcrdist.min_neighbor_overlap,
            };
            run_command(&finder, &args).context("tcrdist analysis failed")
        }
    }
}

fn run_command<F: ClusterFinder>(finder: &F, args: &CommonArgs) -> Result<()> {
    let cores = crate::cores::process_cores(Some(args.cores))?;
    args.validate()?;
    let records = io::input::read_chain_inputs(&args.chains_file, args.hla_file.as_deref())?;
    //let records: Vec<ChainRecord> = records[0..50000].to_vec();
    let reference = io::reference::load_reference(
        &args.reference,
        args.reference_cache_dir.as_deref(),
        args.refresh_reference_cache,
    )?;
    let params = pipeline::ScoringParams {
        algorithm: args.method,
        sim_depth: args.sim_depth,
        hla_threshold: args.hla_threshold,
    };
    let pool = crate::cores::build_pool(cores)?;
    let rows = pool.install(|| pipeline::run(finder, &records, &reference, params))?;
    io::output::write_output(&args.output_file, &rows)
}
