use crate::cores::{thread_default, thread_parser};
use crate::tcr_cluster::blosum::AminoAcidSubMatrix;
use crate::tcr_cluster::model::{Algorithm, ClusterMethod, DistanceMode, Organism, References};
use clap::{Args, Subcommand};
use std::path::PathBuf;

/// Arguments shared by every clustering algorithm.
#[derive(Args)]
pub struct CommonArgs {
    #[arg(help = "TSV containing subject:condition, TcRb, and optional TcRa/V/J/counts columns")]
    pub chains_file: PathBuf,

    #[arg(help = "Output TSV file path; use '-' for stdout")]
    pub output_file: PathBuf,

    #[arg(
        long,
        help = "TSV with a header row, keyed by subject, followed by left/right allele columns for DPA1, DPB1, DQA1, DQB1, DRB1 and DRB3/4/5"
    )]
    pub hla_file: Option<PathBuf>,

    #[arg(long, value_enum, default_value_t = Default::default(), help = "GLIPH algorithm to use"
    )]
    pub method: Algorithm,

    #[arg(
        long,
        default_value_t = 1000,
        help = "Number of simulations to run for enrichment scores. Minimum score is 1/sim-depth."
    )]
    pub sim_depth: usize,

    #[arg(
        long,
        default_value_t = 0.05,
        help = "Mark (with '!') HLA alleles whose cluster enrichment p-value is at or below this threshold"
    )]
    pub hla_threshold: f64,

    #[arg(long, value_enum, default_value_t = Default::default(), help = "Reference dataset to use for global similarity scoring"
    )]
    pub reference: References,

    #[arg(long)]
    pub reference_cache_dir: Option<PathBuf>,

    #[arg(long, action = clap::ArgAction::SetTrue)]
    pub refresh_reference_cache: bool,

    #[arg(
        long,
        value_parser = thread_parser(),
        default_value_t = thread_default(),
        help = "Number of cores to use. 0 = all available; negative values subtract from available cores"
    )]
    pub cores: i32,

    #[arg(
        long,
        value_enum,
        help = "Substitution matrix to use for creating groups"
    )]
    pub sub_matrix: Option<AminoAcidSubMatrix>,
}

#[deny(clippy::missing_docs_in_private_items)]
#[derive(Subcommand)]
pub enum Commands {
    /// A reimplementation of the GLIPH algorithms
    ///
    /// `run-gliph` reads a chain table with required `subject:condition` and `TcRb`
    /// columns (and optional `TcRa`, `V`, `J`, and `counts` columns). Each
    /// `subject:condition` value is split at the first colon during processing. An
    /// optional HLA TSV can be supplied with `--hla-file`; it has a header row and
    /// each row contains the subject followed by left/right allele columns for
    /// `HLA-DPA1`, `HLA-DPB1`, `HLA-DQA1`, `HLA-DQB1`, `HLA-DRB1`, and `HLA-DRB3/4/5`
    /// (`NA` or empty means missing). Output has eight HLA columns (DRB3, DRB4 and DRB5
    /// are split by allele name); distinct alleles are joined with `/`, e.g.
    /// `DPA1*01:03/DPA1*02:01`, and all of them are used for scoring. Alleles whose enrichment p-value within a cluster is at
    /// or below
    /// `--hla-threshold` (default 0.05) are suffixed with `!` in the output. When omitted,
    /// HLA output columns are left blank. It writes one row per chain that belongs
    /// to a GLIPH cluster.
    /// The output is a TSV with the fixed column order
    /// `index`, `pattern`, `number_subject`, `number_unique_cdr3`, `TcRb`,
    /// `Fisher_score`, `Freq`, `subject`, `condition`, `final_score`, `hla_score`,
    /// `vb_score`, `expansion_score`, `length_score`, `cluster_size_score`, `type`,
    /// `ulTcRb`, `V`, `J`, `TcRa`, and the eight HLA columns above. When the input has
    /// CDR3α sequences, a `pattern_alpha` column follows with the motif shared by each
    /// cluster's CDR3αs (found as for `run-tcr-dist`), then `VA` and `JA` when alpha genes
    /// are given.
    RunGliph {
        #[command(flatten)]
        args: CommonArgs,
    },
    /// Cluster chains using tcrdist distances and score the clusters
    ///
    /// Output matches `run-gliph`: `pattern` is the CDR3β motif of each cluster and, when
    /// the input has CDR3α sequences, `pattern_alpha` is the CDR3α motif.
    RunTcrDist {
        #[command(flatten)]
        args: CommonArgs,

        #[command(flatten)]
        tcrdist: TcrDistArgs,
    },
}

/// Arguments specific to `run-tcr-dist`.
#[derive(Args)]
pub struct TcrDistArgs {
    #[arg(
        long,
        value_enum,
        help = "Chains and CDRs used for the distance. V genes without an allele are looked up as `*01` [default: alpha-beta if CDR3α is present, else beta]"
    )]
    pub distance_mode: Option<DistanceMode>,

    #[arg(
        long,
        value_enum,
        default_value_t = Default::default(),
        help = "Organism of the germline V genes for the *-full distance modes"
    )]
    pub organism: Organism,

    #[arg(
        long,
        value_enum,
        default_value_t = Default::default(),
        help = "How records within each CDR3-edit threshold are grouped into clusters"
    )]
    pub cluster_method: ClusterMethod,

    #[arg(
        long,
        default_value_t = crate::tcr_cluster::tcrdist::DEFAULT_MIN_NEIGHBOR_OVERLAP,
        value_parser = parse_fraction,
        help = "For shared-neighbor clustering: link two records within the threshold only if the Jaccard overlap of their neighborhoods (each record plus its neighbors) is at least this. 0 is single linkage; higher values split chains more readily"
    )]
    pub min_neighbor_overlap: f32,
}

fn parse_fraction(value: &str) -> Result<f32, String> {
    let x: f32 = value.parse().map_err(|e| format!("{e}"))?;
    if (0.0..=1.0).contains(&x) {
        Ok(x)
    } else {
        Err("must be between 0 and 1".to_owned())
    }
}

impl CommonArgs {
    pub fn validate(&self) -> anyhow::Result<()> {
        if self.sim_depth == 0 {
            return Err(clap::Error::raw(
                clap::error::ErrorKind::ValueValidation,
                "--sim-depth must be greater than zero",
            )
            .into());
        }
        if !(0.0..=1.0).contains(&self.hla_threshold) {
            return Err(clap::Error::raw(
                clap::error::ErrorKind::ValueValidation,
                "--hla-threshold must be between 0 and 1",
            )
            .into());
        }
        Ok(())
    }
}
