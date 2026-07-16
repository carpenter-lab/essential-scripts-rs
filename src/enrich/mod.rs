#[cfg(feature = "enrichment")]
mod api;
#[cfg(feature = "enrichment")]
mod core;

use clap::{Error, Subcommand, ValueEnum};
use std::fmt;
use std::path::PathBuf;
use tokio;

#[derive(Clone, ValueEnum, Debug)]
pub enum Library {
    ReactomePathways2024,
    Reactome,
    BioCarta2016,
    WikiPathways2024Human,
    GOBiologicalProcess,
}

impl fmt::Display for Library {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s = match self {
            Library::ReactomePathways2024 => String::from("Reactome_Pathways_2024"),
            Library::Reactome => Library::ReactomePathways2024.to_string(),
            Library::BioCarta2016 => String::from("BioCarta_2016"),
            Library::WikiPathways2024Human => String::from("WikiPathways_2024_Human"),
            Library::GOBiologicalProcess => String::from("GO_Biological_Process_2026"),
        };
        write!(f, "{s}")
    }
}

#[cfg(feature = "enrichment")]
fn must_be_none(pb: Option<&PathBuf>) -> Result<Option<PathBuf>, String> {
    match pb {
        None => Ok(None),
        Some(_) => Err("Background is not supported in the API".to_string()),
    }
}

#[derive(Subcommand)]
pub enum Commands {
    #[command(about = "Run Enrichr via API interface.")]
    RunEnrichr {
        #[arg(short, long, help = "Enrichr Library to use")]
        library: Library,

        #[arg(
            short,
            long,
            required = true,
            help = "Input gene list to process. One gene per line"
        )]
        gene_list: PathBuf,

        #[arg(
            short,
            long,
            help = "Input gene list to process as background. One gene per line"
        )]
        background: Option<PathBuf>,

        #[arg(help = "Output file path")]
        output_file: PathBuf,

        #[arg(num_args = 1.., help = "Output file paths")]
        output_plot: Vec<PathBuf>,
    },
}

#[tokio::main]
pub async fn handle_command(cmd: Commands) -> Result<(), Error> {
    #[cfg(feature = "enrichment")]
    {
        match cmd {
            Commands::RunEnrichr {
                library,
                gene_list,
                background,
                output_file,
                output_plot,
            } => {
                if let Err(e) = must_be_none(background.as_ref()) {
                    return Err(Error::raw(clap::error::ErrorKind::ValueValidation, e));
                }
                match core::enrich_command(
                    library.to_string(),
                    gene_list,
                    background,
                    output_file,
                    output_plot,
                )
                .await
                {
                    Ok(()) => Ok(()),
                    Err(e) => Err(Error::raw(
                        clap::error::ErrorKind::Io,
                        format!("Enrichment analysis failed: {e}"),
                    )),
                }
            }
        }
    }

    #[cfg(not(feature = "enrichment"))]
    {
        match cmd {
            Commands::RunEnrichr { .. } => Err(Error::raw(
                clap::error::ErrorKind::MissingSubcommand,
                "This command requires the `enrichment` feature. Rebuild with `--features enrichment`",
            )),
        }
    }
}
