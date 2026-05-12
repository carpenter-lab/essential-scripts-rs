#[cfg(feature = "enrichment")]
mod api;
#[cfg(feature = "enrichment")]
mod core;

use clap::Subcommand;
#[cfg(feature = "enrichment")]
use std::fs;
use std::path::PathBuf;

#[derive(Subcommand)]
pub enum Commands {
    #[command(about = "Run Enrichr via API interface.")]
    RunEnrichr {
        #[arg(short, long, help = "Enrichr Library to use")]
        library: String,

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

#[cfg(feature = "enrichment")]
async fn enrich_command(
    library: String,
    gene_list: PathBuf,
    background: Option<PathBuf>,
    output_file: PathBuf,
    output_plot: Vec<PathBuf>,
) -> Result<(), Box<dyn std::error::Error>> {
    let genes: Vec<String> = fs::read_to_string(&gene_list)
        .map(|line| line.trim().to_string())
        .into_iter()
        .collect();
    let libraries = vec![library.clone()];
    let mut enrich = core::Enrichment::new(genes, libraries);

    if let Some(path) = background {
        let bg_genes: Vec<String> = fs::read_to_string(&path)
            .map(|line| line.trim().to_string())
            .into_iter()
            .collect();
        enrich.with_background(bg_genes);
    }
    enrich.build();
    enrich.run().await?;
    enrich.save_results(output_file)?;
    enrich
        .bar_plot(output_plot, Some(library), None, None)
        .await?;
    Ok(())
}

#[tokio::main]
pub async fn handle_command(cmd: Commands) -> Result<(), clap::Error> {
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
                match enrich_command(library, gene_list, background, output_file, output_plot).await
                {
                    Ok(_) => Ok(()),
                    Err(e) => Err(clap::Error::raw(
                        clap::error::ErrorKind::Io,
                        format!("Enrichment analysis failed: {}", e),
                    )),
                }
            }
        }
    }

    #[cfg(not(feature = "enrichment"))]
    {
        match cmd {
            Commands::RunEnrichr { .. } => Err(clap::Error::raw(
                clap::error::ErrorKind::MissingSubcommand,
                "This command requires the `enrichment` feature. Rebuild with `cargo run --features enrichment -- ...`",
            )),
        }
    }
}
