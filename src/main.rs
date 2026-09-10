use clap::CommandFactory;
use clap::{Parser, Subcommand};
mod docs;

use essential_scripts_rs::{
    aggregate, copy_cellranger_outs, dryad, enrich, geo_submission, plate_reader, split, tcr_align,
    underlying_clap_error_kind,
};

#[derive(Parser)]
#[command(version, about, long_about = None, arg_required_else_help = true)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,

    #[arg(long, hide = true)]
    markdown_help: bool,
}

#[derive(Subcommand)]
enum Commands {
    #[command(flatten)]
    Aggregate(aggregate::Commands),

    #[command(flatten)]
    Split(split::Commands),

    #[command(flatten)]
    ReformatPlateReaderData(plate_reader::Commands),

    #[command(flatten)]
    CopyCellRangerOuts(copy_cellranger_outs::Commands),

    #[command(flatten)]
    TcrAlign(tcr_align::Commands),

    #[command(flatten)]
    MatchFastq(geo_submission::Commands),

    #[command(flatten)]
    RunEnrichr(enrich::Commands),

    #[command(flatten)]
    Dryad(dryad::Commands),
}

#[cfg(feature = "base_cmd")]
pub(crate) fn main_helper(cli: Cli) {
    let result = match cli.command {
        Some(Commands::Aggregate(cmd)) => aggregate::handle_command(cmd),
        Some(Commands::Split(cmd)) => split::handle_command(cmd),
        Some(Commands::ReformatPlateReaderData(cmd)) => plate_reader::handle_command(cmd),
        Some(Commands::CopyCellRangerOuts(cmd)) => copy_cellranger_outs::handle_command(cmd),
        Some(Commands::TcrAlign(cmd)) => tcr_align::handle_command(cmd),
        Some(Commands::MatchFastq(cmd)) => geo_submission::handle_command(&cmd),
        Some(Commands::RunEnrichr(cmd)) => enrich::handle_command(cmd),
        Some(Commands::Dryad(cmd)) => dryad::handle_command(cmd),
        None => Ok(()),
    };
    if let Err(err) = result {
        match underlying_clap_error_kind(&err) {
            Some(clap::error::ErrorKind::MissingSubcommand) => {
                let mut c = Cli::command();
                c.error(
                    clap::error::ErrorKind::MissingSubcommand,
                    err.to_string().replace("error: ", ""),
                )
                .exit();
            }
            _ => {
                eprintln!("{err}");
                std::process::exit(1)
            }
        }
    }
}

#[cfg(not(feature = "base_cmd"))]
pub(crate) fn main_helper(_cli: Cli) {
    println!(
        "Please enable the `base_cmd` feature flag. \
        This should be the default and using --no-default-features \
        is meant only for document generation."
    );
}

fn main() {
    let args = Cli::try_parse().unwrap_or_else(|error| error.exit());

    if args.markdown_help {
        match docs::write_docs_to_file::<Cli>("docs/cli.md") {
            Ok(_) => {
                println!("Markdown help written to docs/cli.md");
                std::process::exit(0);
            }
            Err(err) => {
                eprintln!("Failed to write markdown help: {err}");
                std::process::exit(1);
            }
        }
    }

    main_helper(args);
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use pretty_assertions::assert_eq;

    #[test]
    fn parses_no_command() {
        let cli = Cli::try_parse_from(["essential-scripts-rs"]);
        assert_eq!(
            cli.err().unwrap().kind(),
            clap::error::ErrorKind::DisplayHelpOnMissingArgumentOrSubcommand,
        );
    }

    #[test]
    fn rejects_unknown_command() {
        if let Err(err) = Cli::try_parse_from(["essential-scripts-rs", "does-not-exist"]) {
            assert_eq!(err.kind(), clap::error::ErrorKind::InvalidSubcommand);
        } else {
            panic!("expected error")
        }
    }
}
