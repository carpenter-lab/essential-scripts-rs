use clap::{Parser, Subcommand};

mod aggregate;
mod copy_cellranger_outs;
mod io;
mod plate_reader;
mod split;
mod tcr_align;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
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
}

pub(crate) fn _main(cli: Cli) {
    match cli.command {
        Some(Commands::Aggregate(cmd)) => aggregate::handle_command(cmd),
        Some(Commands::Split(cmd)) => split::handle_command(cmd),
        Some(Commands::ReformatPlateReaderData(cmd)) => plate_reader::handle_command(cmd),
        Some(Commands::CopyCellRangerOuts(cmd)) => copy_cellranger_outs::handle_command(cmd),
        Some(Commands::TcrAlign(cmd)) => tcr_align::handle_command(cmd),
        None => {}
    }
}

fn main() {
    _main(Cli::parse());
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    #[test]
    fn parses_no_command() {
        let cli = Cli::try_parse_from(["essential-scripts-rs"]).unwrap();
        assert!(cli.command.is_none());
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
