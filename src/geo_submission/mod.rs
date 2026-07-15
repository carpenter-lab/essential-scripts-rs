mod build_reports;
mod fastq;
mod md5;
mod traits;

#[cfg(test)]
mod test;

use clap::{Args, Subcommand};
use clap_binary_enum::YesNoArg;
use fastq::*;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use num_cpus;
use std::path::PathBuf;
use std::sync::Arc;

#[derive(Debug, Clone, Copy, PartialEq, YesNoArg)]
pub enum Progress {
    #[yesno(help = "Show progress bar")]
    Progress,
    #[yesno(help = "Hide progress bar")]
    NoProgress,
}

#[derive(Args, Debug)]
pub struct RunArgs {
    #[arg(long, help = "Compute MD5s in parallel", default_value_t = false)]
    parallel_md5: bool,

    #[arg(
            long,
            value_parser = clap::value_parser!(i32).range(1 + -1 * num_cpus::get() as i64..=num_cpus::get() as i64),
            default_value_t = num_cpus::get().try_into().unwrap(),
            help = "Number of threads to use for parallel MD5"
    )]
    threads: i32,

    #[command(flatten)]
    progress: ProgressArg,
}

#[derive(Subcommand)]
pub enum Commands {
    #[command(about = "Match FastQ files by lane and sample, and compute MD5 checksums")]
    GeoFastq {
        #[arg(help = "Directories containing fastq.gz files", num_args = 1..)]
        input_directories: Vec<PathBuf>,

        #[arg(
            long,
            help = "Output file for paired files by lane and sample [default: stdout]"
        )]
        paired_output: Option<PathBuf>,

        #[arg(long, help = "Output file for all files per sample [default: stdout]")]
        sample_output: Option<PathBuf>,

        #[arg(
            long,
            help = "Output file for file paths with MD5 checksums [default: stdout]"
        )]
        md5_output: Option<PathBuf>,

        #[command(flatten)]
        generate_args: RunArgs,
    },
}

fn make_progress_bar(
    total_bytes: u64,
    progress: Progress,
) -> Result<Arc<ProgressBar>, Box<dyn std::error::Error>> {
    let pb = ProgressBar::new(total_bytes);
    pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta}) - {msg}",
        )?
            .progress_chars("=> "),
    );
    if let Progress::NoProgress = progress {
        pb.set_draw_target(ProgressDrawTarget::hidden());
    }
    Ok(Arc::from(pb))
}

const TOO_MANY_CORES_SUBTRACTED_ERROR: &str = "cannot subtract more cores than available";

fn process_cores(requested_cores: Option<i32>) -> Result<usize, String> {
    let available_cores = num_cpus::get();

    match requested_cores {
        None | Some(0) => Ok(available_cores),
        #[allow(clippy::cast_sign_loss)]
        Some(requested) if requested > 0 => Ok(requested as usize),
        Some(requested) => subtract_from_available_cores(available_cores, requested),
    }
}

fn subtract_from_available_cores(available_cores: usize, requested: i32) -> Result<usize, String> {
    let cores_to_subtract = requested.unsigned_abs() as usize;

    if cores_to_subtract > available_cores {
        Err(TOO_MANY_CORES_SUBTRACTED_ERROR.to_string())
    } else {
        Ok(available_cores - cores_to_subtract)
    }
}

pub fn handle_command(cmd: &Commands) {
    match cmd {
        Commands::GeoFastq {
            input_directories,
            paired_output,
            sample_output,
            md5_output,
            generate_args,
        } => {
            let jobs = process_cores(Some(generate_args.threads)).unwrap_or_else(|e| {
                eprintln!("Could not set number of threads: {e}");
                std::process::exit(1);
            });
            match_fastq(
                input_directories,
                Option::from(paired_output),
                Option::from(sample_output),
                Option::from(md5_output),
                &generate_args.parallel_md5,
                &jobs,
                generate_args.progress.get(),
            );
        }
    }
}
