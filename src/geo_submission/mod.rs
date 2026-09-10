#[cfg(feature = "base_cmd")]
mod build_reports;
#[cfg(feature = "base_cmd")]
mod fastq;
#[cfg(feature = "base_cmd")]
mod helper;
#[cfg(feature = "base_cmd")]
mod md5;
#[cfg(feature = "base_cmd")]
mod traits;

#[cfg(test)]
mod test;

use crate::progress::ProgressArg;
use anyhow::Result;
use clap::builder::RangedI64ValueParser;
use clap::{Args, Subcommand};
#[cfg(feature = "base_cmd")]
use fastq::*;
use std::path::PathBuf;

#[cfg(test)]
fn num_cpu() -> usize {
    4
}

#[cfg(not(test))]
fn num_cpu() -> usize {
    num_cpus::get()
}

#[cfg(feature = "base_cmd")]
fn thread_parser() -> RangedI64ValueParser<i32> {
    clap::value_parser!(i32).range(1 + -(num_cpu() as i64)..=num_cpu() as i64)
}
#[cfg(feature = "base_cmd")]
fn thread_default() -> i32 {
    num_cpu().try_into().unwrap()
}

#[cfg(not(feature = "base_cmd"))]
fn thread_default() -> i32 {
    0
}
#[cfg(not(feature = "base_cmd"))]
fn thread_parser() -> RangedI64ValueParser<i32> {
    clap::value_parser!(i32).range(-1..1)
}

#[derive(Args, Debug)]
pub struct RunArgs {
    #[arg(long, help = "Compute MD5s in parallel", default_value_t = false)]
    parallel_md5: bool,

    #[arg(
        long,
        value_parser = thread_parser(),
        default_value_t = thread_default(),
        help = "Number of threads to use for parallel MD5. 0 = all available cores"
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

#[cfg(feature = "base_cmd")]
pub fn handle_command(cmd: &Commands) -> Result<()> {
    match cmd {
        Commands::GeoFastq {
            input_directories,
            paired_output,
            sample_output,
            md5_output,
            generate_args,
        } => {
            let jobs = helper::process_cores(Some(generate_args.threads)).unwrap_or_else(|e| {
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
            )
        }
    }
}
