use crate::geo_submission::traits::*;
use indicatif::ProgressBar;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::fs;
use std::path::PathBuf;
use std::sync::Arc;

pub fn prepare_paths_report<T: HasPath>(item: &T) -> String {
    item.path()
        .file_name()
        .and_then(|n| n.to_str())
        .map_or_else(
            || item.path().to_string_lossy().to_string(),
            |s| s.to_string(),
        )
}

pub(super) fn build_records_from_paths<T: FromPathWithMd5>(
    paths: Vec<PathBuf>,
    pb: Option<&Arc<ProgressBar>>,
    parallel: &bool,
    jobs: &usize,
) -> Vec<T> {
    if *parallel && jobs > &1 {
        let pool = rayon::ThreadPoolBuilder::new().num_threads(*jobs).build();
        match pool {
            Ok(pool) => pool.install(|| {
                paths
                    .into_par_iter()
                    .filter_map(|p| match T::from_path_with_md5(p, pb) {
                        Ok(rec) => Some(rec),
                        Err(e) => {
                            eprintln!("Error building record: {}", e);
                            None
                        }
                    })
                    .collect()
            }),
            Err(e) => {
                eprintln!("Failed to build rayon pool: {}; falling back", e);
                paths
                    .into_par_iter()
                    .filter_map(|p| T::from_path_with_md5(p, pb).ok())
                    .collect()
            }
        }
    } else {
        paths
            .into_par_iter()
            .filter_map(|p| T::from_path_with_md5(p, pb).ok())
            .collect()
    }
}

/// Generate MD5 manifest report
pub fn generate_md5_report<T: Md5Record>(items: &[T]) -> String {
    let mut output = String::new();
    output.push_str("# MD5 Checksums for Files\n");
    output.push_str("# Format: filename\tmd5_hash\n\n");

    for file in items {
        let fname = file.filename();
        output.push_str(&format!("{}\t{}\n", fname, file.md5_str()));
    }

    output
}

/// Write output to file or stdout
pub fn write_output(
    content: &str,
    output_path: Option<&PathBuf>,
) -> Result<(), Box<dyn std::error::Error>> {
    match output_path {
        Some(path) => {
            fs::write(path, content)?;
            println!("Wrote output to: {}", path.display());
        }
        None => {
            print!("{}", content);
        }
    }
    Ok(())
}
