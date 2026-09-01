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
    let build = |paths: Vec<PathBuf>| {
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
    };

    if !*parallel || jobs == &1 {
        return build(paths);
    }

    match rayon::ThreadPoolBuilder::new().num_threads(*jobs).build() {
        Ok(pool) => pool.install(|| build(paths)),
        Err(e) => {
            eprintln!("Failed to build rayon pool: {}; falling back", e);
            build_records_from_paths(paths, pb, &false, &1)
        }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geo_submission::traits::{FromPathWithMd5, HasPath, Md5Record};
    use indicatif::ProgressBar;
    use pretty_assertions::assert_eq;
    use rstest::rstest;
    use std::path::{Path, PathBuf};
    use std::sync::Arc;
    use tempfile::tempdir;

    struct MockItem {
        path: PathBuf,
        md5: String,
    }

    impl HasPath for MockItem {
        fn path(&self) -> &Path {
            &self.path
        }
    }

    impl Md5Record for MockItem {
        fn filename(&self) -> String {
            prepare_paths_report(self)
        }
        fn md5_str(&self) -> &str {
            &self.md5
        }
    }

    impl FromPathWithMd5 for MockItem {
        fn from_path_with_md5(
            path: PathBuf,
            _pb: Option<&Arc<ProgressBar>>,
        ) -> Result<Self, Box<dyn std::error::Error>> {
            if path.to_str().map(|s| s.contains("fail")).unwrap_or(false) {
                return Err("mock failure".into());
            }
            Ok(MockItem {
                path,
                md5: "mock_md5".to_string(),
            })
        }
    }

    #[rstest]
    #[case("file.txt", "file.txt")]
    #[case("/path/to/file.txt", "file.txt")]
    #[case("relative/path/file.txt", "file.txt")]
    #[case("/", "/")]
    #[case(".", ".")]
    #[case("", "")]
    fn test_prepare_paths_report(#[case] path_str: &str, #[case] expected: &str) {
        let item = MockItem {
            path: PathBuf::from(path_str),
            md5: String::new(),
        };
        assert_eq!(prepare_paths_report(&item), expected);
    }

    #[test]
    fn test_generate_md5_report() {
        let items = vec![
            MockItem {
                path: PathBuf::from("file1.txt"),
                md5: "hash1".to_string(),
            },
            MockItem {
                path: PathBuf::from("file2.txt"),
                md5: "hash2".to_string(),
            },
        ];

        let report = generate_md5_report(&items);
        assert!(report.contains("# MD5 Checksums for Files"));
        assert!(report.contains("file1.txt\thash1"));
        assert!(report.contains("file2.txt\thash2"));
    }

    #[test]
    fn test_write_output_to_file() -> Result<(), Box<dyn std::error::Error>> {
        let dir = tempdir()?;
        let file_path = dir.path().join("output.txt");
        let content = "hello world";

        write_output(content, Some(&file_path))?;

        let read_content = fs::read_to_string(file_path)?;
        assert_eq!(read_content, content);
        Ok(())
    }

    #[test]
    fn test_write_output_to_stdout() -> Result<(), Box<dyn std::error::Error>> {
        let content = "hello world";
        // This just checks that it doesn't error out
        write_output(content, None)?;
        Ok(())
    }

    #[rstest]
    #[case(false, 1)]
    #[case(true, 1)]
    #[case(true, 2)]
    fn test_build_records_from_paths(#[case] parallel: bool, #[case] jobs: usize) {
        let paths = vec![
            PathBuf::from("file1.txt"),
            PathBuf::from("file2.txt"),
            PathBuf::from("file3.txt"),
        ];

        let results: Vec<MockItem> = build_records_from_paths(paths, None, &parallel, &jobs);

        assert_eq!(results.len(), 3);
        // results may be out of order due to parallel processing
        let mut paths_found: Vec<String> = results
            .into_iter()
            .map(|i| i.path.to_str().unwrap().to_string())
            .collect();
        paths_found.sort();

        assert_eq!(paths_found, vec!["file1.txt", "file2.txt", "file3.txt"]);
    }

    #[test]
    fn test_build_records_from_paths_with_failure() {
        let paths = vec![
            PathBuf::from("file1.txt"),
            PathBuf::from("fail.txt"),
            PathBuf::from("file2.txt"),
        ];

        // Test parallel with pool
        let results: Vec<MockItem> = build_records_from_paths(paths.clone(), None, &true, &2);
        assert_eq!(results.len(), 2);

        // Test serial (no pool)
        let results: Vec<MockItem> = build_records_from_paths(paths, None, &false, &1);
        assert_eq!(results.len(), 2);
    }
}
