use indicatif::ProgressBar;
use md5;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;
use std::sync::Arc;

/// Compute MD5 and increment a progress bar by the number of bytes read
pub(crate) fn compute_md5_with_progress(
    path: &Path,
    pb: Option<&Arc<ProgressBar>>,
) -> Result<String, std::io::Error> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut buffer = vec![0u8; 8 * 1024 * 1024]; // 8 MiB buffer
    let mut ctx = md5::Context::new();
    loop {
        let n = reader.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        ctx.consume(&buffer[..n]);
        if let Some(pb) = pb {
            pb.inc(n as u64);
        }
    }
    let digest = ctx.finalize();
    Ok(format!("{digest:x}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::{Context, fixture, rstest};
    use std::fs::File;
    use std::io::Write;
    use std::path::PathBuf;
    use tempfile::TempDir;

    /// Fixture: Create a temporary file with test data
    #[fixture]
    fn temp_dir() -> TempDir {
        TempDir::new().unwrap()
    }

    fn create_test_file(dir: &TempDir, filename: &str, content: &[u8]) -> PathBuf {
        let path = dir.path().join(filename);
        let mut file = File::create(&path).unwrap();
        file.write_all(content).unwrap();
        path
    }

    #[rstest]
    #[case(b"hello", "5d41402abc4b2a76b9719d911017c592")] // MD5 of "hello"
    #[case(b"", "d41d8cd98f00b204e9800998ecf8427e")] // MD5 of empty string
    #[case(b"test data", "eb733a00c0c9d336e65691a37ab54293")] // MD5 of "test data"
    fn test_compute_md5_without_progress(
        temp_dir: TempDir,
        #[case] content: &[u8],
        #[case] expected_md5: &str,
        #[context] context: Context,
    ) {
        let file = create_test_file(&temp_dir, context.name, content);
        let result = compute_md5_with_progress(&file, None);

        assert!(result.is_ok());
        assert_eq!(result.unwrap(), expected_md5);
    }

    #[rstest]
    fn test_compute_md5_with_progress_bar(temp_dir: TempDir, #[context] context: Context) {
        let content = b"test content";
        let path = create_test_file(&temp_dir, context.name, content);

        let pb = Arc::new(ProgressBar::new(content.len() as u64));
        let result = compute_md5_with_progress(&path, Some(&pb));

        assert!(result.is_ok());
        assert_eq!(pb.position(), content.len() as u64);
    }

    #[test]
    fn test_compute_md5_nonexistent_file() {
        let path = Path::new("/nonexistent/file");
        let result = compute_md5_with_progress(path, None);

        assert!(result.is_err());
    }

    #[rstest]
    fn test_compute_md5_large_file(temp_dir: TempDir, #[context] context: Context) {
        // Test with large file (> buffer size)
        let content = vec![0u8; 16 * 1024 * 1024]; // 16 MiB
        let path = create_test_file(&temp_dir, context.name, &content);

        let pb = Arc::new(ProgressBar::new(content.len() as u64));
        let result = compute_md5_with_progress(&path, Some(&pb));

        assert!(result.is_ok());
        assert_eq!(pb.position(), content.len() as u64);
    }
}
