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
