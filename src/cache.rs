use anyhow::{Context, bail};
use std::fs;
use std::fs::OpenOptions;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Duration;

fn default_cache_dir() -> anyhow::Result<PathBuf> {
    if let Some(path) = std::env::var_os("XDG_CACHE_HOME") {
        return Ok(PathBuf::from(path)
            .join("essential-scripts-rs")
            .join("gliph"));
    }
    let home = std::env::var_os("HOME").context("HOME is not set; pass --reference-cache-dir")?;
    #[cfg(target_os = "macos")]
    let base = PathBuf::from(home).join("Library").join("Caches");
    #[cfg(not(target_os = "macos"))]
    let base = PathBuf::from(home).join(".cache");
    Ok(base.join("essential-scripts-rs").join("gliph"))
}

pub fn cache_path(cache_dir: Option<&Path>, filename: &str) -> anyhow::Result<PathBuf> {
    let root = cache_dir
        .map(Path::to_path_buf)
        .map_or_else(default_cache_dir, Ok)?;
    fs::create_dir_all(&root).with_context(|| {
        format!(
            "Failed to create reference cache directory {}",
            root.display()
        )
    })?;
    Ok(root.join(filename))
}

pub fn download_file(
    url: &str,
    path: &Path,
    verify_reference: bool,
    file_size: usize,
    md5: &str,
) -> anyhow::Result<()> {
    let client = reqwest::blocking::Client::builder()
        .user_agent(concat!("essential-scripts-rs/", env!("CARGO_PKG_VERSION")))
        .timeout(Duration::from_secs(300))
        .build()
        .context("Failed to create HTTP client for reference data")?;
    let mut response = client
        .get(url)
        .send()
        .with_context(|| format!("Failed to download reference data from {url}"))?
        .error_for_status()
        .with_context(|| format!("Reference server returned an error for {url}"))?;
    let temporary = path.with_extension(format!("download-{}", std::process::id()));
    let mut output = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temporary)
        .with_context(|| {
            format!(
                "Failed to create temporary cache file {}",
                temporary.display()
            )
        })?;
    std::io::copy(&mut response, &mut output)
        .with_context(|| format!("Failed while saving reference data from {url}"))?;
    output.flush()?;
    drop(output);

    if verify_reference {
        let bytes = fs::read(&temporary).with_context(|| {
            format!("Failed to validate downloaded file {}", temporary.display())
        })?;
        let checksum = format!("{:x}", md5::compute(&bytes));
        if bytes.len() != file_size || checksum != md5 {
            fs::remove_file(&temporary).ok();
            bail!(
                "Downloaded GLIPH reference data failed integrity validation (size {}, MD5 {checksum})",
                bytes.len()
            );
        }
    }

    fs::rename(&temporary, path).with_context(|| {
        format!(
            "Failed to move downloaded reference data into cache at {}",
            path.display()
        )
    })?;
    Ok(())
}

pub fn cached_file(
    cache_dir: Option<&Path>,
    filename: &str,
    url: &str,
    force_download: bool,
    verify_reference: bool,
    file_size: usize,
    md5: &str,
) -> anyhow::Result<PathBuf> {
    let path = cache_path(cache_dir, filename)?;
    if force_download && path.exists() {
        fs::remove_file(&path)
            .with_context(|| format!("Failed to remove cached reference {}", path.display()))?;
    }
    if !path.exists() {
        download_file(url, &path, verify_reference, file_size, md5)?;
    }
    Ok(path)
}
