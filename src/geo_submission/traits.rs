use indicatif::ProgressBar;
use std::path::{Path, PathBuf};
use std::sync::Arc;

pub(super) trait Md5Record {
    fn filename(&self) -> String;
    fn md5_str(&self) -> &str;
}

pub(super) trait HasPath {
    fn path(&self) -> &Path;
}

pub(super) trait FromPathWithMd5: Sized + Send {
    /// Build Self from a `PathBuf`, computing MD5 (use pb for progress). Return Err on parse/io failure.
    fn from_path_with_md5(
        path: PathBuf,
        pb: Option<&Arc<ProgressBar>>,
    ) -> Result<Self, Box<dyn std::error::Error>>;
}
