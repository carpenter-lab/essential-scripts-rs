use anyhow::Result;
use clap_binary_enum::YesNoArg;
use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressStyle};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

#[derive(Debug, Clone, Copy, PartialEq, YesNoArg)]
pub enum Progress {
    #[yesno(help = "Show progress bar")]
    Progress,
    #[yesno(help = "Hide progress bar")]
    NoProgress,
}

trait SetDrawTarget<T> {
    fn apply(&self, pb: &T);
}

impl SetDrawTarget<ProgressBar> for Progress {
    fn apply(&self, pb: &ProgressBar) {
        if matches!(self, Self::NoProgress) {
            pb.set_draw_target(ProgressDrawTarget::hidden());
        }
    }
}
impl SetDrawTarget<MultiProgress> for Progress {
    fn apply(&self, pb: &MultiProgress) {
        if matches!(self, Self::NoProgress) {
            pb.set_draw_target(ProgressDrawTarget::hidden());
        }
    }
}
impl SetDrawTarget<FileByteProgress> for Progress {
    fn apply(&self, pb: &FileByteProgress) {
        if matches!(self, Self::NoProgress) {
            pb.pb.set_draw_target(ProgressDrawTarget::hidden());
        }
    }
}

const BYTE_PROGRESS_TEMPLATE: &str =
    "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {msg} ({eta})";
const OUTER_PROGRESS_TEMPLATE: &str =
    "{msg} [{bar:40.cyan/blue}] {pos}/{len} Elapsed: {elapsed_precise} ETA: {eta}";
const INNER_PROGRESS_TEMPLATE: &str = "\x1b[37m{msg}\x1b[0m [{bar:40.cyan/blue}] {pos}/{len}";

pub struct FileByteProgress {
    pb: ProgressBar,
    total_files: usize,
    files_done: AtomicUsize,
}

impl FileByteProgress {
    pub fn new(total_bytes: u64, total_files: usize) -> Result<FileByteProgress> {
        let pb = ProgressBar::new(total_bytes);
        pb.set_style(ProgressStyle::with_template(BYTE_PROGRESS_TEMPLATE)?.progress_chars("=> "));
        pb.set_message(format!("files 0/{}", total_files));
        Ok(FileByteProgress {
            pb,
            total_files,
            files_done: AtomicUsize::new(0),
        })
    }
    pub fn inc_bytes(&self, n: u64) {
        self.pb.inc(n);
    }

    pub fn inc_file(&self) {
        let done = self.files_done.fetch_add(1, Ordering::Relaxed) + 1;
        self.pb
            .set_message(format!("files {done}/{}", self.total_files));
    }

    pub fn position(&self) -> u64 {
        self.pb.position()
    }
    pub fn length(&self) -> Option<u64> {
        self.pb.length()
    }
}

pub fn byte_and_file_progress_bar(
    total_bytes: u64,
    total_files: usize,
    visibility: Progress,
) -> Result<Arc<FileByteProgress>> {
    let pb = FileByteProgress::new(total_bytes, total_files)?;
    visibility.apply(&pb);

    Ok(Arc::new(pb))
}

pub fn dual_progress_bars(
    outer_len: u64,
    inner_len: u64,
    visibility: Progress,
) -> Result<(ProgressBar, ProgressBar)> {
    let outer_style = ProgressStyle::with_template(OUTER_PROGRESS_TEMPLATE)?;
    let inner_style = ProgressStyle::with_template(INNER_PROGRESS_TEMPLATE)?;

    let mp = MultiProgress::new();
    let outer = mp.add(ProgressBar::new(outer_len));
    outer.set_style(outer_style);

    let inner = mp.add(ProgressBar::new(inner_len));
    inner.set_style(inner_style);
    visibility.apply(&mp);

    Ok((outer, inner))
}

#[cfg(test)]
mod test {
    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    #[rstest]
    #[case::progress_bar(Progress::Progress, 1000, 1000)]
    #[case::progress_bar_hidden(Progress::NoProgress, 1000, 1000)]
    #[case::progress_bar_zero_bytes(Progress::Progress, 0, 0)]
    fn test_byte_and_file_progress_bar(
        #[case] progress: Progress,
        #[case] total_bytes: u64,
        #[case] expected_length: u64,
    ) {
        let result = byte_and_file_progress_bar(total_bytes, 1, progress);
        assert!(result.is_ok());

        let pb = result.unwrap();
        assert_eq!(pb.length(), Some(expected_length));
    }

    #[rstest]
    #[case::progress(Progress::Progress)]
    #[case::no_progress(Progress::NoProgress)]
    fn test_dual_progress_bars(#[case] progress: Progress) {
        let (outer, inner) = dual_progress_bars(1000, 1000, progress).unwrap();
        assert_eq!(outer.length(), Some(1000));
        assert_eq!(inner.length(), Some(1000));
    }
}
