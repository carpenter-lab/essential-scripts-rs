use anyhow::Result;
use clap_binary_enum::YesNoArg;
use indicatif::{
    HumanCount, MultiProgress, ProgressBar, ProgressDrawTarget, ProgressState, ProgressStyle,
};
use std::cmp::PartialEq;
use std::str::FromStr;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Duration;

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

const BYTE_PROGRESS_TEMPLATE: &str = "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} {msg} ETA: {eta}";

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
    inner_len: Option<u64>,
    outer_settings: ProgressSettings,
    inner_settings: ProgressSettings,
    visibility: Progress,
) -> Result<(ProgressBar, ProgressBar)> {
    let mp = MultiProgress::new();
    let outer = mp.add(ProgressBar::new(outer_len));

    let template = progress_template(outer_settings, false)?;
    set_progress_style(&outer, &template)?;
    outer.enable_steady_tick(Duration::from_secs(1));

    let inner = if let Some(len) = inner_len {
        mp.add(ProgressBar::new(len))
    } else {
        mp.add(ProgressBar::no_length())
    };

    let template = progress_template(inner_settings, true)?;
    set_progress_style(&inner, &template)?;
    visibility.apply(&mp);

    Ok((outer, inner))
}

#[derive(PartialEq, Eq, Clone, Copy, Debug)]
pub enum ProgressSettings {
    Position,
    Percent,
    PositionSpeed,
    PercentSpeed,
}

fn progress_template(settings: ProgressSettings, grey_msg: bool) -> Result<String> {
    let mut template = String::from_str("{msg:<40} [{bar:40.cyan/blue}]")?;
    if grey_msg {
        template.insert_str(9, "\x1b[0m");
        template.insert_str(0, "\x1b[37m");
    }
    if settings == ProgressSettings::Position || settings == ProgressSettings::PositionSpeed {
        template.push_str(" {human_pos:>7}/{human_len}");
    } else {
        template.push_str(" {percent}%")
    }
    template.push_str(" ETA: {eta}");
    if settings == ProgressSettings::PositionSpeed || settings == ProgressSettings::PercentSpeed {
        template.push_str(" {rate}");
    }
    Ok(template)
}

pub fn progress_bar_with_message(
    length: usize,
    message: &str,
    settings: ProgressSettings,
) -> Result<ProgressBar> {
    let progress = ProgressBar::new(length as u64);
    let template = progress_template(settings, false)?;
    set_progress_style(&progress, &template)?;
    progress.set_message(message.to_owned());
    Ok(progress)
}

fn set_progress_style(progress: &ProgressBar, template: &String) -> Result<()> {
    progress.set_style(
        ProgressStyle::with_template(template)?
            .progress_chars("=> ")
            .with_key(
                "rate",
                |state: &ProgressState, w: &mut dyn std::fmt::Write| {
                    let _ = write!(w, "{}/s", HumanCount(state.per_sec().round() as u64));
                },
            ),
    );
    Ok(())
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
        let (outer, inner) = dual_progress_bars(
            1000,
            Some(1000),
            ProgressSettings::Position,
            ProgressSettings::Position,
            progress,
        )
        .unwrap();
        assert_eq!(outer.length(), Some(1000));
        assert_eq!(inner.length(), Some(1000));
    }
}
