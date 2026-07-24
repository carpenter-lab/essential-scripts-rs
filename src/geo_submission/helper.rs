use crate::geo_submission::Progress;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use num_cpus;
use std::sync::Arc;

pub fn make_progress_bar(
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

pub const TOO_MANY_CORES_SUBTRACTED_ERROR: &str = "cannot subtract more cores than available";

pub fn process_cores_with_available(
    available_cores: usize,
    requested_cores: Option<i32>,
) -> Result<usize, String> {
    match requested_cores {
        None | Some(0) => Ok(available_cores),
        Some(requested) if requested > 0 => Ok(requested as usize),
        Some(requested) => subtract_from_available_cores(available_cores, requested),
    }
}

pub fn process_cores(requested_cores: Option<i32>) -> Result<usize, String> {
    let cores = process_cores_with_available(num_cpus::get(), requested_cores)?;
    if cores == 0 {
        Err(TOO_MANY_CORES_SUBTRACTED_ERROR.to_string())
    } else {
        Ok(cores)
    }
}

pub fn subtract_from_available_cores(
    available_cores: usize,
    requested: i32,
) -> Result<usize, String> {
    let cores_to_subtract = requested.unsigned_abs() as usize;

    if cores_to_subtract > available_cores {
        Err(TOO_MANY_CORES_SUBTRACTED_ERROR.to_string())
    } else {
        Ok(available_cores - cores_to_subtract)
    }
}
