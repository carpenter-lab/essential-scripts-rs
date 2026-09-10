use num_cpus;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("cannot subtract more cores than available ({})", .0)]
    TooManyCoresSubtractedError(usize),
}

pub fn process_cores_with_available(
    available_cores: usize,
    requested_cores: Option<i32>,
) -> std::result::Result<usize, Error> {
    match requested_cores {
        None | Some(0) => Ok(available_cores),
        Some(requested) if requested > 0 => Ok(requested as usize),
        Some(requested) => subtract_from_available_cores(available_cores, requested),
    }
}

pub fn process_cores(requested_cores: Option<i32>) -> std::result::Result<usize, Error> {
    let cores = process_cores_with_available(num_cpus::get(), requested_cores)?;
    if cores == 0 {
        Err(Error::TooManyCoresSubtractedError(num_cpus::get()))
    } else {
        Ok(cores)
    }
}

pub fn subtract_from_available_cores(
    available_cores: usize,
    requested: i32,
) -> std::result::Result<usize, Error> {
    let cores_to_subtract = requested.unsigned_abs() as usize;

    if cores_to_subtract > available_cores {
        Err(Error::TooManyCoresSubtractedError(available_cores))
    } else {
        Ok(available_cores - cores_to_subtract)
    }
}
