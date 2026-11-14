use anyhow::Context;
use clap::builder::RangedI64ValueParser;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("cannot subtract more cores than available ({})", .0)]
    TooManyCoresSubtractedError(usize),
}

pub fn process_cores_with_available(
    available_cores: usize,
    requested_cores: Option<i32>,
) -> Result<usize, Error> {
    match requested_cores {
        None | Some(0) => Ok(available_cores),
        Some(requested) if requested > 0 => Ok(requested as usize),
        Some(requested) => subtract_from_available_cores(available_cores, requested),
    }
}

pub fn process_cores(requested_cores: Option<i32>) -> Result<usize, Error> {
    let cores = process_cores_with_available(num_cpu(), requested_cores)?;
    if cores == 0 {
        Err(Error::TooManyCoresSubtractedError(num_cpu()))
    } else {
        Ok(cores)
    }
}

pub fn subtract_from_available_cores(
    available_cores: usize,
    requested: i32,
) -> Result<usize, Error> {
    let cores_to_subtract = requested.unsigned_abs() as usize;

    if cores_to_subtract > available_cores {
        Err(Error::TooManyCoresSubtractedError(available_cores))
    } else {
        Ok(available_cores - cores_to_subtract)
    }
}

#[cfg(test)]
fn num_cpu() -> usize {
    4
}

#[cfg(not(test))]
fn num_cpu() -> usize {
    num_cpus::get()
}

pub fn thread_parser() -> RangedI64ValueParser<i32> {
    clap::value_parser!(i32).range(1 + -(num_cpu() as i64)..=num_cpu() as i64)
}
pub fn thread_default() -> i32 {
    num_cpu().try_into().unwrap()
}

pub fn build_pool(cores: usize) -> anyhow::Result<rayon::ThreadPool> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(cores)
        .build()
        .context("Failed to build thread pool")
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::builder::TypedValueParser;
    use pretty_assertions::assert_eq;
    use rstest::{fixture, rstest};
    use std::ffi::OsStr;

    #[fixture]
    fn available_cores() -> usize {
        4
    }

    #[rstest]
    #[case::return_avail_cores(None, 4)]
    #[case::return_avail_cores_0(Some(0), 4)]
    #[case::return_requested_cores(Some(2), 2)]
    #[case::return_requested_cores_1(Some(1), 1)]
    fn test_process_cores(
        available_cores: usize,
        #[case] requested: Option<i32>,
        #[case] expected: usize,
    ) {
        let cores = process_cores_with_available(available_cores, requested)
            .expect("process_cores_with_available should succeed for this test case");
        assert_eq!(cores, expected);
    }

    #[rstest]
    #[case::subtract_one(4, -1, Ok(3))] // Subtract 1 from 4 cores = 3
    #[case::subtract_two(4, -2, Ok(2))] // Subtract 2 from 4 cores = 2
    #[case::subtract_all(4, -4, Ok(0))] // Subtract all cores = 0
    #[case::too_many_cores(4, -5, Err("cannot subtract more cores than available")
    )] // Should error: subtracting more than available
    fn test_subtract_from_available_cores(
        #[case] available: usize,
        #[case] requested: i32,
        #[case] expected: Result<usize, &str>,
    ) {
        let result = subtract_from_available_cores(available, requested);

        match expected {
            Ok(exp) => assert_eq!(result.unwrap(), exp),
            Err(_) => assert!(result.is_err()),
        }
    }

    #[test]
    fn test_subtract_too_many_cores_error_message() {
        let result = subtract_from_available_cores(4, -10);
        assert!(result.is_err());
    }

    #[rstest]
    fn test_thread_default() {
        assert_eq!(thread_default(), 4);
    }

    #[rstest]
    #[case("1", 1)]
    #[case("4", 4)]
    #[case("-3", -3)]
    fn test_thread_parser_accepts_valid_values(#[case] raw: &str, #[case] expected: i32) {
        let parser = thread_parser();
        let value = parser
            .parse_ref(&clap::Command::new("test"), None, OsStr::new(raw))
            .expect("value should parse");
        assert_eq!(value, expected);
    }

    #[rstest]
    #[case("-4")]
    #[case("5")]
    fn test_thread_parser_rejects_invalid_values(#[case] raw: &str) {
        let parser = thread_parser();
        assert!(
            parser
                .parse_ref(&clap::Command::new("test"), None, OsStr::new(raw))
                .is_err()
        );
    }
}
