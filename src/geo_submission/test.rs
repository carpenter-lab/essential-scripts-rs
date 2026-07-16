use super::*;
use rstest::{fixture, rstest};

#[rstest]
#[case::progress_bar(Progress::Progress, 1000, 1000)]
#[case::progress_bar_hidden(Progress::NoProgress, 1000, 1000)]
#[case::progress_bar_zero_bytes(Progress::Progress, 0, 0)]
fn test_make_progress_bar(
    #[case] progress: Progress,
    #[case] total_bytes: u64,
    #[case] expected_length: u64,
) {
    let result = make_progress_bar(total_bytes, progress);
    assert!(result.is_ok());

    let pb = result.unwrap();
    assert_eq!(pb.length(), Some(expected_length));
}

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
#[case::too_many_cores(4, -5, Err(TOO_MANY_CORES_SUBTRACTED_ERROR)
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
    assert_eq!(result.unwrap_err(), TOO_MANY_CORES_SUBTRACTED_ERROR);
}
