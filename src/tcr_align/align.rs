use parasail_rs::aligner::Aligner;
use parasail_rs::matrix::Matrix;
use polars::error::{PolarsError, PolarsResult};
use std::sync::OnceLock;

static BLOSUM62: OnceLock<Result<Matrix, String>> = OnceLock::new();

fn blosum62() -> PolarsResult<&'static Matrix> {
    let cached = BLOSUM62.get_or_init(|| {
        Matrix::from("blosum62").map_err(|e| format!("failed to load blosum62 matrix: {e}"))
    });

    match cached {
        Ok(m) => Ok(m),
        Err(msg) => Err(PolarsError::ComputeError(msg.clone().into())),
    }
}

pub fn build_protein_aligner(gap_open: i32, gap_extend: i32) -> PolarsResult<Aligner> {
    // Clone the cached matrix handle into this aligner.
    // (Matrix is typically a thin handle; cloning avoids lifetime issues.)
    let matrix = blosum62()?.clone();

    Ok(Aligner::new()
        .global()
        .matrix(matrix)
        .gap_open(gap_open)
        .gap_extend(gap_extend)
        .build())
}

pub fn run_parasail<S: AsRef<str>>(
    aligner: &Aligner,
    a: &[S],
    reference: Option<&[S]>,
) -> PolarsResult<f64> {
    if a.is_empty() {
        return Err(PolarsError::ComputeError(
            "Cannot run alignment: input sequence list is empty".into(),
        ));
    }

    let mut sum: i64 = 0;
    let mut n: i64 = 0;

    match reference {
        Some(b) => {
            for query_seq in a.iter() {
                let q_bytes = query_seq.as_ref().as_bytes();
                for ref_seq in b.iter() {
                    let r_bytes = ref_seq.as_ref().as_bytes();
                    let result = aligner.align(Some(q_bytes), r_bytes).map_err(|e| {
                        PolarsError::ComputeError(format!("parasail align failed: {e}").into())
                    })?;
                    sum += result.get_score() as i64;
                    n += 1;
                }
            }
        }
        None => {
            let n_a = a.len();
            for i in 0..n_a {
                let q_bytes = a[i].as_ref().as_bytes();
                for j in (i + 1)..n_a {
                    let r_bytes = a[j].as_ref().as_bytes();
                    let result = aligner.align(Some(q_bytes), r_bytes).map_err(|e| {
                        PolarsError::ComputeError(format!("parasail align failed: {e}").into())
                    })?;
                    sum += result.get_score() as i64;
                    n += 1;
                }
            }
        }
    }

    Ok(if n == 0 {
        f64::NAN
    } else {
        (sum as f64) / (n as f64)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use polars::error::PolarsError;
    use rstest::{fixture, rstest};

    const GAP_OPEN: i32 = 10;
    const GAP_EXTEND: i32 = 1;

    #[fixture]
    fn protien_align_fixture() -> Aligner {
        build_protein_aligner(GAP_OPEN, GAP_EXTEND).expect("aligner build")
    }

    #[rstest]
    fn empty_input_returns_err(#[from(protien_align_fixture)] aligner: Aligner) {
        let res = run_parasail(&aligner, &[] as &[&str], None);
        assert!(res.is_err(), "expected error for empty input");
        match res {
            Err(PolarsError::ComputeError(err)) => {
                let msg = format!("{}", err);
                assert!(msg.contains("Cannot run alignment"));
            }
            Err(e) => panic!("unexpected error variant: {:?}", e),
            Ok(_) => panic!("expected Err, got Ok"),
        }
    }

    #[rstest]
    fn empty_reference_produces_nan(#[from(protien_align_fixture)] aligner: Aligner) {
        let a = ["A", "C"];
        let res = run_parasail(&aligner, &a, Some(&[] as &[&str])).expect("call ok");
        assert!(res.is_nan());
    }

    #[rstest]
    fn single_sequence_no_reference_produces_nan(#[from(protien_align_fixture)] aligner: Aligner) {
        let a = ["ACDEF"];
        let res = run_parasail(&aligner, &a, None).expect("call ok");
        assert!(res.is_nan());
    }

    #[rstest]
    fn reference_pairing_mean_matches_manual_sum(#[from(protien_align_fixture)] aligner: Aligner) {
        let a = ["ACD", "WQ"];
        let b = ["ACD", "WQ", "MN"];

        let mean = run_parasail(&aligner, &a, Some(&b)).expect("call ok");

        // Manual calculation
        let mut sum: i64 = 0;
        let mut n: i64 = 0;
        for q in a.iter() {
            for r in b.iter() {
                let result = aligner
                    .align(Some(q.as_bytes()), r.as_bytes())
                    .expect("align ok");
                sum += result.get_score() as i64;
                n += 1;
            }
        }
        let expected = if n == 0 {
            f64::NAN
        } else {
            (sum as f64) / (n as f64)
        };
        if expected.is_nan() {
            assert!(mean.is_nan());
        } else {
            assert!(
                (mean - expected).abs() < 1e-9,
                "mean {} != expected {}",
                mean,
                expected
            );
        }
    }

    #[rstest]
    fn no_reference_pairing_mean_matches_manual_sum(
        #[from(protien_align_fixture)] aligner: Aligner,
    ) {
        let a = ["ACD", "WQ", "MN"];

        let mean = run_parasail(&aligner, &a, None).expect("call ok");

        let mut sum: i64 = 0;
        let mut n: i64 = 0;
        for i in 0..a.len() {
            for j in (i + 1)..a.len() {
                let result = aligner
                    .align(Some(a[i].as_bytes()), a[j].as_bytes())
                    .expect("align ok");
                sum += result.get_score() as i64;
                n += 1;
            }
        }
        let expected = if n == 0 {
            f64::NAN
        } else {
            (sum as f64) / (n as f64)
        };
        if expected.is_nan() {
            assert!(mean.is_nan());
        } else {
            assert!(
                (mean - expected).abs() < 1e-9,
                "mean {} != expected {}",
                mean,
                expected
            );
        }
    }

    // Test helper to set the BLOSUM62 OnceLock to a failure. This is only available in tests.
    #[cfg(test)]
    pub(crate) fn set_blosum62_to_err_for_test(msg: &str) -> Result<(), String> {
        BLOSUM62
            .set(Err(msg.to_string()))
            .map_err(|_| "BLOSUM62 already initialized".to_string())
    }

    #[test]
    #[ignore = "for manual testing only"]
    fn test_blosum62_error_handling() {
        // Force the BLOSUM62 to an error state
        set_blosum62_to_err_for_test("test error").expect("set error");

        // Now, building the aligner should return an error
        let res = build_protein_aligner(GAP_OPEN, GAP_EXTEND);
        assert!(res.is_err(), "expected error when building aligner");

        // Optionally, you can check the specific error message
        if let Err(e) = res {
            match e {
                PolarsError::ComputeError(msg) => {
                    assert!(msg.contains("test error"), "unexpected error message");
                }
                _ => panic!("unexpected error variant: {:?}", e),
            }
        }
    }
}
