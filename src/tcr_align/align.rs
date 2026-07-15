use parasail_rs::aligner::Aligner;
use parasail_rs::matrix::Matrix;
use polars::error::{PolarsError, PolarsResult};
use std::sync::OnceLock;

static BLOSUM62: OnceLock<Result<Matrix, String>> = OnceLock::new();

/// Retrieves the BLOSUM62 substitution matrix, used in bioinformatics for scoring alignments
/// between sequences such as proteins. The matrix is cached for efficiency to prevent
/// reloading it multiple times.
///
/// # Returns
/// - `Ok(&'static Matrix)` containing a reference to the BLOSUM62 matrix if successfully loaded.
/// - `Err(PolarsError::ComputeError)` if there was an error loading the matrix.
///
/// The cache is initialized using `BLOSUM62.get_or_init` which attempts to load the matrix
/// from file or another source. On failure, the error is wrapped into a `PolarsError::ComputeError`.
///
/// # Errors
/// Returns a `PolarsError::ComputeError` if the matrix could not be loaded, along with the
/// underlying error message.
///
/// # Example
/// ```ignore
/// let matrix = blosum62();
/// // Use the matrix for sequence alignment or scoring
/// ```
fn blosum62() -> PolarsResult<&'static Matrix> {
    let cached = BLOSUM62.get_or_init(|| {
        Matrix::from("blosum62").map_err(|e| format!("failed to load blosum62 matrix: {e}"))
    });

    match cached {
        Ok(m) => Ok(m),
        Err(msg) => Err(PolarsError::ComputeError(msg.clone().into())),
    }
}

/// Constructs a protein sequence aligner with specified gap penalty parameters.
///
/// This function initializes and builds a global protein sequence aligner using the BLOSUM62
/// substitution matrix. The gap penalty parameters `gap_open` and `gap_extend` control the
/// penalties for opening and extending gaps in the alignment, respectively.
///
/// # Arguments
///
/// - `gap_open` - An `i32` representing the penalty for opening a gap in the alignment.
/// - `gap_extend` - An `i32` representing the penalty for extending an existing gap in the alignment.
///
/// # Returns
///
/// - `PolarsResult<Aligner>` - A result containing the constructed `Aligner` instance if successful,
///   or an error if the matrix initialization fails.
///
/// # Errors
///
/// This function may return an error in the event that the BLOSUM62 substitution matrix cannot
/// be retrieved or cloned.
///
/// # Example
///
/// ```ignore
/// let gap_open = -10;
/// let gap_extend = -1;
/// let aligner = build_protein_aligner(gap_open, gap_extend).expect("Failed to build aligner");
/// ```
///
/// # Notes
///
/// - This implementation ensures that the matrix is cloned into the aligner to avoid any
///   lifetime issues, as the substitution matrix is typically a lightweight handle.
///
/// - The aligner is configured for global sequence alignment by default.
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

/// Executes pairwise sequence alignment using the given `Aligner` over a set of input sequences
/// and optionally a reference set of sequences. Computes the average alignment score.
///
/// # Parameters
/// - `aligner`: A reference to an `Aligner` object that performs the sequence alignment.
/// - `a`: A slice of query input sequences, where each sequence implements the `AsRef<str>` trait.
/// - `reference`: An optional slice of reference sequences. If `None`, the function performs
///   pairwise alignment between all sequences in `a`.
///
/// # Returns
/// - `PolarsResult<f64>`:
///   - On success, returns the average alignment score as a `f64`.
///   - Returns `NaN` if there are no valid alignments (`n == 0`).
///   - On failure, returns a `PolarsError` if an error occurs during the alignment process or other
///     computation-related issues.
///
/// # Errors
/// - Returns `PolarsError::ComputeError` in the following cases:
///   - If the input sequence list `a` is empty.
///   - If the alignment operation in the `Aligner` fails, encapsulating the specific error message.
///
/// # Alignment Procedure
/// - If `reference` is `Some(b)`:
///   - Each sequence in `a` is aligned against every sequence in `b`.
/// - If `reference` is `None`:
///   - Pairwise alignment is performed within the list of sequences in `a`. Each sequence is aligned
///     with every other sequence in `a` once (excluding self-alignments).
///
/// # Example
/// ```ignore
/// use some_alignment_lib::{Aligner, run_parasail};
///
/// let aligner = Aligner::new(); // Aligner initialization
/// let sequences = vec!["ATCG", "GCTA", "TACG"];
/// let reference_sequences = vec!["GCTAGC", "ATCGAT"];
///
/// // Align against reference sequences
/// let score_with_ref = run_parasail(&aligner, &sequences, Some(&reference_sequences)).unwrap();
///
/// // Perform pairwise alignment within `sequences` only
/// let score_without_ref = run_parasail(&aligner, &sequences, None).unwrap();
/// ```
///
/// # Notes
/// - This function leverages parallel alignment computations via the `Aligner`.
/// - The alignment scores are summed, and the average is computed as `(sum of scores) / (number of alignments)`.
///   If no alignments are performed, the function safely returns `NaN`.
///
/// # Panics
/// This function does not explicitly panic but may propagate panics from the `Aligner`'s `align` method
/// if it does not properly handle internal errors.
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

    if let Some(b) = reference {
        for query_seq in a {
            let q_bytes = query_seq.as_ref().as_bytes();
            for ref_seq in b {
                let r_bytes = ref_seq.as_ref().as_bytes();
                let result = aligner.align(Some(q_bytes), r_bytes).map_err(|e| {
                    PolarsError::ComputeError(format!("parasail align failed: {e}").into())
                })?;
                sum += i64::from(result.get_score());
                n += 1;
            }
        }
    } else {
        let n_a = a.len();
        for i in 0..n_a {
            let q_bytes = a[i].as_ref().as_bytes();
            for item in a.iter().take(n_a).skip(i + 1) {
                let r_bytes = item.as_ref().as_bytes();
                let result = aligner.align(Some(q_bytes), r_bytes).map_err(|e| {
                    PolarsError::ComputeError(format!("parasail align failed: {e}").into())
                })?;
                sum += i64::from(result.get_score());
                n += 1;
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
