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
