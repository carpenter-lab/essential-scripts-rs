use crate::io::WriteToCsvOrStdout;
use crate::tcr_align::dataframe;
use polars::error::PolarsResult;
use polars::frame::DataFrame;
use polars::prelude::{
    DataFrameJoinOps, JoinArgs, JoinType, LazyCsvReader, LazyFileListReader, PlRefPath,
};
use std::path::PathBuf;

/// A trait that represents a runnable operation or computation.
///
/// Types implementing this trait are expected to define a custom behavior
/// for running a task or computation a specified number of times
/// (`n_replicates`) and returning the result.
///
/// # Associated Types
/// - `Output`: The type of the result produced by the `run` method.
///
/// # Required Methods
/// - `fn run(&self, n_replicates: usize) -> PolarsResult<Self::Output>`
///
///   Executes the operation or computation `n_replicates` times and
///   returns the result wrapped in a `PolarsResult`. The exact behavior
///   of this function depends on the implementation.
///
/// # Parameters
/// - `n_replicates`: Number of times the operation or computation should
///   be executed.
///
/// # Returns
/// - A `PolarsResult<Self::Output>` containing the result of the computation,
///   or an error value if the computation fails.
///
/// # Examples
/// ```ignore
/// use your_crate::Run; // Replace `your_crate` with the actual crate name.
///
/// struct MyRunner;
///
/// impl Run for MyRunner {
///     type Output = i32;
///
///     fn run(&self, n_replicates: usize) -> PolarsResult<Self::Output> {
///         // Example implementation here
///         Ok((n_replicates * 2) as i32)  // Example logic
///     }
/// }
///
/// let runner = MyRunner;
/// let result = runner.run(5).unwrap();
/// assert_eq!(result, 10);
/// ```
///
/// # Errors
/// - The method can return an error wrapped in `PolarsResult` if any issue
///   occurs during the execution of the computation.
pub trait Run {
    type Output;
    fn run(&self, n_replicates: usize) -> PolarsResult<Self::Output>;
}

impl Run for DataFrame {
    type Output = DataFrame;

    /// Executes the primary logic for the `DataFrame` instance, ensuring that specific data manipulations
    /// and computations are performed. This method processes groups prepared from the input `DataFrame`,
    /// checks uniqueness of certain alpha values, and calculates results based on fractional comparisons.
    ///
    /// # Parameters
    /// - `self`: A reference to the `DataFrame` instance on which the method is invoked.
    /// - `n_replicates`: The number of replicates used in the fractional comparison process.
    ///
    /// # Returns
    /// - `PolarsResult<Self::Output>`: A `Result` type that either contains the computed output of the method
    ///   or an error if the processing fails.
    ///
    /// # Workflow
    /// 1. Calls `dataframe::all_unique_cdr3_alpha` to ensure all CDR3 alpha chains are unique within the
    ///    `DataFrame`.
    /// 2. Prepares parasail groups using `dataframe::prepare_parasail_groups`. If an error occurs during group
    ///    preparation, it logs the error and returns it.
    /// 3. If parasail groups are successfully prepared, invokes `dataframe::fraction_self_greater` to compute
    ///    results based on the prepared data, the unique alpha chains, the given number of replicates, and
    ///    additional hardcoded parameters.
    ///
    /// # Errors
    /// - Returns an error if:
    ///   - The uniqueness check for CDR3 alpha chains fails.
    ///   - The preparation of parasail groups fails during processing.
    ///   - Any downstream operations invoked within `fraction_self_greater` result in an error.
    ///
    /// # Example
    /// ```ignore
    /// use my_crate::dataframe::DataFrame;
    ///
    /// let df = DataFrame::new(); // Example initialization
    /// let n_replicates = 100;
    /// match df.run(n_replicates) {
    ///     Ok(output) => println!("Operation completed successfully: {:?}", output),
    ///     Err(err) => eprintln!("Error during processing: {:?}", err),
    /// }
    /// ```
    fn run(self: &DataFrame, n_replicates: usize) -> PolarsResult<Self::Output> {
        let all_unique_alpha = dataframe::all_unique_cdr3_alpha(self)?;

        //let groups = dataframe::prepare_parasail_groups(self);
        match dataframe::prepare_parasail_groups(self) {
            Err(e) => {
                eprintln!("Error preparing parasail groups: {e}");
                Err(e)
            }
            Ok(groups) => {
                dataframe::fraction_self_greater(&groups, &all_unique_alpha, n_replicates, 7, 1)
            }
        }
    }
}

/// Scores a given `DataFrame` by running a specified number of replicates and
/// performs a right join between the original `DataFrame` and the resulting `DataFrame`.
///
/// # Arguments
///
/// * `df` - A reference to the input `DataFrame` to be scored.
/// * `replicates` - The number of replicates to be run for scoring the `DataFrame`.
///
/// # Returns
///
/// Returns a `PolarsResult` containing the new `DataFrame` after scoring and joining
/// the results. If an error occurs during scoring or while performing the join,
/// it returns the error wrapped in a `PolarsResult`.
///
/// # Errors
///
/// This function may return a `PolarsResult::Err` in the following cases:
/// - If the `df.run` method fails to execute with the given replicates.
/// - If the join operation fails due to incompatible `DataFrame` structures or other issues.
///
/// # Example
///
/// ```ignore
/// use polars::prelude::*;
///
/// let df: DataFrame = /* create or load a DataFrame */;
/// let replicates = 10;
///
/// match score_df(&df, replicates) {
///     Ok(result_df) => {
///         println!("Scored DataFrame: {:?}", result_df);
///     }
///     Err(e) => {
///         eprintln!("Error scoring DataFrame: {:?}", e);
///     }
/// }
/// ```
///
/// # Notes
///
/// - The join operation is performed on the `pattern` column from both the original
///   and resulting `DataFrame` using a `Right Join`.
/// - The `JoinArgs::new(JoinType::Right)` specifies the type of join to be used.
fn score_df(df: &DataFrame, replicates: usize) -> PolarsResult<DataFrame> {
    let res = df.run(replicates)?;
    df.join(
        &res,
        ["pattern"],
        ["pattern"],
        JoinArgs::new(JoinType::Right),
        None,
    )
}

/// Compute the TCR score for alignments from an input CSV file and save the results to an output file.
///
/// This function reads a CSV file containing T-cell receptor (TCR) data, computes scoring information
/// using the provided number of replicates, and writes the output to a specified file.
///
/// # Arguments
///
/// * `input_file` - A `PathBuf` representing the path to the input CSV file containing TCR alignments.
/// * `output_file` - A `PathBuf` representing the destination for the computed TCR score output.
///   If `output_file` points to `stdout`, the output will be printed to the standard output.
/// * `replicates` - The number of replicates to use in performing the scoring computation.
///
/// # Panics
///
/// This function will panic for the following reasons:
/// - If the input file cannot be read or parsed as a `DataFrame` due to invalid format or missing content.
/// - If the scoring operation fails (e.g., due to issues with the provided replicates).
///
/// # Returns
///
/// This function does not return a value. Results are written to the specified `output_file`.
///
/// # Example
///
/// ```ignore
/// use std::path::PathBuf;
///
/// let input = PathBuf::from("tcr_input.csv");
/// let output = PathBuf::from("tcr_output.csv");
/// let replicates = 3;
///
/// tcr_score(input, output, replicates);
/// ```
///
/// Ensure the input CSV file follows the expected format for TCR alignments, or the computation will fail.
pub(crate) fn tcr_score(input_file: PathBuf, output_file: PathBuf, replicates: usize) {
    let df = LazyCsvReader::new(PlRefPath::try_from_pathbuf(input_file).unwrap())
        .finish()
        .unwrap()
        .collect()
        .unwrap();

    score_df(&df, replicates)
        .expect("Failed to score TCR alignments. Please check the input file format and try again.")
        .write_to_csv_or_stdout(output_file);
}

#[cfg(test)]
mod tests {
    use super::*;
    use polars::prelude::*;
    use std::collections::HashSet;

    #[test]
    fn score_df_returns_expected_columns() {
        let df = df![
            "pattern" => ["p1", "p2"],
            "TcRa" => [Some("A;B"), Some("C")],
            "TcRb" => [Some("X;Y"), Some("Z")],
        ]
        .unwrap();

        let scored = score_df(&df, 1).unwrap();

        let cols: HashSet<_> = scored
            .get_column_names()
            .into_iter()
            .map(|s| s.as_str().to_string())
            .collect();

        assert!(cols.contains("pattern"));
        assert!(cols.contains("TcRa"));
        assert!(cols.contains("TcRb"));
        assert!(cols.contains("TcRa_alignment_score"));
        assert!(cols.contains("TcRb_alignment_score"));
        assert!(cols.contains("TcRa_alignment_score_v_background"));
    }
}
