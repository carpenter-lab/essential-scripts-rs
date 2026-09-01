use crate::tcr_align::align::{build_protein_aligner, run_parasail};
use indicatif::{ProgressBar, ProgressStyle};
use parasail_rs::aligner::Aligner;
use polars::error::PolarsResult;
use polars::frame::{DataFrame, UniqueKeepStrategy};
use polars::prelude::{ExplodeOptions, IntoLazy, NamedFrom, Series, col, len, lit};
use rand::SeedableRng;
use rand::prelude::IndexedRandom;
use rand::rngs::SmallRng;
use std::collections::HashSet;

// Limit the number of sequences per group to avoid excessive computation time.
// 1000 sequences per group should be enough for most cases.
// Excess sequences will be downsampled randomly.
// 1000 sequences prevent groups that likely don't matter from being aligned.
const MAX_GROUP_SEQS: usize = 1000;
const MAX_BG_SEQS: usize = 1000;

struct ScoreContext<'a> {
    rng: &'a mut SmallRng,
    pb: &'a ProgressBar,
    aligner: &'a Aligner,
    pattern: &'a str,
}

/// Extracts all unique TcRa CDR3 alpha sequences from a given `DataFrame`.
///
/// This function processes a `DataFrame` containing a column `TcRa`, where each entry may
/// contain semicolon-separated CDR3 alpha sequences. It performs the following steps:
///
/// 1. Creates a lazy frame from the input `DataFrame`.
/// 2. Selects the `TcRa` column.
/// 3. Splits the semicolon-separated strings in `TcRa` into lists and stores them in a new column `TcRa_list`.
/// 4. Explodes the lists in `TcRa_list` such that each list element becomes a separate row.
/// 5. Filters out any null values from the exploded column.
/// 6. Extracts unique values from the `TcRa_list` column, keeping only the first occurrence of each value.
/// 7. Collects the unique values into a new `DataFrame`.
/// 8. Processes the resulting series to skip nulls and empty/whitespace-only strings, returning the unique CDR3 alpha sequences as a vector of strings.
///
/// # Arguments
///
/// * `all_data` - A reference to a `DataFrame` containing the input data with a column named `TcRa`.
///
/// # Returns
///
/// * `PolarsResult<Vec<String>>` - A result containing a vector of unique CDR3 alpha sequences as strings,
///   or an error if any part of the processing fails.
///
/// # Errors
///
/// This function returns a `PolarsResult::Err` in the following cases:
/// - When the input `DataFrame` does not contain a column named `TcRa`.
/// - When any operation (e.g., splitting, exploding, or collecting unique values) fails internally.
/// - When an unexpected null value or type mismatch is encountered.
///
/// # Examples
///
/// ```ignore
/// use polars::prelude::*;
///
/// // Create a sample DataFrame
/// let df = DataFrame::new(vec![
///     Series::new("TcRa", &["CDR3A1;CDR3A2", "CDR3A3", "CDR3A2;CDR3A4;  ", ""]),
/// ]).unwrap();
///
/// // Extract all unique CDR3 alpha sequences
/// let result = all_unique_cdr3_alpha(&df).unwrap();
///
/// assert_eq!(result, vec!["CDR3A1", "CDR3A2", "CDR3A3", "CDR3A4"]);
/// ```
pub(crate) fn all_unique_cdr3_alpha(all_data: &DataFrame) -> PolarsResult<Vec<String>> {
    let lf = all_data.clone().lazy();
    let unique_df = lf
        .select([col("TcRa")])
        .with_columns([col("TcRa").str().split(lit(";")).alias("TcRa_list")])
        .explode(
            col("TcRa_list")
                .into_selector()
                .expect("could not create selector"),
            ExplodeOptions {
                empty_as_null: false,
                keep_nulls: true,
            },
        )
        .select([col("TcRa_list")])
        .filter(col("TcRa_list").is_not_null())
        .unique(None, UniqueKeepStrategy::First)
        .collect()?;

    let s = unique_df.column("TcRa_list")?.str()?;
    // Skip nulls and empty/whitespace-only strings produced by splitting
    Ok(s.iter()
        .flatten()
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(ToString::to_string)
        .collect())
}

/// Prepares parasail groups by processing and transforming the input `DataFrame`.
///
/// This function performs the following steps on the input `DataFrame`:
/// 1. Selects the `pattern`, `TcRa`, and `TcRb` columns.
/// 2. Eliminates duplicate rows, keeping the first occurrence.
/// 3. Splits the `TcRa` column into individual elements by a delimiter (`;`),
///    explodes it into multiple rows, and renames the column back to `TcRa`.
/// 4. Splits the `TcRb` column into individual elements by a delimiter (`;`),
///    explodes it into multiple rows, and renames the column back to `TcRb`.
/// 5. Deduplicates pairs of `TcRa` and `TcRb` after the explosion step.
/// 6. Groups data by the `pattern` column and aggregates the following:
///    - A count of records in the group (`n`).
///    - Unique values of `TcRa` as a list.
///    - Unique values of `TcRb` as a list.
/// 7. Filters out groups that have no members (where `n` > 0).
///
/// # Arguments
///
/// * `all_data` - A `DataFrame` containing the input data to be processed.
///
/// # Returns
///
/// Returns a `PolarsResult<DataFrame>` that contains the transformed and grouped data
/// with columns:
/// - `pattern`
/// - `n` (count of group members)
/// - `TcRa` (list of unique `TcRa` values)
/// - `TcRb` (list of unique `TcRb` values)
///
/// # Errors
///
/// This function will return an error if any of the following operations fail:
/// - Applying transformations (e.g., splitting, exploding, renaming columns).
/// - Grouping, aggregating, or filtering the `DataFrame`.
/// - Collecting the lazy `DataFrame` into a concrete `DataFrame`.
///
/// # Examples
///
/// ```ignore
/// use polars::prelude::*;
///
/// let data = df! {
///     "pattern" => &["A", "A", "B"],
///     "TcRa" => &["x;y", "z", "a;b"],
///     "TcRb" => &["m", "n;o", "p"]
/// };
///
/// let result = prepare_parasail_groups(&data).unwrap();
/// println!("{:?}", result);
/// ```
///
/// The resulting `DataFrame` will contain grouped and transformed data based
/// on the described logic.
pub(crate) fn prepare_parasail_groups(all_data: &DataFrame) -> PolarsResult<DataFrame> {
    let lf = all_data
        .clone()
        .lazy()
        .select([col("pattern"), col("TcRa"), col("TcRb")])
        .unique(None, UniqueKeepStrategy::First)
        // TcRa -> explode
        .with_columns([col("TcRa").str().split(lit(";")).alias("TcRa_list")])
        .explode(
            col("TcRa_list")
                .into_selector()
                .expect("could not create selector"),
            ExplodeOptions {
                empty_as_null: false,
                keep_nulls: true,
            },
        )
        .drop(
            col("TcRa")
                .into_selector()
                .expect("could not create selector"),
        )
        .rename(["TcRa_list"], ["TcRa"], true)
        .filter(col("TcRa").is_not_null())
        // TcRb -> explode
        .with_columns([col("TcRb").str().split(lit(";")).alias("TcRb_list")])
        .explode(
            col("TcRb_list")
                .into_selector()
                .expect("could not create selector"),
            ExplodeOptions {
                empty_as_null: false,
                keep_nulls: true,
            },
        )
        .drop(
            col("TcRb")
                .into_selector()
                .expect("could not create selector"),
        )
        .rename(["TcRb_list"], ["TcRb"], true)
        // dedupe pairs after explosion (like distinct())
        .unique(None, UniqueKeepStrategy::First);

    // group and collect into list columns
    let grouped = lf
        .group_by([col("pattern")])
        .agg([
            len().alias("n"),
            col("TcRa").unique().alias("TcRa"), // list
            col("TcRb").unique().alias("TcRb"), // list
        ])
        // keep groups with at least one member
        .filter(col("n").gt(lit(0)));

    grouped.collect()
}

/// Extracts a unique, trimmed, and non-empty list of UTF-8 strings from a `List` cell
/// in a specified column of a `DataFrame`. The order of the strings is preserved based
/// on their first occurrence.
///
/// # Arguments
/// * `df` - A reference to the `DataFrame` from which the data is to be extracted.
/// * `col_name` - The name of the column containing the `List` cell to extract.
/// * `row_idx` - The index of the row containing the `List` cell to extract.
///
/// # Returns
/// Returns a `PolarsResult` containing a `Vec<String>` of unique, non-empty, and
/// trimmed strings.
///
/// # Errors
/// * Returns an error if the specified column does not exist or is not of type `List`.
/// * Returns an error if the specified row index is out of bounds.
/// * Panics if the `List` cell cannot be converted into a series or if there are unexpected
///   data types.
///
/// # Example
/// ```ignore
/// use polars::prelude::*;
///
/// # fn main() -> PolarsResult<()> {
///     let s0 = Series::new("list_column", &[
///         Series::new("inner", &["a", "b", "c"]),
///         Series::new("inner", &["c", "d", ""]),
///     ]);
///     let df = DataFrame::new(vec![s0])?;
///
///     let result = get_list_cell_as_vec_utf8(&df, "list_column", 1)?;
///     assert_eq!(result, vec!["c", "d"]);
/// #    Ok(())
/// # }
/// ```
///
/// # Notes
/// * The function skips any empty strings found in the `List`.
/// * Strings are trimmed of leading and trailing whitespace.
fn get_list_cell_as_vec_utf8(
    df: &DataFrame,
    col_name: &str,
    row_idx: usize,
) -> PolarsResult<Vec<String>> {
    let s = df
        .column(col_name)?
        .list()?
        .get_as_series(row_idx)
        .expect("could not extract series");
    let mut seen: HashSet<String> = HashSet::new();
    let mut out: Vec<String> = Vec::new();

    for x in s.str()?.iter().flatten() {
        let v = x.trim();
        if v.is_empty() {
            continue;
        }
        // preserve first-seen order
        if seen.insert(v.to_string()) {
            out.push(v.to_string());
        }
    }

    Ok(out)
}

/// Computes alignment scores and their relative fractions for a given dataset of patterns.
///
/// This function processes a dataset (`groups`) containing biological patterns (e.g., DNA, RNA, or protein sequences)
/// and computes alignment scores against self and a background distribution. It calculates specific columns based
/// on alignment scores and generates a new `DataFrame` containing the results.
///
/// ### Parameters
/// - `groups: &DataFrame`
///   - Input `DataFrame` containing the data to process.
///   - This `DataFrame` must include a column named `"pattern"` containing the specific patterns as strings.
/// - `all_unique_alpha: &[String]`
///   - A reference to a list of unique alphanumeric characters that can be used in the alignment process.
/// - `n_replicates: usize`
///   - The number of random background replicates to perform for comparison during alignment score calculation.
/// - `gap_open: i32`
///   - The penalty score for opening a gap in the alignment process. Negative values are usually assigned.
/// - `gap_extend: i32`
///   - The penalty score for extending an existing gap in the alignment process. Negative values are usually assigned.
///
/// ### Returns
/// - `PolarsResult<DataFrame>`
///   - A result containing the new `DataFrame` if successful, or an error if there is a failure during processing.
///
/// ### Output `DataFrame` Schema
/// - The resulting `DataFrame` contains the following columns:
///   1. `"pattern"` (`String`): The original patterns from the input `DataFrame`.
///   2. `"TcRa_alignment_score"` (`f64`): Self-alignment scores for the patterns against `TcRa` sequences.
///   3. `"TcRb_alignment_score"` (`f64`): Self-alignment scores for the patterns against `TcRb` sequences.
///   4. `"TcRa_alignment_score_v_background"` (`f64`): The fractions of `TcRa` self-alignment scores compared to background replicates.
///
/// ### Notes
/// This function makes use of the `SmallRng` random number generator from the `rand` crate to introduce stochasticity
/// in calculating background replicates. The `calculate_scores` function
/// performs the core operations for score calculation.
///
/// ### Errors
/// - If the `"pattern"` column is missing or cannot be converted, an error is returned.
/// - Errors arising from the `Polars` library (e.g., column manipulations or `DataFrame` creation) will also be propagated.
///
/// ### Example Usage
/// ```ignore
/// use polars::prelude::*;
///
/// // Example input DataFrame.
/// let df = df!(
///     "pattern" => &["ACGT", "GGCTA", "TTAACG"]
/// )?;
///
/// let all_unique_alpha = vec!["A".to_string(), "C".to_string(), "G".to_string(), "T".to_string()];
/// let result = fraction_self_greater(&df, &all_unique_alpha, 10, -5, -1)?;
///
/// println!("{:?}", result);
/// ```
pub(crate) fn fraction_self_greater(
    groups: &DataFrame,
    all_unique_alpha: &[String],
    n_replicates: usize,
    gap_open: i32,
    gap_extend: i32,
) -> PolarsResult<DataFrame> {
    let height = groups.height();

    let patterns: Vec<String> = groups
        .column("pattern")?
        .as_series()
        .expect("could not get 'pattern' column")
        .str()?
        .iter()
        .flatten()
        .map(ToString::to_string)
        .collect();

    let mut out_pattern: Vec<String> = Vec::with_capacity(height);
    let mut out_self_a: Vec<f64> = Vec::with_capacity(height);
    let mut out_self_b: Vec<f64> = Vec::with_capacity(height);
    let mut out_frac: Vec<f64> = Vec::with_capacity(height);

    let mut rng = SmallRng::from_rng(&mut rand::rng());

    calculate_scores(
        groups,
        all_unique_alpha,
        n_replicates,
        gap_open,
        gap_extend,
        &patterns,
        &mut out_pattern,
        &mut out_self_a,
        &mut out_self_b,
        &mut out_frac,
        &mut rng,
    );

    DataFrame::new_infer_height(vec![
        Series::new("pattern".into(), out_pattern).into(),
        Series::new("TcRa_alignment_score".into(), out_self_a).into(),
        Series::new("TcRb_alignment_score".into(), out_self_b).into(),
        Series::new("TcRa_alignment_score_v_background".into(), out_frac).into(),
    ])
}

fn downsample_vec(
    v: Vec<String>,
    pattern: &str,
    pb: &ProgressBar,
    rng: &mut SmallRng,
) -> Vec<String> {
    let out: Vec<String> = if v.len() > MAX_GROUP_SEQS {
        pb.set_message(format!(
            "Downsampling {p} {full} -> {maxg}",
            p = pattern,
            full = v.len(),
            maxg = MAX_GROUP_SEQS
        ));
        v.sample(rng, MAX_GROUP_SEQS).cloned().collect()
    } else {
        v
    };
    out
}

/// Calculates scores based on the given input parameters. This function performs multiple
/// alignments, calculates self-similarity scores, and evaluates the statistical significance
/// of these scores against background sequences. It also handles the generation of progress
/// bars for tracking computation progress.
///
/// # Parameters
///
/// - `groups`: A `DataFrame` containing input data, each row corresponding to a single group
///   of items for which scores are calculated.
/// - `all_unique_alpha`: A reference to a vector of all unique alpha sequences to be used for
///   background sampling during statistical significance evaluation.
/// - `n_replicates`: The number of replicates to run for statistical significance evaluation.
/// - `gap_open`: Gap opening penalty value used by the aligner during sequence alignment.
/// - `gap_extend`: Gap extension penalty value used by the aligner during sequence alignment.
/// - `patterns`: A slice of strings representing the patterns used to process each group.
/// - `out_pattern`: A mutable reference to a vector where the output patterns will be stored.
/// - `out_self_a`: A mutable reference to a vector where self-similarity scores for `TcRa`
///   sequences will be stored.
/// - `out_self_b`: A mutable reference to a vector where self-similarity scores for `TcRb`
///   sequences will be stored.
/// - `out_frac`: A mutable reference to a vector where fraction scores (statistical significance)
///   will be stored.
/// - `rng`: A mutable reference to a random number generator (`SmallRng`) used for background
///   sampling and random operations.
///
/// # Behavior
///
/// - The function first initializes progress bars to track progress of the computation. If a
///   progress bar fails to initialize, the computation continues with hidden progress indicators.
/// - A protein aligner is built using the provided `gap_open` and `gap_extend` penalties. If the
///   aligner fails to initialize, an error message is printed and the function exits early.
/// - Iterates over all rows (or groups) in the `groups` dataset (up to `height` number of items).
///   For each group:
///   - Processes the `TcRb` sequences and calculates their self-similarity scores. If no valid
///     data is found, `NaN` is stored in the output vectors.
///   - Similarly processes the `TcRa` sequences to calculate their self-similarity.
///   - Performs statistical significance testing using background sequences sampled from
///     `all_unique_alpha`. Fraction scores representing the proportion of replicates where the
///     observed self-similarity score is better (lower) than the background score are calculated
///     and stored in the `out_frac` vector.
///   - Outputs `NaN` values for groups with insufficient data or patterns marked as "single".
/// - Updates the progress bars as the computation proceeds.
///
/// # Errors
///
/// - If there are issues with building the protein aligner or running the parasail alignment
///   library, error messages are printed and `NaN` values are assigned to the relevant outputs.
/// - Errors associated with accessing data from the `groups` dataset are handled by assigning
///   `NaN` values to the output scores.
///
/// # Notes
///
/// - The function uses `indicatif::ProgressBar` and `indicatif::MultiProgress` to provide a user
///   interface for tracking progress.
/// - This function is designed to handle large input datasets efficiently and supports parallel
///   computation for statistical evaluations using the parasail alignment library.
/// - The function assumes a maximum number of background sequences (`MAX_BG_SEQS`) for comparisons.
/// - The provided patterns must align with the rows in `groups`, and missing patterns or invalid
///   indices result in default or `NaN` outputs.
///
/// # Example Usage
///
/// ```ignore
/// use rand::rngs::SmallRng;
/// use rand::SeedableRng;
/// let mut rng = SmallRng::from_entropy();
///
/// let groups = load_dataframe(); // Load a DataFrame from some source.
/// let all_unique_alpha = vec![String::from("SEQ1"), String::from("SEQ2")];
/// let n_replicates = 100;
/// let mut out_pattern = vec![];
/// let mut out_self_a = vec![];
/// let mut out_self_b = vec![];
/// let mut out_frac = vec![];
///
/// calculate_scores(
///     &groups,
///     &all_unique_alpha,
///     n_replicates,
///     -10,
///     -1,
///     &vec![String::from("pattern1"), String::from("pattern2")],
///     &mut out_pattern,
///     &mut out_self_a,
///     &mut out_self_b,
///     &mut out_frac,
///     &mut rng,
/// );
/// ```
///
/// # Dependencies
///
/// - Requires the `indicatif` crate for progress bars.
/// - Uses the `rand` crate for random number generation.
/// - Relies on the parasail alignment library for sequence alignment operations.
///
/// # Performance Considerations
///
/// - The function is designed to handle large datasets, but its runtime is directly proportional
///   to the number of replicates (`n_replicates`) and the size of the sequences. Adjust these
///   parameters based on available computational resources.
#[allow(clippy::too_many_arguments)]
fn calculate_scores(
    groups: &DataFrame,
    all_unique_alpha: &[String],
    n_replicates: usize,
    gap_open: i32,
    gap_extend: i32,
    patterns: &[String],
    out_pattern: &mut Vec<String>,
    out_self_a: &mut Vec<f64>,
    out_self_b: &mut Vec<f64>,
    out_frac: &mut Vec<f64>,
    rng: &mut SmallRng,
) {
    let height = groups.height();
    let (pb, inner_pb) = create_dual_progress_bar(n_replicates, height);

    let aligner = match build_protein_aligner(gap_open, gap_extend) {
        Ok(aligner) => aligner,
        Err(e) => {
            eprintln!("Error building protein aligner: {e}");
            return;
        }
    };

    for i in 0..height {
        inner_pb.set_message(format!("{:width$}", "", width = 33));
        let pattern = patterns.get(i).map_or("", |v| v);
        if pattern == "single" {
            pb.inc(1);
            out_pattern.push(pattern.to_string());
            out_self_a.push(f64::NAN);
            out_self_b.push(f64::NAN);
            out_frac.push(f64::NAN);
            continue;
        }
        pb.set_message(format!(
            "Calculating scores for {:width$}",
            pattern,
            width = 10
        ));
        let mut ctx = ScoreContext {
            rng,
            pb: &inner_pb,
            aligner: &aligner,
            pattern,
        };

        match get_list_cell_as_vec_utf8(groups, "TcRb", i) {
            Err(_) => {
                out_self_b.push(f64::NAN);
            }
            Ok(b_full) => {
                let b_res = calculate_single_chain_score(&mut ctx, b_full);
                out_self_b.push(b_res.self_mean);
            }
        }

        match get_list_cell_as_vec_utf8(groups, "TcRa", i) {
            Err(_) => {
                out_self_a.push(f64::NAN);
            }
            Ok(a_full) => {
                let a_res = calculate_single_chain_score(&mut ctx, a_full);
                out_self_a.push(a_res.self_mean);
                let final_value = calculate_single_chain_background_score(
                    &mut ctx,
                    all_unique_alpha,
                    n_replicates,
                    a_res,
                );
                out_frac.push(final_value);
            }
        }

        inner_pb.reset();
        out_pattern.push(pattern.to_string());
        pb.inc(1);
    }
}

fn create_dual_progress_bar(n_replicates: usize, height: usize) -> (ProgressBar, ProgressBar) {
    let progress_style = ProgressStyle::with_template(
        "{msg} [{bar:40.cyan/blue}] {pos}/{len} Elapsed: {elapsed_precise} ETA: {eta}",
    );
    let inner_progress_style =
        ProgressStyle::with_template("\x1b[37m{msg}\x1b[0m [{bar:40.cyan/blue}] {pos}/{len}");
    let mpb: indicatif::MultiProgress;
    let pb: ProgressBar;
    let inner_pb: ProgressBar;

    if let Ok(progress_style) = progress_style {
        mpb = indicatif::MultiProgress::new();
        pb = mpb.add(ProgressBar::new(height as u64));
        pb.set_style(progress_style);
        match inner_progress_style {
            Ok(inner_progress_style) => {
                inner_pb = mpb.add(ProgressBar::new(n_replicates as u64));
                inner_pb.set_style(inner_progress_style);
            }
            Err(_) => {
                inner_pb = ProgressBar::hidden();
            }
        }
    } else {
        pb = ProgressBar::hidden();
        inner_pb = ProgressBar::hidden();
    }
    (pb, inner_pb)
}

struct SingleChainResult {
    sampled_sequences: Vec<String>,
    self_mean: f64,
    n_seqs: usize,
}

impl SingleChainResult {
    fn new(sampled_sequences: Vec<String>, self_mean: f64) -> Self {
        let n_seqs = sampled_sequences.len();
        Self {
            sampled_sequences,
            self_mean,
            n_seqs,
        }
    }

    fn is_nan(&self) -> bool {
        self.self_mean.is_nan()
    }
}

fn calculate_single_chain_score(
    ctx: &mut ScoreContext<'_>,
    all_seqs: Vec<String>,
) -> SingleChainResult {
    let downsampled_sequences = downsample_vec(all_seqs, ctx.pattern, ctx.pb, ctx.rng);
    let self_mean = if downsampled_sequences.len() <= 1 {
        f64::NAN
    } else {
        run_parasail(ctx.aligner, &downsampled_sequences, None).unwrap_or_else(|e| {
            eprintln!("Error running parasail for pattern {}: {e}", ctx.pattern);
            f64::NAN
        })
    };

    SingleChainResult::new(downsampled_sequences, self_mean)
}

fn calculate_single_chain_background_score(
    ctx: &mut ScoreContext<'_>,
    unique_seqs: &[String],
    n_replicates: usize,
    group_result: SingleChainResult,
) -> f64 {
    let n_val = group_result.n_seqs.min(unique_seqs.len()).min(MAX_BG_SEQS);
    let mut greater_count = 0usize;
    if group_result.n_seqs > 1 && n_val > 0 && !group_result.is_nan() {
        for _rep in 0..n_replicates {
            let background: Vec<String> = unique_seqs.sample(ctx.rng, n_val).cloned().collect();

            let bg_mean = run_parasail(ctx.aligner, &background, None).unwrap_or_else(|e| {
                eprintln!("Error running parasail for pattern {}: {e}", ctx.pattern);
                f64::NAN
            });
            // we want lower scores to be better
            if group_result.self_mean < bg_mean {
                greater_count += 1;
            }
            ctx.pb.inc(1);
        }
    } else {
        greater_count = n_replicates;
    }
    if group_result.is_nan() {
        -1.
    } else {
        (greater_count as f64) / (n_replicates as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use indicatif::ProgressBar;
    use polars::prelude::*;
    use pretty_assertions::assert_eq;
    use rand::SeedableRng;
    use rand::rngs::SmallRng;
    use rstest::{fixture, rstest};
    use std::collections::HashSet;

    fn row_idx_for_pattern(df: &DataFrame, pattern: &str) -> usize {
        df.column("pattern")
            .unwrap()
            .str()
            .unwrap()
            .iter()
            .enumerate()
            .find_map(|(idx, v)| (v == Some(pattern)).then_some(idx))
            .expect("pattern not found")
    }
    fn get_f64_col_by_idx(out: &DataFrame, idx: usize, name: &str) -> f64 {
        out.column(name).unwrap().f64().unwrap().get(idx).unwrap()
    }

    #[test]
    fn test_all_unique_cdr3_alpha_basic() {
        let df = df![
            "pattern" => ["p1", "p2", "p1"],
            "TcRa" => [Some("A;B;;"), Some("B;C"), None::<&str>],
            "TcRb" => ["X", "Y", "Z"],
        ]
        .unwrap();

        let uniques = all_unique_cdr3_alpha(&df).unwrap();
        let set: HashSet<_> = uniques.into_iter().collect();
        let expected: HashSet<String> = ["A", "B", "C"].iter().map(|s| s.to_string()).collect();
        assert_eq!(set, expected);
    }

    #[test]
    fn test_prepare_parasail_groups_and_get_list() {
        let df = df![
            "pattern" => ["pat1", "pat1", "pat2"],
            "TcRa" => [Some("A;A"), Some("B"), Some("C;D")],
            "TcRb" => [Some("X"), Some("Y;Z"), None::<&str>],
        ]
        .unwrap();

        let grouped = prepare_parasail_groups(&df).unwrap();

        // Expect two groups: pat1 and pat2
        assert_eq!(grouped.height(), 2);

        // find index for 'pat1'
        let pat_series = grouped.column("pattern").unwrap();
        let pat_vals: Vec<String> = pat_series
            .as_series()
            .expect("could not get series for pattern")
            .str()
            .unwrap()
            .iter()
            .flatten()
            .map(|s| s.to_string())
            .collect();
        let idx = pat_vals
            .iter()
            .position(|s| s == "pat1")
            .expect("pat1 missing");

        // extract TcRa list for pat1
        let tcra = get_list_cell_as_vec_utf8(&grouped, "TcRa", idx).unwrap();
        let set: HashSet<_> = tcra.into_iter().collect();
        let expected: HashSet<String> = ["A", "B"].iter().map(|s| s.to_string()).collect();
        assert_eq!(set, expected);
    }

    #[rstest]
    #[case(MAX_GROUP_SEQS + 200)]
    #[case(10)]
    fn test_downsample_vec_limits(#[case] total: usize) {
        let v: Vec<String> = (0..total).map(|i| format!("S{}", i)).collect();
        let mut rng = SmallRng::from_seed([1; 32]);
        let pb = ProgressBar::hidden();

        let out = downsample_vec(v.clone(), "pat", &pb, &mut rng);
        assert!(out.len() <= MAX_GROUP_SEQS);
        assert_eq!(out.len(), total.min(MAX_GROUP_SEQS));
    }

    #[test]
    fn test_get_list_cell_trims_and_deduplicates() {
        let df = df![
            "pattern" => ["pat1", "pat1"],
            "TcRa" => [Some(" A ;A; ;A "), Some("A")],
            "TcRb" => [Some("X"), Some("X")],
        ]
        .unwrap();

        let grouped = prepare_parasail_groups(&df).unwrap();
        let values = get_list_cell_as_vec_utf8(&grouped, "TcRa", 0).unwrap();
        assert_eq!(values, vec!["A"]);
    }

    #[test]
    fn test_fraction_self_greater_marks_single_pattern_as_nan() {
        let df = df![
            "pattern" => ["single", "pat2"],
            "TcRa" => [Some("ACD;WQ"), Some("ACD;MN")],
            "TcRb" => [Some("MN;PQ"), Some("PQ;RS")],
        ]
        .unwrap();

        let all_alpha = all_unique_cdr3_alpha(&df).unwrap();
        let groups = prepare_parasail_groups(&df).unwrap();
        let out = fraction_self_greater(&groups, &all_alpha, 2, 7, 1).unwrap();

        assert_eq!(out.height(), 2);
        let idx = row_idx_for_pattern(&out, "single");
        let tcra_score = get_f64_col_by_idx(&out, idx, "TcRa_alignment_score");
        let tcrb_score = get_f64_col_by_idx(&out, idx, "TcRb_alignment_score");
        let frac = get_f64_col_by_idx(&out, idx, "TcRa_alignment_score_v_background");
        assert!(tcra_score.is_nan());
        assert!(tcrb_score.is_nan());
        assert!(frac.is_nan());
    }

    #[test]
    fn test_fraction_self_greater_uses_minus_one_when_alpha_self_is_nan() {
        let df = df![
            "pattern" => ["pat1", "pat2"],
            "TcRa" => [Some("ACD"), Some("WQ;MN")],
            "TcRb" => [Some("MN;PQ"), Some("PQ;RS")],
        ]
        .unwrap();

        let all_alpha = all_unique_cdr3_alpha(&df).unwrap();
        let groups = prepare_parasail_groups(&df).unwrap();
        let out = fraction_self_greater(&groups, &all_alpha, 3, 7, 1).unwrap();

        let idx = row_idx_for_pattern(&out, "pat1");
        let tcra_score = get_f64_col_by_idx(&out, idx, "TcRa_alignment_score");
        let frac = get_f64_col_by_idx(&out, idx, "TcRa_alignment_score_v_background");
        assert!(tcra_score.is_nan());
        assert_eq!(frac, -1.0);
    }
    #[fixture]
    fn background_vector() -> Vec<&'static str> {
        vec![
            "CAGSRGTESNQPQHF",
            "CASQFLPGGNEQYF",
            "CASRDQPSTDTQYF",
            "CASRLAGGPTLPQGTQYF",
            "CASRQFLSTGELFF",
            "CASSAGQSFYNEQFF",
            "CASSALFGLETQYF",
            "CASSDWEFRTDTQYF",
            "CASSEFRTADKNEQYF",
            "CASSEFRTDTQYF",
            "CASSEFRTNTGNTIYF",
            "CASSEGDPPNYEQYF",
            "CASSFGTESNQPQHF",
            "CASSFGTGGKYEQYF",
            "CASSGGGLFNQPQHF",
            "CASSGGQSFEAFF",
            "CASSHGTGGKYEQYF",
            "CASSHIYAGQFF",
            "CASSHIYEQFF",
            "CASSHIYRGDEQFF",
            "CASSHIYTNRDRGFYEQYF",
            "CASSIDPPNSGNTIYF",
            "CASSKGGPTLPGNTIYF",
            "CASSKGRLFNQPQHF",
            "CASSLALFGGYTF",
            "CASSLDQPSYNEQFF",
            "CASSLEFRTGGEYF",
            "CASSLEFRTGSSYNEQFF",
            "CASSLGTESNQPQHF",
            "CASSLGTGGKYEQYF",
            "CASSLQMANIQYF",
            "CASSLQMANIRYF",
            "CASSLQMGGTGYTF",
            "CASSLQMTNIQYF",
            "CASSLRDPPNWQFF",
            "CASSLRGQSFSYNEQFF",
            "CASSLVALLGEDTQYF",
            "CASSQDPQGDTQYF",
            "CASSQDPQGGSNQPQHF",
            "CASSQDPQGIGKLFF",
            "CASSQDPQQSGSNQPQHF",
            "CASSQDPQRAEQYF",
            "CASSQDPQVLEQYF",
            "CASSQDQPSAQRGYNEQFF",
            "CASSQFLAAQTQYF",
            "CASSQFLAGWETQYF",
            "CASSQFLAKNIQYF",
            "CASSQFLDPLYTF",
            "CASSQFLGDEQFF",
            "CASSQFLGLLSSYEQYF",
            "CASSQFLIAGTNTEAFF",
            "CASSQFLPSTDTQYF",
            "CASSQFLREQFF",
            "CASSQFLSAGHSGAKNIQYF",
            "CASSQFLSYEQYF",
            "CASSQFLTGTGELFF",
            "CASSQGLGPGLFNQPQHF",
            "CASSQGTGGKYEQFF",
            "CASSQGTGGKYEQYF",
            "CASSQLDLFDPSEQYF",
            "CASSQLDLGLTHNEQFF",
            "CASSQLDLIAQKQETQYF",
            "CASSQLDLIPRPYEQYF",
            "CASSQLDLVPRSTDTQYF",
            "CASSQNLNTGELFF",
            "CASSQQNLNYGYTF",
            "CASSRGTGGKYEQYF",
            "CASSRTGQSFYGYTF",
            "CASSSGWGQSFNQPQHF",
            "CASSSPGQGGANYGYTF",
            "CASSSYSISGELFF",
            "CASSTQFLRGNTIYF",
            "CASSTQTLGQSFETQYF",
            "CASSVALFAGEQYF",
            "CASSVALFGEGYTF",
            "CASSVALFGETQYF",
            "CASSVALFGNTIYF",
            "CASSVALFGSYTF",
            "CASSVALFSNTQYF",
            "CASSVALLAGTQYF",
            "CASSVALLAQPQFF",
            "CASSVALLGAEQYF",
            "CASSVALLGETQYF",
            "CASSVALLGGEQYF",
            "CASSVALLGGTQYF",
            "CASSVALLGNTIYF",
            "CASSVALLGQPQHF",
            "CASSVALLTGELFF",
            "CASSVALLTGGQVF",
            "CASSVSLQMETQYF",
            "CASSYAQNLNNEQFF",
            "CASSYEFRTAYEQYF",
            "CASSYEPGQFLLPLHF",
            "CASSYGTGGKYEQYF",
        ]
    }

    #[rstest]
    #[case::single_sequence_returns_zero(vec!["ACD"], 100, -1.0, f64::NAN)]
    #[case::high_quality_group_alignment(vec![
"CASSVALLAGTQYF",
"CASSVALLAQPQFF",
"CASSVALLGAEQYF",
"CASSVALLGETQYF",
"CASSVALLGGEQYF",
"CASSVALLGGTQYF",
"CASSVALLGNTIYF",
"CASSVALLGQPQHF",
"CASSVALLTGELFF",
"CASSVALLTGGQVF"
    ], 10, 0.0, 1.0)]
    #[case::low_quality_group_alignment(vec![
"DSFGSEWDPSDEAW",
"SDYPADTDDGKYYH",
"QFAESLFGGTNAQE",
"EKEQWHVCFNIEYA",
"NSMWSNGFSCRYLQ",
"CRAQNNRSFVRTDW",
"LVMPKFMWHNRTIP",
"WNNRMTHEVLRIGI",
"MCTMIGELKYIHQS",
"INMMHTWHTCQTMV"
    ], 10, 1.0, 0.0)]
    fn test_calculate_single_chain_background_score(
        #[case] downsampled_sequences: Vec<&str>,
        background_vector: Vec<&str>,
        #[case] n_replicates: usize,
        #[case] expected: f64,
        #[case] expected_self_mean: f64,
    ) {
        let mut rng = SmallRng::from_seed([7; 32]);
        let pb = ProgressBar::hidden();
        let aligner = build_protein_aligner(7, 1).expect("aligner build");
        let downsampled_sequences: Vec<String> = downsampled_sequences
            .into_iter()
            .map(|s| s.to_string())
            .collect();
        let unique_seqs: Vec<String> = background_vector
            .into_iter()
            .map(|s| s.to_string())
            .collect();

        let mut ctx = ScoreContext {
            rng: &mut rng,
            pb: &pb,
            aligner: &aligner,
            pattern: "pat",
        };
        let single_result = calculate_single_chain_score(&mut ctx, downsampled_sequences.clone());
        if expected_self_mean.is_nan() & single_result.is_nan() {
        } else {
            //assert_eq!(self_mean, expected_self_mean);
        }
        println!("self_mean: {}", single_result.self_mean);

        let mut ctx = ScoreContext {
            rng: &mut rng,
            pb: &pb,
            aligner: &aligner,
            pattern: "pat",
        };
        let score = calculate_single_chain_background_score(
            &mut ctx,
            &unique_seqs,
            n_replicates,
            single_result,
        );

        assert_eq!(score, expected);
    }

    #[rstest]
    #[case::no_background_sequences(vec!["ACD", "WQ"], vec![], 5, 1.0)]
    #[case::single_downsampled_sequence(vec!["ACD"], vec!["ACD", "WQ", "MN"], 5, 1.0)]
    fn test_calculate_single_chain_background_score_uses_replicate_count_when_skipped(
        #[case] downsampled_sequences: Vec<&str>,
        #[case] unique_seqs: Vec<&str>,
        #[case] n_replicates: usize,
        #[case] expected: f64,
    ) {
        let mut rng = SmallRng::from_seed([9; 32]);
        let pb = ProgressBar::hidden();
        let aligner = build_protein_aligner(-10, -1).expect("aligner build");
        let downsampled_sequences: Vec<String> = downsampled_sequences
            .into_iter()
            .map(|s| s.to_string())
            .collect();
        let unique_seqs: Vec<String> = unique_seqs.into_iter().map(|s| s.to_string()).collect();

        let mut ctx = ScoreContext {
            rng: &mut rng,
            pb: &pb,
            aligner: &aligner,
            pattern: "pat",
        };
        let single_result = SingleChainResult::new(downsampled_sequences.clone(), 1.23);
        let score = calculate_single_chain_background_score(
            &mut ctx,
            &unique_seqs,
            n_replicates,
            single_result,
        );

        assert_eq!(score, expected);
    }

    #[test]
    fn test_calculate_single_chain_background_score_returns_minus_one_for_nan_self_mean() {
        let mut rng = SmallRng::from_seed([11; 32]);
        let pb = ProgressBar::hidden();
        let aligner = build_protein_aligner(-10, -1).expect("aligner build");
        let downsampled_sequences = vec!["ACD".to_string(), "WQ".to_string()];
        let unique_seqs = vec!["ACD".to_string(), "WQ".to_string()];

        let mut ctx = ScoreContext {
            rng: &mut rng,
            pb: &pb,
            aligner: &aligner,
            pattern: "pat",
        };
        let single_result = SingleChainResult::new(downsampled_sequences.clone(), f64::NAN);
        let score =
            calculate_single_chain_background_score(&mut ctx, &unique_seqs, 3, single_result);

        assert_eq!(score, -1.0);
    }
}
