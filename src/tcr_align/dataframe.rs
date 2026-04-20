use crate::tcr_align::align::{build_protein_aligner, run_parasail};
use indicatif::ProgressStyle;
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
                keep_nulls: false,
            },
        )
        .select([col("TcRa_list")])
        .filter(col("TcRa_list").is_not_null())
        .unique(None, UniqueKeepStrategy::First)
        .collect()?;

    let s = unique_df.column("TcRa_list")?.str()?;
    // Skip nulls and empty/whitespace-only strings produced by splitting
    Ok(s.into_iter()
        .flatten()
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect())
}

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
                keep_nulls: false,
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
                keep_nulls: false,
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

    for x in s.str()?.into_iter().flatten() {
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

pub(crate) fn fraction_self_greater(
    groups: &DataFrame,
    all_unique_alpha: &[String],
    n_replicates: usize,
    gap_open: i32,
    gap_extend: i32,
) -> PolarsResult<DataFrame> {
    let height = &groups.height();

    let patterns: Vec<String> = groups
        .column("pattern")?
        .as_series()
        .expect("could not get 'pattern' column")
        .str()?
        .into_no_null_iter()
        .map(|s| s.to_string())
        .collect();

    let mut out_pattern: Vec<String> = Vec::with_capacity(*height);
    let mut out_self_a: Vec<f64> = Vec::with_capacity(*height);
    let mut out_self_b: Vec<f64> = Vec::with_capacity(*height);
    let mut out_frac: Vec<f64> = Vec::with_capacity(*height);

    let mut rng = SmallRng::from_rng(&mut rand::rng());

    calculate_scores(
        groups,
        all_unique_alpha,
        n_replicates,
        &gap_open,
        &gap_extend,
        *height,
        patterns,
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
    pb: &indicatif::ProgressBar,
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

fn calculate_scores(
    groups: &DataFrame,
    all_unique_alpha: &[String],
    n_replicates: usize,
    gap_open: &i32,
    gap_extend: &i32,
    height: usize,
    patterns: Vec<String>,
    out_pattern: &mut Vec<String>,
    out_self_a: &mut Vec<f64>,
    out_self_b: &mut Vec<f64>,
    out_frac: &mut Vec<f64>,
    rng: &mut SmallRng,
) {
    let progress_style = ProgressStyle::with_template(
        "{msg} [{bar:40.cyan/blue}] {pos}/{len} Elapsed: {elapsed_precise} ETA: {eta}",
    );
    let inner_progress_style =
        ProgressStyle::with_template("\x1b[37m{msg}\x1b[0m [{bar:40.cyan/blue}] {pos}/{len}");
    let mpb: indicatif::MultiProgress;
    let pb: indicatif::ProgressBar;
    let inner_pb: indicatif::ProgressBar;

    match progress_style {
        Ok(progress_style) => {
            mpb = indicatif::MultiProgress::new();
            pb = mpb.add(indicatif::ProgressBar::new(height as u64));
            pb.set_style(progress_style);
            match inner_progress_style {
                Ok(inner_progress_style) => {
                    inner_pb = mpb.add(indicatif::ProgressBar::new(n_replicates as u64));
                    inner_pb.set_style(inner_progress_style);
                }
                Err(_) => {
                    inner_pb = indicatif::ProgressBar::hidden();
                }
            }
        }
        Err(_) => {
            pb = indicatif::ProgressBar::hidden();
            inner_pb = indicatif::ProgressBar::hidden();
        }
    }

    let aligner = match build_protein_aligner(*gap_open, *gap_extend) {
        Ok(aligner) => aligner,
        Err(e) => {
            eprintln!("Error building protein aligner: {}", e);
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

        match get_list_cell_as_vec_utf8(&groups, "TcRb", i) {
            Err(_) => {
                out_self_b.push(f64::NAN);
            }
            Ok(b_full) => {
                let b = downsample_vec(b_full, pattern, &inner_pb, rng);
                let self_b_mean = if b.len() <= 1 {
                    f64::NAN
                } else {
                    match run_parasail(&aligner, &b, None) {
                        Ok(mean) => mean,
                        Err(e) => {
                            eprintln!("Error running parasail for pattern {}: {}", pattern, e);
                            f64::NAN
                        }
                    }
                };
                out_self_b.push(self_b_mean);
            }
        }

        match get_list_cell_as_vec_utf8(&groups, "TcRa", i) {
            Err(_) => {
                out_self_a.push(f64::NAN);
            }
            Ok(a_full) => {
                let a = downsample_vec(a_full, pattern, &inner_pb, rng);
                let self_a_mean = if a.len() <= 1 {
                    f64::NAN
                } else {
                    match run_parasail(&aligner, &a, None) {
                        Ok(mean) => mean,
                        Err(e) => {
                            eprintln!("Error running parasail for pattern {}: {}", pattern, e);
                            f64::NAN
                        }
                    }
                };
                out_self_a.push(self_a_mean);
                let n_val = a.len().min(all_unique_alpha.len()).min(MAX_BG_SEQS);
                let mut greater_count = 0usize;
                if a.len() > 1 && n_val > 0 && !self_a_mean.is_nan() {
                    for _rep in 0..n_replicates {
                        let background: Vec<String> =
                            all_unique_alpha.sample(rng, n_val).cloned().collect();

                        let bg_mean = match run_parasail(&aligner, &a, Some(&background)) {
                            Ok(mean) => mean,
                            Err(e) => {
                                eprintln!("Error running parasail for pattern {}: {}", pattern, e);
                                f64::NAN
                            }
                        };
                        // we want lower scores to be better
                        if self_a_mean < bg_mean {
                            greater_count += 1;
                        }
                        inner_pb.inc(1);
                    }
                } else {
                    greater_count = n_replicates;
                }
                out_frac.push((greater_count as f64) / (n_replicates as f64));
            }
        }

        inner_pb.reset();
        out_pattern.push(pattern.to_string());
        pb.inc(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use indicatif::ProgressBar;
    use polars::prelude::*;
    use rand::SeedableRng;
    use rand::rngs::SmallRng;
    use std::collections::HashSet;

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
            .into_no_null_iter()
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

    #[test]
    fn test_downsample_vec_limits() {
        // create > MAX_GROUP_SEQS items
        let total = MAX_GROUP_SEQS + 200;
        let v: Vec<String> = (0..total).map(|i| format!("S{}", i)).collect();
        let mut rng = SmallRng::from_seed([1; 32]);
        let pb = ProgressBar::hidden();

        let out = downsample_vec(v.clone(), "pat", &pb, &mut rng);
        assert_eq!(out.len(), MAX_GROUP_SEQS);

        // small vector should be returned unchanged
        let small = vec!["a".to_string()];
        let mut rng2 = SmallRng::from_seed([1; 32]);
        let out2 = downsample_vec(small.clone(), "pat", &pb, &mut rng2);
        assert_eq!(out2, small);
    }
}
