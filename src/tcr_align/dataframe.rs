use crate::tcr_align::align::{build_protein_aligner, run_parasail};
use indicatif::ProgressStyle;
use polars::error::PolarsResult;
use polars::frame::{DataFrame, UniqueKeepStrategy};
use polars::prelude::{IntoLazy, NamedFrom, Series, col, len, lit};
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
        .explode(col("TcRa_list").into_selector().unwrap())
        .select([col("TcRa_list")])
        .filter(col("TcRa_list").is_not_null())
        .unique(None, UniqueKeepStrategy::First)
        .collect()?;

    let s = unique_df.column("TcRa_list")?.str()?;
    Ok(s.into_no_null_iter().map(|s| s.to_string()).collect())
}

pub(crate) fn prepare_parasail_groups(all_data: &DataFrame) -> DataFrame {
    let lf = all_data
        .clone()
        .lazy()
        .select([col("pattern"), col("TcRa"), col("TcRb")])
        .unique(None, UniqueKeepStrategy::First)
        // TcRa -> explode
        .with_columns([col("TcRa").str().split(lit(";")).alias("TcRa_list")])
        .explode(col("TcRa_list").into_selector().unwrap())
        .drop(col("TcRa").into_selector().unwrap())
        .rename(["TcRa_list"], ["TcRa"], true)
        .filter(col("TcRa").is_not_null())
        // TcRb -> explode
        .with_columns([col("TcRb").str().split(lit(";")).alias("TcRb_list")])
        .explode(col("TcRb_list").into_selector().unwrap())
        .drop(col("TcRb").into_selector().unwrap())
        .rename(["TcRb_list"], ["TcRb"], true)
        .filter(col("TcRb").is_not_null())
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
        .filter(col("n").gt(lit(1)));

    grouped.collect().unwrap()
}

fn get_list_cell_as_vec_utf8(
    df: &DataFrame,
    col_name: &str,
    row_idx: usize,
) -> PolarsResult<Vec<String>> {
    let s = df.column(col_name)?.list()?.get_as_series(row_idx).unwrap();
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
        .unwrap()
        .str()?
        .into_no_null_iter()
        .map(|s| s.to_string())
        .collect();

    let mut out_pattern: Vec<String> = Vec::with_capacity(*height);
    let mut out_self_a: Vec<f64> = Vec::with_capacity(*height);
    let mut out_self_b: Vec<f64> = Vec::with_capacity(*height);
    let mut out_frac: Vec<f64> = Vec::with_capacity(*height);

    let mut rng = SmallRng::from_rng(&mut rand::rng());

    let progress_style = ProgressStyle::with_template(
        "{msg} [{bar:40.cyan/blue}] {pos}/{len} Elapsed: {elapsed_precise} ETA: {eta}",
    )
    .unwrap();
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
        progress_style,
    );

    DataFrame::new(vec![
        Series::new("pattern".into(), out_pattern).into(),
        Series::new("TcRa_alignment_score".into(), out_self_a).into(),
        Series::new("TcRb_alignment_score".into(), out_self_b).into(),
        Series::new("TcRa_alignment_score_background".into(), out_frac).into(),
    ])
}

fn downsample_vec(
    v: Vec<String>,
    pattern: &str,
    pb: &indicatif::ProgressBar,
    mut rng: &mut SmallRng,
) -> Vec<String> {
    let out: Vec<String> = if v.len() > MAX_GROUP_SEQS {
        pb.set_message(format!(
            "Downsampling {p} {full} -> {maxg}",
            p = pattern,
            full = v.len(),
            maxg = MAX_GROUP_SEQS
        ));
        v.choose_multiple(&mut rng, MAX_GROUP_SEQS)
            .cloned()
            .collect()
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
    mut rng: &mut SmallRng,
    progress_style: ProgressStyle,
) {
    let mpb = indicatif::MultiProgress::new();

    let pb = mpb.add(indicatif::ProgressBar::new(height as u64));
    pb.set_style(progress_style.clone());
    let inner_pb = mpb.add(indicatif::ProgressBar::new(n_replicates as u64));
    inner_pb.set_style(
        ProgressStyle::with_template("\x1b[37m{msg}\x1b[0m [{bar:40.cyan/blue}] {pos}/{len}")
            .unwrap(),
    );

    let aligner = build_protein_aligner(*gap_open, *gap_extend).unwrap();
    for i in 0..height {
        inner_pb.set_message(format!("{:width$}", "", width = 33).to_string());
        let pattern = patterns.get(i).map_or("", |v| v);
        if pattern == "single" {
            pb.inc(1);
            out_pattern.push(pattern.to_string());
            out_self_a.push(f64::NAN);
            out_self_b.push(f64::NAN);
            out_frac.push(f64::NAN);
            continue;
        }
        pb.set_message(
            format!("Calculating scores for {:width$}", pattern, width = 10).to_string(),
        );

        let a_full = get_list_cell_as_vec_utf8(&groups, "TcRa", i).unwrap();
        let a = downsample_vec(a_full, pattern, &inner_pb, rng);

        let self_a_mean = if a.len() <= 1 {
            f64::NAN
        } else {
            run_parasail(&aligner, &a, None).unwrap()
        };

        let b_full = get_list_cell_as_vec_utf8(&groups, "TcRb", i).unwrap();
        let b: Vec<String> = downsample_vec(b_full, pattern, &inner_pb, rng);

        let self_b_mean = if b.len() <= 1 {
            f64::NAN
        } else {
            run_parasail(&aligner, &b, None).unwrap()
        };

        // Keep background size bounded too.
        let n_val = a.len().min(all_unique_alpha.len()).min(MAX_BG_SEQS);

        let mut greater_count = 0usize;
        if a.len() > 1 && n_val > 0 {
            for _rep in 0..n_replicates {
                let background: Vec<String> = all_unique_alpha
                    .choose_multiple(&mut rng, n_val)
                    .cloned()
                    .collect();

                let bg_mean = run_parasail(&aligner, &a, Some(&background)).unwrap();

                if self_a_mean > bg_mean {
                    greater_count += 1;
                }
                inner_pb.inc(1);
            }
        }
        inner_pb.reset();

        out_pattern.push(pattern.to_string());
        out_self_a.push(self_a_mean);
        out_self_b.push(self_b_mean);
        out_frac.push((greater_count as f64) / (n_replicates as f64));
        pb.inc(1);
    }
}
