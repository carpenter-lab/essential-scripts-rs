use crate::io;
use crate::io::WriteToCsvOrStdout;
use clap::Subcommand;
use polars::prelude::*;
use std::path::PathBuf;

#[derive(Subcommand)]
pub enum Commands {
    #[command(about = "Split sample ID into subject and condition from GLIPH2 output")]
    SplitSampleId {
        #[arg(required = true, num_args = 1.., help = "Input CSV file to process")]
        input_file: PathBuf,

        #[arg(required = true, help = "Output file path")]
        output_file: PathBuf,

        #[arg(
            short,
            long,
            default_value = "subject:condition",
            help = "Column to split"
        )]
        column_name: String,
    },
    #[command(about = "Split CDR3 sequences and genes if a semicolon is present")]
    SplitCdr3Seq {
        #[arg(required = true, help = "Input CSV file to process")]
        input_file: PathBuf,

        #[arg(required = true, help = "Output file path")]
        output_file: PathBuf,

        #[arg(
            short,
            long,
            default_value = "subject:condition",
            help = "Columns to group by"
        )]
        group: Vec<String>,
    },
}

fn get_gene(gene_series: &Column, fct_series: &Series) -> PolarsResult<Series> {
    let gene_ca = gene_series.str()?;
    let fct_ca = fct_series.i64()?;

    let result: Vec<Option<String>> = gene_ca
        .into_iter()
        .zip(fct_ca.into_iter())
        .map(|(gene_opt, fct_opt)| match (gene_opt, fct_opt) {
            (Some(gene), Some(fct)) if gene.contains(';') => {
                let parts: Vec<&str> = gene.split(';').collect();
                parts.get(fct as usize).map(|s| s.to_string())
            }
            (Some(gene), _) => Some(gene.to_string()),
            _ => None,
        })
        .collect();

    Ok(Series::new(gene_series.name().clone(), result))
}

fn split_col<'a>(
    df: DataFrame,
    col_name: &PlSmallStr,
    fct_series: &Series,
) -> PolarsResult<DataFrame> {
    if !df.get_column_names().contains(&col_name) {
        return Ok(df);
    }

    let gene_series = df.column(col_name)?;
    let new_series = get_gene(gene_series, fct_series)?;

    df.drop(col_name)?.with_column(new_series.into()).cloned()
}

fn factor_groups(df: &LazyFrame, group: Vec<&str>) -> Series {
    let bind = df
        .clone()
        .with_column(
            Expr::from(cols(group))
                .str()
                .join("_", false)
                .alias("_group_key"),
        )
        .collect()
        .unwrap();
    let s = bind.column("_group_key").unwrap().as_series().unwrap();
    let categories: Arc<Categories> = s.unique().unwrap().0.as_arc_any().downcast().unwrap();
    let categories_m: Arc<CategoricalMapping> =
        s.unique().unwrap().0.as_arc_any().downcast().unwrap();
    s.cast(&DataType::Categorical(categories.clone(), categories_m))
        .unwrap()
}

fn split_cdr3_seq(df: LazyFrame, group: &[&str], chain: &str) -> PolarsResult<LazyFrame> {
    if chain != "alpha" && chain != "beta" {
        return Err(PolarsError::ComputeError(
            format!("chain must be alpha or beta, not {}", chain).into(),
        ));
    }

    let (cdr3_col, gene) = if chain == "alpha" {
        ("CDR3a", "A")
    } else {
        ("CDR3b", "B")
    };

    let cols_to_split: Vec<PlSmallStr> = vec![
        format!("TR{}V", gene).into(),
        format!("TR{}J", gene).into(),
        cdr3_col.to_string().into(),
    ];

    // Build the complete group columns list
    let mut all_group_cols: Vec<&str> = group.to_vec();
    let cols_to_split_refs: Vec<&str> = cols_to_split.iter().map(|s| s.as_str()).collect();
    all_group_cols.extend(&cols_to_split_refs);

    // Calculate factorization
    let fct = factor_groups(&df, all_group_cols);

    // Split each column
    let mut df2 = df.clone().collect()?;
    for col_name in &cols_to_split {
        df2 = split_col(df2, col_name, &fct)?;
    }

    Ok(df2.lazy())
}

pub(crate) fn split_cdr3_seq_main(
    input_file: PathBuf,
    output_file: PathBuf,
    group: Vec<String>,
) -> () {
    let mut df: LazyFrame = io::read_from_file(input_file, None);

    // Convert group strings to string slices
    let group_refs: Vec<&str> = group.iter().map(|s| s.as_str()).collect();

    // Process both chains
    for chain in &["alpha", "beta"] {
        df = split_cdr3_seq(df, &group_refs, chain).expect("Failed to split CDR3 sequences");
    }

    // Write output
    df.collect()
        .expect("Failed to collect dataframe")
        .write_to_flat_or_stdout(output_file, None)
}

pub(crate) fn split_sample_id(
    input_file: PathBuf,
    output_file: PathBuf,
    column_name: &String,
) -> () {
    let df = io::read_from_file(input_file, None);
    let df = df
        .with_column(col(column_name).str().split(lit(":")))
        .with_columns([
            col(column_name).list().first().alias("subject"),
            col(column_name).list().last().alias("condition"),
        ])
        .select([all().exclude_cols([column_name]).as_expr()]);

    df.collect()
        .expect("Failed to collect dataframe")
        .write_to_flat_or_stdout(output_file, None)
}

pub fn handle_command(cmd: Commands) -> () {
    match cmd {
        Commands::SplitSampleId {
            input_file,
            output_file,
            column_name,
        } => {
            split_sample_id(input_file, output_file, &column_name);
        }
        Commands::SplitCdr3Seq {
            input_file,
            output_file,
            group,
        } => {
            split_cdr3_seq_main(input_file, output_file, group);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_file_path(prefix: &str, suffix: &str) -> PathBuf {
        let mut path = std::env::temp_dir();
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        path.push(format!(
            "{}_{}_{}{}",
            prefix,
            std::process::id(),
            nanos,
            suffix
        ));
        path
    }

    #[test]
    fn get_gene_splits_by_factor_index() {
        let gene_series = Series::new("TRAV".into(), &[Some("A;B"), Some("C;D"), None]);
        let gene_col: Column = gene_series.into();
        let fct_series = Series::new("fct".into(), &[0i64, 1i64, 0i64]);

        let out = get_gene(&gene_col, &fct_series).unwrap();
        let out_ca = out.str().unwrap();

        assert_eq!(out_ca.get(0), Some("A"));
        assert_eq!(out_ca.get(1), Some("D"));
        assert_eq!(out_ca.get(2), None);
    }

    #[test]
    fn split_col_noop_when_missing() {
        let df =
            DataFrame::new_infer_height(vec![Series::new("foo".into(), &["x", "y"]).into_column()])
                .unwrap();
        let fct_series = Series::new("fct".into(), &[0i64, 1i64]);

        let out = split_col(df.clone(), &PlSmallStr::from("bar"), &fct_series).unwrap();

        assert_eq!(out.get_column_names(), df.get_column_names());
        assert_eq!(out.height(), df.height());
    }

    #[test]
    fn split_col_replaces_values() {
        let df = DataFrame::new_infer_height(vec![
            Series::new("TRAV".into(), &["AV1;AV2", "AV3;AV4"]).into_column(),
        ])
        .unwrap();
        let fct_series = Series::new("fct".into(), &[0i64, 1i64]);

        let out = split_col(df, &PlSmallStr::from("TRAV"), &fct_series).unwrap();
        let out_ca = out.column("TRAV").unwrap().str().unwrap();

        assert_eq!(out_ca.get(0), Some("AV1"));
        assert_eq!(out_ca.get(1), Some("AV4"));
        assert!(out_ca.into_iter().flatten().all(|v| !v.contains(';')));
    }

    #[test]
    fn split_cdr3_seq_rejects_invalid_chain() {
        let df = DataFrame::new_infer_height(vec![
            Series::new("subject".into(), &["s1", "s2"]).into_column(),
            Series::new("TRAV".into(), &["AV1;AV2", "AV3;AV4"]).into_column(),
            Series::new("TRAJ".into(), &["AJ1;AJ2", "AJ3;AJ4"]).into_column(),
            Series::new("CDR3a".into(), &["CAA;CAB", "CCA;CCD"]).into_column(),
        ])
        .unwrap();
        let lf = df.lazy();

        match split_cdr3_seq(lf, &["subject"], "gamma") {
            Err(PolarsError::ComputeError(_)) => {}
            Err(_) => panic!("unexpected error type"),
            Ok(_) => panic!("expected error for invalid chain"),
        }
    }

    #[test]
    fn split_sample_id_end_to_end() {
        let input_path = temp_file_path("split_sample_id_input", ".csv");
        let output_path = temp_file_path("split_sample_id_output", ".csv");

        let csv = "sample\nsub1:condA\nsub2:condB\n";
        fs::write(&input_path, csv).unwrap();

        split_sample_id(input_path, output_path.clone(), &"sample".to_string());

        let out_df = io::read_from_csv(output_path).collect().unwrap();
        let columns = out_df.get_column_names();

        assert!(columns.iter().any(|name| name.as_str() == "subject"));
        assert!(columns.iter().any(|name| name.as_str() == "condition"));
        assert!(!columns.iter().any(|name| name.as_str() == "sample"));

        let subject = out_df.column("subject").unwrap().str().unwrap();
        let condition = out_df.column("condition").unwrap().str().unwrap();
        assert_eq!(subject.get(0), Some("sub1"));
        assert_eq!(condition.get(0), Some("condA"));
        assert_eq!(subject.get(1), Some("sub2"));
        assert_eq!(condition.get(1), Some("condB"));
    }
}
