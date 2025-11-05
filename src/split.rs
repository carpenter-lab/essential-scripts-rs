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

    df.drop(col_name)?.with_column(new_series).cloned()
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
    input_file: &PathBuf,
    output_file: &PathBuf,
    group: Vec<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut df: LazyFrame = LazyCsvReader::new(PlPath::from_str(input_file.to_str().expect("")))
        .finish()
        .unwrap();

    // Convert group strings to string slices
    let group_refs: Vec<&str> = group.iter().map(|s| s.as_str()).collect();

    // Process both chains
    for chain in &["alpha", "beta"] {
        df = split_cdr3_seq(df, &group_refs, chain)?;
    }

    // Write output
    let mut output = std::fs::File::create(output_file)?;
    let mut df = df.collect()?;
    let fin = CsvWriter::new(&mut output).finish(&mut df);

    match fin {
        Ok(_) => Ok(()),
        Err(e) => Err(e.into()),
    }
}

pub(crate) fn split_sample_id(
    input_file: &PathBuf,
    output_file: &PathBuf,
) -> Result<(), Box<dyn std::error::Error>> {
    let df = LazyCsvReader::new(PlPath::from_str(input_file.to_str().expect("")))
        .finish()
        .unwrap();
    let df = df.with_column(col("subject:condition").str().split(lit(":")));
    let output = std::fs::File::create(output_file).expect("Failed to create output file");
    let fin = CsvWriter::new(output).finish(&mut df.collect()?);

    match fin {
        Ok(_) => Ok(()),
        Err(e) => Err(e.into()),
    }
}

pub fn handle_command(cmd: Commands) -> () {
    match cmd {
        Commands::SplitSampleId {
            input_file,
            output_file,
        } => {
            split_sample_id(&input_file, &output_file).unwrap();
        }
        Commands::SplitCdr3Seq {
            input_file,
            output_file,
            group,
        } => {
            split_cdr3_seq_main(&input_file, &output_file, group).unwrap();
        }
    }
}
