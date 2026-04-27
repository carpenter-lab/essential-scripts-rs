mod align;
mod dataframe;

use crate::io::WriteToCsvOrStdout;
use clap::Subcommand;
use polars::prelude::*;
use std::path::PathBuf;

#[derive(Subcommand)]
pub enum Commands {
    #[command(about = "Score a GLIPH2 output file with a TCR alignment pipeline")]
    ScoreTCRAlignments {
        #[arg(required = true, help = "Input CSV file to process")]
        input_file: PathBuf,

        #[arg(required = true, help = "Output file path")]
        output_file: PathBuf,

        #[arg(short, long, default_value_t = 1000)]
        replicates: usize,
    },
}

pub trait Run {
    type Output;
    fn run(&self, n_replicates: usize) -> PolarsResult<Self::Output>;
}

impl Run for DataFrame {
    type Output = DataFrame;
    fn run(self: &DataFrame, n_replicates: usize) -> PolarsResult<Self::Output> {
        let all_unique_alpha = dataframe::all_unique_cdr3_alpha(self)?;

        //let groups = dataframe::prepare_parasail_groups(self);
        match dataframe::prepare_parasail_groups(self) {
            Err(e) => {
                eprintln!("Error preparing parasail groups: {}", e);
                Err(e)
            }
            Ok(groups) => {
                dataframe::fraction_self_greater(&groups, &all_unique_alpha, n_replicates, 7, 1)
            }
        }
    }
}

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

pub(crate) fn tcr_score(input_file: PathBuf, output_file: PathBuf, replicates: usize) -> () {
    let df = LazyCsvReader::new(PlRefPath::try_from_pathbuf(input_file).unwrap())
        .finish()
        .unwrap()
        .collect()
        .unwrap();

    let res = score_df(&df, replicates).unwrap();
    res.write_to_csv_or_stdout(output_file)
}

pub fn handle_command(cmd: Commands) -> () {
    match cmd {
        Commands::ScoreTCRAlignments {
            input_file,
            output_file,
            replicates,
        } => {
            tcr_score(input_file, output_file, replicates);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
