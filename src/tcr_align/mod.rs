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
        #[arg(required = true, num_args = 1.., help = "Input CSV file to process")]
        input_file: String,

        #[arg(required = true, help = "Output file path")]
        output_file: PathBuf,

        #[arg(short, long, default_value_t = 1000)]
        replicates: usize,
    },
}

pub trait Run {
    type Output;
    fn run(&self, n_replicates: usize) -> PolarsResult<DataFrame>;
}

impl Run for DataFrame {
    type Output = DataFrame;
    fn run(self: &DataFrame, n_replicates: usize) -> PolarsResult<DataFrame> {
        let all_unique_alpha = dataframe::all_unique_cdr3_alpha(self)?;

        let groups = dataframe::prepare_parasail_groups(self);

        let result =
            dataframe::fraction_self_greater(&groups, &all_unique_alpha, n_replicates, 7, 1)?;

        Ok(result)
    }
}

pub(crate) fn tcr_score(input_file: &String, output_file: &PathBuf, replicates: usize) -> () {
    let df = LazyCsvReader::new(PlPath::new(input_file))
        .finish()
        .unwrap()
        .collect()
        .unwrap();

    let res = df.run(replicates).unwrap();
    let res = df
        .join(
            &res,
            ["pattern"],
            ["pattern"],
            JoinArgs::new(JoinType::Right),
            None,
        )
        .unwrap();
    res.write_to_csv_or_stdout(output_file)
}

pub fn handle_command(cmd: Commands) -> () {
    match cmd {
        Commands::ScoreTCRAlignments {
            input_file,
            output_file,
            replicates,
        } => {
            tcr_score(&input_file, &output_file, replicates);
        }
    }
}
