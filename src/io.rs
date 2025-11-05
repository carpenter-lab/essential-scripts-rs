use polars::prelude::*;
use std::error::Error;
use std::path::PathBuf;

pub trait WriteToCsvOrStdout {
    fn write_to_csv_or_stdout(self, output_file: &PathBuf) -> Result<(), Box<dyn Error>>;
}

impl WriteToCsvOrStdout for DataFrame {
    fn write_to_csv_or_stdout(mut self, output_file: &PathBuf) -> Result<(), Box<dyn Error>> {
        if *output_file == PathBuf::from("-") {
            let stdout = std::io::stdout();
            let output = stdout.lock();
            let fin = CsvWriter::new(output).finish(&mut self);
            match fin {
                Ok(_) => Ok(()),
                Err(e) => Err(e.into()),
            }
        } else {
            let output = std::fs::File::create(output_file).expect("Failed to create output file");
            let fin = CsvWriter::new(output).finish(&mut self);
            match fin {
                Ok(_) => Ok(()),
                Err(e) => Err(e.into()),
            }
        }
    }
}

impl WriteToCsvOrStdout for LazyFrame {
    fn write_to_csv_or_stdout(self, output_file: &PathBuf) -> Result<(), Box<dyn Error>> {
        if *output_file == PathBuf::from("-") {
            let stdout = std::io::stdout();
            let output = stdout.lock();
            let fin = CsvWriter::new(output).finish(&mut self.clone().collect()?);
            match fin {
                Ok(_) => Ok(()),
                Err(e) => Err(e.into()),
            }
        } else {
            let output = std::fs::File::create(output_file).expect("Failed to create output file");
            let fin = CsvWriter::new(output).finish(&mut self.clone().collect()?);
            match fin {
                Ok(_) => Ok(()),
                Err(e) => Err(e.into()),
            }
        }
    }
}
