use polars::prelude::*;
use std::path::PathBuf;

pub trait WriteToCsvOrStdout {
    fn write_to_csv_or_stdout(self, output_file: &PathBuf) -> ();
}

impl WriteToCsvOrStdout for DataFrame {
    fn write_to_csv_or_stdout(mut self, output_file: &PathBuf) -> () {
        if *output_file == PathBuf::from("-") {
            let stdout = std::io::stdout();
            let output = stdout.lock();
            CsvWriter::new(output)
                .finish(&mut self)
                .expect("Failed to write to stdout");
        } else {
            let output = std::fs::File::create(output_file).expect("Failed to create output file");
            CsvWriter::new(output)
                .finish(&mut self)
                .expect("Failed to write to file");
        }
    }
}

impl WriteToCsvOrStdout for LazyFrame {
    fn write_to_csv_or_stdout(self, output_file: &PathBuf) -> () {
        if *output_file == PathBuf::from("-") {
            let stdout = std::io::stdout();
            let output = stdout.lock();
            CsvWriter::new(output)
                .finish(
                    &mut self
                        .clone()
                        .collect()
                        .expect("Failed to collect lazy frame prior to writing to stdout"),
                )
                .expect("Failed to write to stdout");
        } else {
            self.sink_csv(
                SinkTarget::Path(PlPath::new(output_file.to_str().unwrap())),
                CsvWriterOptions::default(),
                None,
                SinkOptions::default(),
            )
            .expect("Failed to open to CSV file for writing")
            .collect()
            .expect("Failed to collect lazy frame prior to writing to file");
        }
    }
}
