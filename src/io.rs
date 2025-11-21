use polars::prelude::*;
use std::ffi::OsStr;
use std::path::{Path, PathBuf};

fn get_extension_from_filename(filename: &PathBuf) -> Option<&str> {
    Path::new(filename).extension().and_then(OsStr::to_str)
}

fn write_df_output(df: &mut DataFrame, output_file: &PathBuf, separator: u8) {
    if output_file == std::path::Path::new("-") {
        let stdout = std::io::stdout();
        let output = stdout.lock();
        CsvWriter::new(output)
            .with_separator(separator)
            .finish(df)
            .expect("Failed to write to stdout");
    } else {
        let output = std::fs::File::create(output_file).expect("Failed to create output file");
        CsvWriter::new(output)
            .with_separator(separator)
            .finish(df)
            .expect("Failed to write to file");
    }
}

fn write_lazy_output(lf: LazyFrame, output_file: &PathBuf, separator: u8) {
    if output_file == std::path::Path::new("-") {
        let mut df = lf
            .collect()
            .expect("Failed to collect lazy frame prior to writing to stdout");
        write_df_output(&mut df, output_file, separator);
    } else {
        let mut writer_opts = CsvWriterOptions::default();
        writer_opts.serialize_options.separator = separator;
        lf.sink_csv(
            SinkTarget::Path(PlPath::new(output_file.to_str().unwrap())),
            writer_opts,
            None,
            SinkOptions::default(),
        )
        .expect("Failed to open to CSV file for writing")
        .collect()
        .expect("Failed to collect lazy frame prior to writing to file");
    }
}

pub trait WriteToCsvOrStdout {
    fn write_to_csv_or_stdout(self, output_file: &PathBuf) -> ();
    fn write_to_tsv_or_stdout(self, output_file: &PathBuf) -> ();
    fn write_to_flat_or_stdout(self, output_file: &PathBuf, separator: Option<u8>) -> ()
    where
        Self: Sized,
    {
        match separator {
            Some(sep) => match sep {
                b'\t' => self.write_to_tsv_or_stdout(output_file),
                b',' => self.write_to_csv_or_stdout(output_file),
                _ => panic!("Unsupported separator"),
            },
            _ => match get_extension_from_filename(output_file) {
                Some(v) => match v {
                    "tsv" => self.write_to_tsv_or_stdout(output_file),
                    _ => self.write_to_csv_or_stdout(output_file),
                },
                None => self.write_to_tsv_or_stdout(output_file),
            },
        }
    }
}

impl WriteToCsvOrStdout for DataFrame {
    fn write_to_csv_or_stdout(mut self, output_file: &PathBuf) -> () {
        write_df_output(&mut self, output_file, b',');
    }
    fn write_to_tsv_or_stdout(mut self, output_file: &PathBuf) -> () {
        write_df_output(&mut self, output_file, b'\t');
    }
}
impl WriteToCsvOrStdout for LazyFrame {
    fn write_to_csv_or_stdout(self, output_file: &PathBuf) -> () {
        write_lazy_output(self, output_file, b',');
    }
    fn write_to_tsv_or_stdout(self, output_file: &PathBuf) -> () {
        write_lazy_output(self, output_file, b'\t');
    }
}
