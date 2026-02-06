use calamine::{Reader, Xlsx, open_workbook};
use polars::prelude::*;
use std::path::{Path, PathBuf};

fn get_extension_from_filename(filename: &PathBuf) -> Option<&str> {
    Path::new(filename).extension().and_then(|s| s.to_str())
}

fn write_df_output(df: &mut DataFrame, output_file: &PathBuf, separator: u8) {
    if output_file == Path::new("-") {
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
    if output_file == Path::new("-") {
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

pub fn read_from_csv(input_file: &PathBuf) -> LazyFrame {
    match LazyCsvReader::new(PlPath::new(input_file.to_str().unwrap())).finish() {
        Ok(lf) => lf,
        Err(e) => panic!("Failed to read CSV file: {}", e),
    }
}

pub fn read_from_tsv(input_file: &PathBuf) -> LazyFrame {
    match LazyCsvReader::new(PlPath::new(input_file.to_str().unwrap())).finish() {
        Ok(lf) => lf,
        Err(e) => panic!("Failed to read TSV file: {}", e),
    }
}

pub fn read_from_file(input_file: &PathBuf, separator: Option<u8>) -> LazyFrame {
    match separator {
        Some(sep) => match sep {
            b'\t' => read_from_tsv(input_file),
            b',' => read_from_csv(input_file),
            _ => panic!("Unsupported separator"),
        },
        _ => match get_extension_from_filename(input_file) {
            Some(v) => match v {
                "tsv" => read_from_tsv(input_file),
                "xls" | "xlsx" => read_excel(input_file, None, 0, None)
                    .expect("Failed to read Excel file")
                    .lazy(),
                _ => read_from_csv(input_file),
            },
            None => read_from_tsv(input_file),
        },
    }
}

pub fn strip_quotes(s: &str) -> String {
    if (s.starts_with('"') && s.ends_with('"')) || (s.starts_with('\'') && s.ends_with('\'')) {
        s[1..s.len() - 1].to_string()
    } else {
        s.to_string()
    }
}

pub fn read_excel(
    input_path: &PathBuf,
    sheet_name: Option<&str>,
    skiprows: usize,
    header: Option<bool>,
) -> PolarsResult<DataFrame> {
    // Read Excel file using calamine
    let mut workbook: Xlsx<_> = open_workbook(input_path).map_err(|e| {
        PolarsError::ComputeError(format!("Failed to open Excel file: {}", e).into())
    })?;

    let sheet_names = workbook.sheet_names();
    let sheet_name = match sheet_name {
        Some(name) => name,
        None => sheet_names
            .get(0)
            .ok_or_else(|| PolarsError::ComputeError("Excel file contains no sheets".into()))?
            .as_str(),
    };

    // Get the specified sheet
    let range = workbook.worksheet_range(sheet_name).map_err(|e| {
        PolarsError::ComputeError(format!("Failed to get sheet '{}': {}", sheet_name, e).into())
    })?;

    // Convert the range to a DataFrame
    let (height, width) = range.get_size();
    let col_names = match range.headers() {
        Some(headers) => match header {
            Some(true) | None => headers
                .iter()
                .map(|s| strip_quotes(&s.to_string()))
                .collect::<Vec<String>>(),
            Some(false) => (0..width)
                .map(|i| format!("column_{}", i))
                .collect::<Vec<String>>(),
        },
        None => (0..width)
            .map(|i| format!("column_{}", i))
            .collect::<Vec<String>>(),
    };
    // check if repeated values are in col_names
    if col_names.len()
        != col_names
            .iter()
            .collect::<std::collections::HashSet<_>>()
            .len()
    {
        return Err(PolarsError::ComputeError(
            "Duplicate column names found in Excel sheet".into(),
        ));
    }

    // Extract all data as strings first
    let mut columns: Vec<Vec<Option<String>>> = vec![Vec::new(); width];

    for row_idx in skiprows..height {
        if let Some(_h) = range.headers() {
            if row_idx == skiprows {
                continue; // Skip header row
            }
        };
        for col_idx in 0..width {
            let cell_value = range.get((row_idx, col_idx));
            let value = cell_value.map(|v| v.to_string());
            columns[col_idx].push(value);
        }
    }

    // Create Series from columns
    let mut series_vec: Vec<Column> = Vec::new();
    for (col_idx, col_data) in columns.into_iter().enumerate() {
        let col_name = &col_names[col_idx];
        let s = Series::new(PlSmallStr::from_str(&col_name), col_data);
        series_vec.push(s.into_column());
    }

    let df = DataFrame::new(series_vec);
    df
}
