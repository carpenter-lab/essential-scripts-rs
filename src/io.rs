use calamine::{Reader, Xlsx, open_workbook};
use polars::prelude::*;
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::sync::Arc;

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

fn write_lazy_output(lf: LazyFrame, output_file: PathBuf, separator: u8) {
    if output_file == Path::new("-") {
        let mut df = lf
            .collect()
            .expect("Failed to collect lazy frame prior to writing to stdout");
        write_df_output(&mut df, &output_file, separator);
    } else {
        let writer_opts = CsvWriterOptions {
            include_bom: false,
            compression: ExternalCompression::default(),
            check_extension: false,
            include_header: false,
            batch_size: NonZero::new(1024).unwrap(),
            serialize_options: Arc::from(SerializeOptions {
                date_format: None,
                time_format: None,
                datetime_format: None,
                float_scientific: None,
                float_precision: None,
                decimal_comma: false,
                separator,
                quote_char: 0,
                null: Default::default(),
                line_terminator: PlSmallStr::from_str("\n"),
                quote_style: Default::default(),
            }),
        };
        lf.sink(
            SinkDestination::File {
                target: SinkTarget::Path(PlRefPath::try_from_pathbuf(output_file).unwrap()),
            },
            FileWriteFormat::Csv(writer_opts),
            UnifiedSinkArgs::default(),
        )
        .expect("Failed to open to CSV file for writing")
        .collect()
        .expect("Failed to collect lazy frame prior to writing to file");
    }
}

pub trait WriteToCsvOrStdout {
    fn write_to_csv_or_stdout(self, output_file: PathBuf) -> ();
    fn write_to_tsv_or_stdout(self, output_file: PathBuf) -> ();
    fn write_to_flat_or_stdout(self, output_file: PathBuf, separator: Option<u8>) -> ()
    where
        Self: Sized,
    {
        match separator {
            Some(sep) => match sep {
                b'\t' => self.write_to_tsv_or_stdout(output_file),
                b',' => self.write_to_csv_or_stdout(output_file),
                _ => panic!("Unsupported separator"),
            },
            _ => match get_extension_from_filename(&output_file) {
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
    fn write_to_csv_or_stdout(mut self, output_file: PathBuf) -> () {
        write_df_output(&mut self, &output_file, b',');
    }
    fn write_to_tsv_or_stdout(mut self, output_file: PathBuf) -> () {
        write_df_output(&mut self, &output_file, b'\t');
    }
}
impl WriteToCsvOrStdout for LazyFrame {
    fn write_to_csv_or_stdout(self, output_file: PathBuf) -> () {
        write_lazy_output(self, output_file, b',');
    }
    fn write_to_tsv_or_stdout(self, output_file: PathBuf) -> () {
        write_lazy_output(self, output_file, b'\t');
    }
}

pub fn read_from_csv(input_file: PathBuf) -> LazyFrame {
    match LazyCsvReader::new(PlRefPath::try_from_pathbuf(input_file).unwrap()).finish() {
        Ok(lf) => lf,
        Err(e) => panic!("Failed to read CSV file: {}", e),
    }
}

pub fn read_from_tsv(input_file: PathBuf) -> LazyFrame {
    match LazyCsvReader::new(PlRefPath::try_from_pathbuf(input_file).unwrap())
        .with_separator(b'\t')
        .finish()
    {
        Ok(lf) => lf,
        Err(e) => panic!("Failed to read TSV file: {}", e),
    }
}

pub fn read_from_file(input_file: PathBuf, separator: Option<u8>) -> LazyFrame {
    match separator {
        Some(sep) => match sep {
            b'\t' => read_from_tsv(input_file),
            b',' => read_from_csv(input_file),
            _ => panic!("Unsupported separator"),
        },
        _ => match get_extension_from_filename(&input_file) {
            Some(v) => match v {
                "tsv" => read_from_tsv(input_file),
                "xls" | "xlsx" => read_excel(&input_file, None, 0, None)
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

    let df = DataFrame::new_infer_height(series_vec);
    df
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::path::PathBuf;

    #[rstest]
    #[case("data.csv", Some("csv"))]
    #[case("archive.tar.gz", Some("gz"))]
    #[case("no_extension", None)]
    #[case(".bashrc", None)]
    #[case("/path/to/file.txt", Some("txt"))]
    fn test_get_extension_from_filename(#[case] filename: String, #[case] ext: Option<&str>) {
        assert_eq!(get_extension_from_filename(&PathBuf::from(filename)), ext);
    }

    #[rstest]
    #[case("\"quoted\"", "quoted")]
    #[case("'single_quoted'", "single_quoted")]
    #[case("not_quoted", "not_quoted")]
    #[case("'mismatched\"", "'mismatched\"")]
    #[case("", "")]
    #[case("''", "")]
    fn test_strip_quotes(#[case] s: &str, #[case] exp: String) {
        assert_eq!(strip_quotes(s), exp);
    }
}
