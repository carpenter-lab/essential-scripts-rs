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
            include_header: true,
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
                null: PlSmallStr::default(),
                line_terminator: PlSmallStr::from_str("\n"),
                quote_style: QuoteStyle::default(),
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
    fn write_to_csv_or_stdout(self, output_file: PathBuf);
    fn write_to_tsv_or_stdout(self, output_file: PathBuf);
    fn write_to_flat_or_stdout(self, output_file: PathBuf, separator: Option<u8>)
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
    fn write_to_csv_or_stdout(mut self, output_file: PathBuf) {
        write_df_output(&mut self, &output_file, b',');
    }
    fn write_to_tsv_or_stdout(mut self, output_file: PathBuf) {
        write_df_output(&mut self, &output_file, b'\t');
    }
}
impl WriteToCsvOrStdout for LazyFrame {
    fn write_to_csv_or_stdout(self, output_file: PathBuf) {
        write_lazy_output(self, output_file, b',');
    }
    fn write_to_tsv_or_stdout(self, output_file: PathBuf) {
        write_lazy_output(self, output_file, b'\t');
    }
}

pub fn read_from_csv(input_file: PathBuf) -> LazyFrame {
    match LazyCsvReader::new(PlRefPath::try_from_pathbuf(input_file).unwrap()).finish() {
        Ok(lf) => lf,
        Err(e) => panic!("Failed to read CSV file: {e}"),
    }
}

/// Reads a TSV (Tab-Separated Values) file into a `LazyFrame`.
///
/// This function takes the path to a TSV file as input and uses a `LazyCsvReader`
/// to parse the file. The separator is explicitly set to a tab character (`\t`).
///
/// # Arguments
///
/// * `input_file` - A `PathBuf` representing the path to the TSV file to be read.
///
/// # Returns
///
/// * A `LazyFrame` containing the data from the TSV file.
///
/// # Panics
///
/// This function will panic if:
/// - The provided path cannot be converted to a `PlRefPath`.
/// - The TSV file cannot be read or parsed for any reason.
///
/// # Example
///
/// ```rust
/// use std::path::PathBuf;
/// use essential_scripts_rs::io::read_from_tsv;
///
/// let input_file = PathBuf::from("data/example.tsv");
/// let lazy_frame = read_from_tsv(input_file);
/// // Use `lazy_frame` for further processing...
/// ```
///
/// # Dependencies
///
/// This function relies on the `LazyCsvReader` to load the file, and expects the
/// file to use tab characters as the delimiter.
///
/// # Notes
///
/// Ensure that the input file exists and is accessible, and the format strictly
/// follows the TSV structure (tab-separated values).
pub fn read_from_tsv(input_file: PathBuf) -> LazyFrame {
    match LazyCsvReader::new(PlRefPath::try_from_pathbuf(input_file).unwrap())
        .with_separator(b'\t')
        .finish()
    {
        Ok(lf) => lf,
        Err(e) => panic!("Failed to read TSV file: {e}"),
    }
}

/// Reads data from a file and returns a `LazyFrame` based on the file's content and specified format.
///
/// # Parameters
/// - `input_file`: A `PathBuf` representing the path to the input file that is to be read.
/// - `separator`: An `Option<u8>` representing the optional byte separator to specify the file format:
///   - `Some(b'\t')`: Treat the file as a tab-separated values (TSV) file.
///   - `Some(b',')`: Treat the file as a comma-separated values (CSV) file.
///   - Any other `Some` value will result in a panic indicating an unsupported separator.
///   - `None`: Automatically detect the file format based on its extension.
///
/// # Behavior
/// - If `separator` is provided:
///   - Handles TSV files for a tab (`\t`) separator.
///   - Handles CSV files for a comma (`,`) separator.
///   - Panics if an unsupported separator is provided.
/// - If `separator` is not provided:
///   - Determines the file format based on the file's extension:
///     - `"tsv"`: Reads the file as a TSV file.
///     - `"xls"`/`"xlsx"`: Reads the file as an Excel file (**Note**: The `read_excel` function expects the following parameters: no sheet name, sheet index `0`, and no additional options).
///     - Any other or unknown extension: Reads the file as a CSV file.
///   - Defaults to reading the file as TSV if no extension is found.
///
/// # Returns
/// - A `LazyFrame` generated from the data in the file, ready for further manipulation.
///
/// # Panics
/// - Panics if the separator provided is unsupported.
/// - Panics if the Excel file reading fails.
///
/// # Examples
/// ```no_run
/// use std::path::PathBuf;
/// use essential_scripts_rs::io::read_from_file;
///
/// // Example: Reading a file with a specified separator
/// let input_path = PathBuf::from("data.csv");
/// let separator = Some(b',');
/// let lazy_frame = read_from_file(input_path, separator);
///
/// // Example: Auto-detecting file format using the file extension
/// let input_path = PathBuf::from("data.xlsx");
/// let lazy_frame = read_from_file(input_path, None);
/// ```
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

#[must_use]
pub fn strip_quotes(s: &str) -> String {
    if (s.starts_with('"') && s.ends_with('"')) || (s.starts_with('\'') && s.ends_with('\'')) {
        s[1..s.len() - 1].to_string()
    } else {
        s.to_string()
    }
}

/// Reads an Excel file and parses its content into a Polars `DataFrame`.
///
/// This function leverages the `calamine` crate to read an Excel file
/// and provides the ability to specify the sheet to load, rows to skip,
/// and whether to treat the first row as a header row.
///
/// # Arguments
///
/// - `input_path` - A `PathBuf` representing the path to the Excel file to be read.
/// - `sheet_name` - An optional `&str` specifying the name of the sheet to read. If `None`, the first sheet is used.
/// - `skiprows` - A `usize` specifying the number of rows to skip when parsing the sheet's data.
/// - `header` - An optional `bool` indicating whether to consider the first row as a header.
///   - If `Some(true)`, the first row is used as column names.
///   - If `Some(false)`, standard column names (`"column_0", "column_1"`, etc.) are used.
///   - If `None`, the behavior defaults to considering the first row as a header.
///
/// # Returns
///
/// - `PolarsResult<DataFrame>` - A `PolarsResult` containing the resulting `DataFrame` if successful,
///   or an error message if any errors occur during parsing.
///
/// # Errors
///
/// This function can return the following errors:
/// - `PolarsError::ComputeError` if:
///   - The file cannot be opened.
///   - A specified sheet cannot be found or accessed.
///   - Duplicate column names are found (when using headers).
///   - Other issues arise while reading or converting the Excel data.
///
/// # Examples
///
/// ```rust
/// # use polars::prelude::*;
/// # use std::path::PathBuf;
/// # use essential_scripts_rs::io::read_excel;
///
/// let path = PathBuf::from("example.xlsx");
/// let sheet_name = Some("Sheet1");
/// let skiprows = 0;
/// let header = Some(true);
///
/// match read_excel(&path, sheet_name, skiprows, header) {
///     Ok(df) => println!("{:?}", df),
///     Err(e) => eprintln!("Failed to read Excel file: {}", e),
/// }
/// ```
///
/// # Notes
///
/// - Duplicate column names are not allowed when treating the first row as a header,
///   and will result in an error.
/// - Each cell in the sheet is initially parsed as a `String`.
///
pub fn read_excel(
    input_path: &PathBuf,
    sheet_name: Option<&str>,
    skiprows: usize,
    header: Option<bool>,
) -> PolarsResult<DataFrame> {
    // Read Excel file using calamine
    let mut workbook: Xlsx<_> = open_workbook(input_path)
        .map_err(|e| PolarsError::ComputeError(format!("Failed to open Excel file: {e}").into()))?;

    let sheet_names = workbook.sheet_names();
    let sheet_name = match sheet_name {
        Some(name) => name,
        None => sheet_names
            .first()
            .ok_or_else(|| PolarsError::ComputeError("Excel file contains no sheets".into()))?
            .as_str(),
    };

    // Get the specified sheet
    let range = workbook.worksheet_range(sheet_name).map_err(|e| {
        PolarsError::ComputeError(format!("Failed to get sheet '{sheet_name}': {e}").into())
    })?;

    // Convert the range to a DataFrame
    let (height, width) = range.get_size();
    let col_names = match range.headers() {
        Some(headers) => match header {
            Some(true) | None => headers
                .iter()
                .map(|s| strip_quotes(&s.clone()))
                .collect::<Vec<String>>(),
            Some(false) => (0..width)
                .map(|i| format!("column_{i}"))
                .collect::<Vec<String>>(),
        },
        None => (0..width)
            .map(|i| format!("column_{i}"))
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
        if let Some(_h) = range.headers()
            && row_idx == skiprows
        {
            continue; // Skip header row
        }
        for (col_idx, item) in columns.iter_mut().enumerate().take(width) {
            let cell_value = range.get((row_idx, col_idx));
            let value = cell_value.map(ToString::to_string);
            item.push(value);
        }
    }

    // Create Series from columns
    let mut series_vec: Vec<Column> = Vec::new();
    for (col_idx, col_data) in columns.into_iter().enumerate() {
        let col_name = &col_names[col_idx];
        let s = Series::new(PlSmallStr::from_str(col_name), col_data);
        series_vec.push(s.into_column());
    }

    DataFrame::new_infer_height(series_vec)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_helpers::*;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use rstest_reuse::{apply, template};
    use std::fs;
    use std::path::PathBuf;
    use temp_testdir::TempDir;

    fn collect_df(lf: LazyFrame) -> DataFrame {
        lf.collect().expect("failed to collect lazyframe")
    }

    #[fixture]
    fn test_df() -> DataFrame {
        df![
            "a" => &[1],
            "b" => &[2]
        ]
        .expect("failed to build dataframe")
    }

    #[fixture]
    fn test_lf(test_df: DataFrame) -> LazyFrame {
        test_df.lazy()
    }

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

    #[template]
    #[rstest]
    #[case::parse_csv(Some(b','), "a,b\n1,2\n", Some("csv"))]
    #[case::parse_tsv(Some(b'\t'), "a\tb\n1\t2\n", Some("tsv"))]
    #[case::parse_tsv_default_sep(None, "a\tb\n1\t2\n", Some("tsv"))]
    #[should_panic(expected = "Unsupported separator")]
    #[case::fail_unsupported_sep(Some(b';'), "a\tb\n1\t2\n", Some("csv"))]
    #[case::parse_tsv_no_ext(None, "a\tb\n1\t2\n", None)]
    fn file_case_template(#[case] sep: Option<u8>, #[case] text: &str, #[case] ext: Option<&str>) {}

    #[apply(file_case_template)]
    fn read_from_file_with_sep(
        #[case] sep: Option<u8>,
        #[case] text: &str,
        #[case] ext: Option<&str>,
    ) {
        let temp = TempDir::default();
        let path = temp
            .with_file_name(get_timestamp())
            .with_extension(ext.unwrap_or(""));
        write_text_file(&path, text);

        let df = collect_df(read_from_file(path, sep));

        assert_eq!(df.shape(), (1, 2));
        assert_eq!(df.get_column_names()[0].as_str(), "a");
        assert_eq!(df.get_column_names()[1].as_str(), "b");
    }

    #[apply(file_case_template)]
    fn write_to_flat(
        #[case] sep: Option<u8>,
        #[case] text: &str,
        #[case] ext: Option<&str>,
        test_df: DataFrame,
    ) {
        let temp = TempDir::default();
        let path = temp
            .with_file_name(get_timestamp())
            .with_extension(ext.unwrap_or(""));

        test_df.write_to_flat_or_stdout(path.clone(), sep);
        let written = fs::read_to_string(&path).expect("failed to read output file");
        assert!(
            written.contains(text),
            "expected row {text}, got: {written}"
        );
    }

    #[apply(file_case_template)]
    fn write_to_flat_lazy(
        #[case] sep: Option<u8>,
        #[case] text: &str,
        #[case] ext: Option<&str>,
        test_lf: LazyFrame,
    ) {
        let temp = TempDir::default();
        let path = temp
            .with_file_name(get_timestamp())
            .with_extension(ext.unwrap_or(""));

        test_lf.write_to_flat_or_stdout(path.clone(), sep);
        let written = fs::read_to_string(&path).expect("failed to read output file");
        assert!(
            written.contains(text),
            "expected row {text}, got: {written}"
        );
    }

    #[rstest]
    fn tsv_parsing_regression_tab_delimited_file_produces_two_columns() {
        let temp = TempDir::default();
        let path = temp.with_file_name(get_timestamp()).with_extension("tsv");
        write_text_file(&path, "col1\tcol2\nv1\tv2\n");

        let df = collect_df(read_from_file(path.clone(), Some(b'\t')));

        // This test is intentionally strict: if TSV is parsed as CSV, you'll get one column.
        assert_eq!(
            df.width(),
            2,
            "TSV should parse into 2 columns; got {} columns with names {:?}",
            df.width(),
            df.get_column_names()
        );
    }

    #[test]
    fn read_excel_header_false_generates_column_names() {
        let path = PathBuf::from("tests/data/plate_reader_data.xlsx");
        let df = read_excel(&path, None, 0, Some(false))
            .expect("expected read with header=false to succeed");

        let names = df.get_column_names();
        assert!(
            !names.is_empty(),
            "expected at least one generated column name"
        );
        assert_eq!(names[0].as_str(), "column_0");
    }

    #[test]
    fn read_excel_skiprows_reduces_height() {
        let path = PathBuf::from("tests/data/plate_reader_data_expected.xlsx");

        let df_no_skip = read_excel(&path, None, 0, None).expect("read without skiprows failed");
        let df_skip_1 = read_excel(&path, None, 1, None).expect("read with skiprows=1 failed");

        assert!(
            df_skip_1.height() < df_no_skip.height(),
            "skiprows should not increase row count"
        );
    }

    #[rstest]
    #[should_panic(
        expected = "Failed to get sheet 'this_sheet_does_not_exist': Worksheet 'this_sheet_does_not_exist' not found"
    )]
    #[case::nonexisting_sheet(
        "tests/data/plate_reader_data.xlsx",
        Some("this_sheet_does_not_exist"),
        None
    )]
    #[should_panic(
        expected = "Failed to open Excel file: I/O error: No such file or directory (os error 2)"
    )]
    #[case::nonexisting_file("tests/data/this_file_does_not_exist.xlsx", None, None)]
    #[should_panic(expected = "Duplicate column names found in Excel sheet")]
    #[case::duplicate_header_none("tests/data/duplicate_headers.xlsx", None, None)]
    #[should_panic(expected = "Duplicate column names found in Excel sheet")]
    #[case::duplicate_header_true("tests/data/duplicate_headers.xlsx", None, Some(true))]
    #[case::duplicate_header_false("tests/data/duplicate_headers.xlsx", None, Some(false))]
    #[should_panic(expected = "Duplicate column names found in Excel sheet")]
    #[case::normal_file("tests/data/plate_reader_data.xlsx", None, None)]
    fn test_read_excel(
        #[case] path: PathBuf,
        #[case] sheet: Option<&str>,
        #[case] header: Option<bool>,
    ) {
        let df = read_excel(&path, sheet, 0, header).unwrap();

        assert!(df.width() > 0, "expected at least one column");
        assert!(df.height() > 0, "expected at least one row");
    }
}
