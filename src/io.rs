use anyhow::{Context as AnyhowContext, Result};
use calamine::{Reader, Xlsx, open_workbook};
use polars::prelude::*;
use std::io::{BufRead, BufReader, Cursor, Read};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::sync::Arc;

fn get_extension_from_filename(filename: &PathBuf) -> Option<&str> {
    Path::new(filename).extension().and_then(|s| s.to_str())
}

fn write_df_output(df: &mut DataFrame, output_file: &PathBuf, separator: u8) -> Result<()> {
    if output_file == Path::new("-") {
        let stdout = std::io::stdout();
        let mut output = stdout.lock();
        AnyhowContext::context(
            CsvWriter::new(&mut output)
                .with_separator(separator)
                .finish(df),
            "Failed to write to stdout",
        )?;
    } else {
        let mut output = AnyhowContext::context(
            std::fs::File::create(output_file),
            "Failed to create output file",
        )?;
        AnyhowContext::context(
            CsvWriter::new(&mut output)
                .with_separator(separator)
                .finish(df),
            "Failed to write to file",
        )?;
    }
    Ok(())
}

fn write_lazy_output(lf: LazyFrame, output_file: PathBuf, separator: u8) -> Result<()> {
    if output_file == Path::new("-") {
        let mut df = AnyhowContext::context(
            lf.collect(),
            "Failed to collect lazy frame prior to writing to stdout",
        )?;
        write_df_output(&mut df, &output_file, separator)?;
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

        let sink = AnyhowContext::context(
            lf.sink(
                SinkDestination::File {
                    target: SinkTarget::Path(PlRefPath::try_from_pathbuf(output_file).unwrap()),
                },
                FileWriteFormat::Csv(writer_opts),
                UnifiedSinkArgs::default(),
            ),
            "Failed to open CSV file for writing",
        )?;

        AnyhowContext::context(
            sink.collect(),
            "Failed to collect lazy frame prior to writing to file",
        )?;
    }
    Ok(())
}

pub trait WriteToCsvOrStdout {
    fn write_to_csv_or_stdout(self, output_file: PathBuf) -> Result<()>;
    fn write_to_tsv_or_stdout(self, output_file: PathBuf) -> Result<()>;
    fn write_to_flat_or_stdout(self, output_file: PathBuf, separator: Option<u8>) -> Result<()>
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
    fn write_to_csv_or_stdout(mut self, output_file: PathBuf) -> Result<()> {
        write_df_output(&mut self, &output_file, b',')
    }
    fn write_to_tsv_or_stdout(mut self, output_file: PathBuf) -> Result<()> {
        write_df_output(&mut self, &output_file, b'\t')
    }
}
impl WriteToCsvOrStdout for LazyFrame {
    fn write_to_csv_or_stdout(self, output_file: PathBuf) -> Result<()> {
        write_lazy_output(self, output_file, b',')
    }
    fn write_to_tsv_or_stdout(self, output_file: PathBuf) -> Result<()> {
        write_lazy_output(self, output_file, b'\t')
    }
}

pub fn read_from_csv(input_file: PathBuf) -> Result<LazyFrame> {
    AnyhowContext::with_context(
        LazyCsvReader::new(PlRefPath::try_from_pathbuf(input_file.clone())?).finish(),
        || format!("Failed to read CSV file: {}", input_file.to_string_lossy()),
    )
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
pub fn read_from_tsv(input_file: PathBuf) -> Result<LazyFrame> {
    AnyhowContext::with_context(
        LazyCsvReader::new(PlRefPath::try_from_pathbuf(input_file.clone())?)
            .with_separator(b'\t')
            .finish(),
        || format!("Failed to read TSV file: {}", input_file.to_string_lossy()),
    )
}

/// Reads data from a file and returns a `LazyFrame` based on the file's content and specified format.
///
/// # Parameters
/// - `input_file`: A `PathBuf` representing the path to the input file that is to be read.
/// - `separator`: An `Option<u8>` representing the optional field separator:
///   - `Some(sep)`: Read the file as delimited text using `sep`, skipping detection.
///   - `None`: Detect the format from the file's extension and content.
///
/// # Behavior
/// - If `separator` is provided, the file is read as delimited text with that separator.
/// - If `separator` is not provided:
///   - `"xls"`/`"xlsx"`: Reads the file as an Excel file (first sheet, no skipped rows).
///   - Otherwise the separator is guessed from the first lines of the file (see
///     [`guess_separator`]). If guessing is inconclusive, it falls back to tab for
///     `"tsv"` files or files without an extension, and comma for everything else.
///
/// # Returns
/// - A `LazyFrame` generated from the data in the file, ready for further manipulation.
///
/// # Errors
/// - Returns an error if the file cannot be opened or parsed.
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
pub fn read_from_file(input_file: PathBuf, separator: Option<u8>) -> Result<LazyFrame> {
    if input_file == Path::new("-") {
        return read_from_reader(std::io::stdin().lock(), separator)
            .context("Failed to read from stdin");
    }
    if let Some(sep) = separator {
        return read_delimited(input_file, sep);
    }
    match get_extension_from_filename(&input_file) {
        Some("xls" | "xlsx") => Ok(read_excel(&input_file, None, 0, None)?.lazy()),
        ext => {
            let fallback = if matches!(ext, None | Some("tsv")) {
                b'\t'
            } else {
                b','
            };
            let sep = guess_separator(&input_file)?.unwrap_or(fallback);
            read_delimited(input_file, sep)
        }
    }
}

const CANDIDATE_SEPARATORS: [u8; 4] = *b"\t,;|";
const SEP_GUESS_LINES: usize = 20;

/// Guesses the field separator of a delimited text file from its first lines.
///
/// Each candidate in [`CANDIDATE_SEPARATORS`] is counted on every sampled line,
/// ignoring occurrences inside double quotes. A candidate qualifies when it occurs
/// the same, non-zero number of times on every line; among those, the one yielding
/// the most fields wins.
///
/// # Returns
/// - `Ok(Some(sep))` when a separator was detected.
/// - `Ok(None)` when the file is empty or no candidate occurs consistently
///   (e.g. a single-column file).
///
/// # Errors
/// - Returns an error if the file cannot be opened or read.
pub fn guess_separator(path: &Path) -> Result<Option<u8>> {
    let file = AnyhowContext::with_context(std::fs::File::open(path), || {
        format!("Failed to open {}", path.display())
    })?;
    Ok(guess_separator_from_reader(BufReader::new(file))?)
}

/// Reads delimited text from any reader (e.g. stdin) into a `LazyFrame`.
///
/// The whole input is buffered in memory, since a stream can only be read once.
/// With `separator` set to `None`, the separator is guessed from the content
/// (see [`guess_separator`]) and defaults to tab when that is inconclusive.
///
/// # Errors
/// - Returns an error if the input cannot be read or parsed.
pub fn read_from_reader(mut reader: impl Read, separator: Option<u8>) -> Result<LazyFrame> {
    let mut bytes = Vec::new();
    reader.read_to_end(&mut bytes)?;
    let sep = match separator {
        Some(sep) => sep,
        None => guess_separator_from_reader(bytes.as_slice())?.unwrap_or(b'\t'),
    };
    let df = CsvReadOptions::default()
        .with_parse_options(CsvParseOptions::default().with_separator(sep))
        .into_reader_with_file_handle(Cursor::new(bytes))
        .finish()?;
    Ok(df.lazy())
}

/// Content-based separator detection shared by [`guess_separator`] and
/// [`read_from_reader`]. Only the first [`SEP_GUESS_LINES`] lines are consumed.
fn guess_separator_from_reader(reader: impl BufRead) -> std::io::Result<Option<u8>> {
    let lines: Vec<Vec<u8>> = reader
        .split(b'\n')
        .take(SEP_GUESS_LINES)
        .collect::<std::io::Result<_>>()?;
    let lines: Vec<&[u8]> = lines
        .iter()
        .map(|l| l.strip_suffix(b"\r").unwrap_or(l))
        .filter(|l| !l.is_empty())
        .collect();
    let Some(first) = lines.first() else {
        return Ok(None);
    };

    let count = |line: &[u8], sep: u8| {
        let mut in_quotes = false;
        line.iter()
            .filter(|&&b| {
                if b == b'"' {
                    in_quotes = !in_quotes;
                }
                !in_quotes && b == sep
            })
            .count()
    };

    Ok(CANDIDATE_SEPARATORS
        .iter()
        .copied()
        .filter_map(|sep| {
            let n = count(first, sep);
            (n > 0 && lines.iter().all(|l| count(l, sep) == n)).then_some((sep, n))
        })
        .max_by_key(|&(_, n)| n)
        .map(|(sep, _)| sep))
}

/// Reads a delimited text file with the given separator into a `LazyFrame`.
///
/// # Errors
/// - Returns an error if the file cannot be read or parsed.
pub fn read_delimited(input_file: PathBuf, separator: u8) -> Result<LazyFrame> {
    AnyhowContext::with_context(
        LazyCsvReader::new(PlRefPath::try_from_pathbuf(input_file.clone())?)
            .with_separator(separator)
            .finish(),
        || {
            format!(
                "Failed to read delimited file: {}",
                input_file.to_string_lossy()
            )
        },
    )
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

    PolarsContext::context(
        DataFrame::new_infer_height(series_vec),
        "Failed to read Excel file",
    )
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

    #[rstest]
    #[case::parse_csv(Some(b','), "a,b\n1,2\n", Some("csv"))]
    #[case::parse_tsv(Some(b'\t'), "a\tb\n1\t2\n", Some("tsv"))]
    #[case::parse_tsv_default_sep(None, "a\tb\n1\t2\n", Some("tsv"))]
    #[case::parse_tsv_no_ext(None, "a\tb\n1\t2\n", None)]
    #[case::explicit_semicolon(Some(b';'), "a;b\n1;2\n", Some("csv"))]
    #[case::guess_tabs_in_csv(None, "a\tb\n1\t2\n", Some("csv"))]
    #[case::guess_commas_in_tsv(None, "a,b\n1,2\n", Some("tsv"))]
    #[case::guess_semicolon_txt(None, "a;b\n1;2\n", Some("txt"))]
    #[case::guess_pipe_no_ext(None, "a|b\n1|2\n", None)]
    #[case::guess_crlf(None, "a\tb\r\n1\t2\r\n", Some("csv"))]
    #[case::guess_ignores_quoted(None, "a\tb\n\"x,y,z\"\t2\n", Some("csv"))]
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

        let df = collect_df(read_from_file(path, sep).unwrap());

        assert_eq!(df.shape(), (1, 2));
        assert_eq!(df.get_column_names()[0].as_str(), "a");
        assert_eq!(df.get_column_names()[1].as_str(), "b");
    }

    #[rstest]
    #[case::single_column(Some("csv"), "a\n1\n2\n", None)]
    #[case::empty(Some("csv"), "", None)]
    #[case::inconsistent_counts(Some("csv"), "a,b\n1,2,3\n", None)]
    #[case::prefers_most_fields(Some("txt"), "a;b,c;d\n1;2,3;4\n", Some(b';'))]
    fn guess_separator_cases(
        #[case] ext: Option<&str>,
        #[case] text: &str,
        #[case] expected: Option<u8>,
    ) {
        let temp = TempDir::default();
        let path = temp
            .with_file_name(get_timestamp())
            .with_extension(ext.unwrap_or(""));
        write_text_file(&path, text);

        assert_eq!(guess_separator(&path).unwrap(), expected);
    }

    #[rstest]
    #[case::guess_comma(None, "a,b\n1,2\n")]
    #[case::guess_semicolon(None, "a;b\n1;2\n")]
    #[case::fallback_tab(None, "a\tb\n1\t2\n")]
    #[case::explicit_sep(Some(b'|'), "a|b\n1|2\n")]
    fn read_from_reader_cases(#[case] sep: Option<u8>, #[case] text: &str) {
        let df = collect_df(read_from_reader(text.as_bytes(), sep).unwrap());

        assert_eq!(df.shape(), (1, 2));
        assert_eq!(df.get_column_names()[0].as_str(), "a");
        assert_eq!(df.get_column_names()[1].as_str(), "b");
    }

    #[rstest]
    fn read_from_reader_reads_past_guessed_lines() {
        let text: String = std::iter::once("a,b\n".to_string())
            .chain((0..100).map(|i| format!("{i},{i}\n")))
            .collect();

        let df = collect_df(read_from_reader(text.as_bytes(), None).unwrap());

        assert_eq!(df.shape(), (100, 2));
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

        test_df.write_to_flat_or_stdout(path.clone(), sep).unwrap();
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

        test_lf.write_to_flat_or_stdout(path.clone(), sep).unwrap();
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

        let df = collect_df(read_from_file(path.clone(), Some(b'\t')).unwrap());

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

    #[rstest]
    #[case::normal_file("tests/data/plate_reader_data_expected.xlsx")]
    fn test_read_excel_from_generic(#[case] path: PathBuf) {
        let df = read_from_file(path, None).unwrap().collect().unwrap();

        assert!(df.width() > 0, "expected at least one column");
        assert!(df.height() > 0, "expected at least one row");
    }
}
