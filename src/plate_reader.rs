use crate::io;
use crate::io::WriteToCsvOrStdout;
use calamine::{Reader, Xlsx, open_workbook};
use clap::Subcommand;
use polars::prelude::*;
use std::path::{Path, PathBuf};

#[derive(Subcommand)]
#[command(about = "Reformat plate reader data into useful format")]
pub enum Commands {
    ReformatPlateReaderData {
        #[arg(required = true, help = "Input Excel file to process")]
        input_file: PathBuf,

        #[arg(
            required = true,
            help = "Output directory path. Will create one CSV per sheet."
        )]
        output_path: PathBuf,
    },
}

const LETTERS: [char; 8] = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'];

/// Split a block of rows into 'stride' slices (every nth row starting with offset i),
/// convert each to vector-of-vectors for flattening later.
fn flatten_by_stride(
    block: &DataFrame,
    stride: usize,
) -> PolarsResult<Vec<Vec<Vec<AnyValue<'_>>>>> {
    if stride == 0 {
        return Ok(Vec::new());
    }

    let columns = block.columns();
    let ncols = columns.len();
    let height = block.height();

    // Helper to extract a single row's values across all columns
    let get_row = |row_idx: usize| -> PolarsResult<Vec<AnyValue<'_>>> {
        let mut row_values = Vec::with_capacity(ncols);
        for col in columns {
            row_values.push(col.get(row_idx)?);
        }
        Ok(row_values)
    };

    // If stride > height, pad with empty vectors
    if stride > height {
        let mut group = Vec::with_capacity(stride);
        for row_idx in 0..height {
            group.push(get_row(row_idx)?);
        }
        for _ in height..stride {
            group.push(Vec::new());
        }
        return Ok(vec![group]);
    }

    // Split rows by stride: group[i] contains rows at positions i, i+stride, i+2*stride, ...
    let mut result = Vec::with_capacity(stride);
    for offset in 0..stride {
        let mut group = Vec::new();
        let mut row_idx = offset;
        while row_idx < height {
            group.push(get_row(row_idx)?);
            row_idx += stride;
        }
        result.push(group);
    }

    Ok(result)
}

fn generate_well_ids(rows: &[char], cols: &[u32]) -> Vec<String> {
    rows.iter()
        .flat_map(|r| cols.iter().map(move |c| format!("{r}{c}")))
        .collect()
}

/// Generate well IDs for a standard 96-well plate (A1-H12)
fn generate_default_well_ids() -> Vec<String> {
    let cols: Vec<u32> = (1..=12).collect();
    generate_well_ids(&LETTERS, &cols)
}

/// Drop rows that have fewer than `threshold` non-null values
fn drop_mostly_empty_rows(df: &DataFrame, threshold: usize) -> PolarsResult<DataFrame> {
    let height = df.height();
    let mut keep_rows = Vec::with_capacity(height);

    for row_idx in 0..height {
        let mut non_null_count = 0;

        for col in df.columns() {
            let val = col.get(row_idx)?;

            if !val.is_null() && !val.is_nan() && val != AnyValue::from("") {
                non_null_count += 1;
            }
        }

        if non_null_count >= threshold {
            keep_rows.push(u32::try_from(row_idx).map_err(|e| {
                PolarsError::ComputeError(
                    format!("Row index {row_idx} exceeds u32 max: {e}").into(),
                )
            })?);
        }
    }

    // Create a mask and filter
    let indices = UInt32Chunked::from_vec("idx".into(), keep_rows);
    df.take(&indices)
}

/// Clean '[Concentration]' style string column:
///   - values containing '<' become 0.0
///   - numeric parts are parsed as float
///   - values containing '>' become NaN (unknown upper bound)
fn clean_concentration(col: &Column) -> PolarsResult<Series> {
    // Convert to string type if not already
    let str_series = col.cast(&DataType::String)?;
    let str_chunked = str_series.str()?;

    let len = str_chunked.len();
    let mut result = Vec::with_capacity(len);

    for opt_str in str_chunked {
        match opt_str {
            None => result.push(None),
            Some(s) => {
                let contains_lt = s.contains('<');
                let contains_gt = s.contains('>');

                if contains_gt {
                    // Values with '>' become NaN
                    result.push(None);
                } else if contains_lt {
                    // Values with '<' become 0.0
                    result.push(Some(0.0));
                } else {
                    // Extract numeric part using regex-like pattern
                    // Match pattern: optional digits, optional decimal point, digits
                    let numeric_part = extract_numeric(s);
                    result.push(numeric_part);
                }
            }
        }
    }

    Ok(Series::new(col.name().clone(), result))
}

/// Extract the first numeric value from a string
/// Matches pattern: [0-9]*\.?[0-9]+
fn extract_numeric(s: &str) -> Option<f64> {
    let mut num_str = String::new();
    let mut has_digit = false;
    let mut has_decimal = false;

    for ch in s.chars() {
        if ch.is_ascii_digit() {
            num_str.push(ch);
            has_digit = true;
        } else if ch == '.' && !has_decimal && has_digit {
            num_str.push(ch);
            has_decimal = true;
        } else if ch == '.' && !has_decimal {
            // Decimal before any digit
            num_str.push(ch);
            has_decimal = true;
        } else if has_digit {
            // Stop at first non-numeric after we've found digits
            break;
        }
    }

    if has_digit {
        num_str.parse::<f64>().ok()
    } else {
        None
    }
}

fn prepare_data(
    df: &DataFrame,
    value_col_start: usize,
    value_col_end: Option<usize>,
    setup_stride: usize,
    data_stride: usize,
) -> PolarsResult<(DataFrame, Vec<String>, Vec<String>)> {
    let n_cols = df.width();

    // Determine the range of columns with plate values
    let col_end = value_col_end.unwrap_or(n_cols - 1);

    // Extract setup slice (first setup_stride * 8 rows)
    let setup_end_row = setup_stride * 8;
    let col_names: Vec<String> = df.get_column_names()[value_col_start..col_end]
        .iter()
        .map(|s| io::strip_quotes(s.as_str()))
        .collect::<Vec<String>>();
    let col_names: Vec<&str> = col_names.iter().map(String::as_str).collect();
    // Build setup slice by extracting columns and slicing rows
    let mut setup_cols = Vec::new();
    for &col_name in &col_names {
        let col = df.column(col_name)?;
        let sliced = col.slice(0, setup_end_row);
        setup_cols.push(sliced);
    }
    let setup_slice = DataFrame::new_infer_height(setup_cols)?;

    // Extract data slice (after setup + 5 rows gap)
    let data_row_start = setup_end_row + 5;
    let mut data_cols = Vec::new();
    for &col_name in &col_names {
        let col = df.column(col_name)?;
        let sliced = col.slice(data_row_start as i64, data_stride * 8);
        data_cols.push(sliced);
    }
    let data_slice = DataFrame::new_infer_height(data_cols)?;

    // Extract labels from the last column
    let last_col_name = df.get_column_names()[n_cols - 1];
    let last_col = df.column(last_col_name)?;

    // Get setup names (first setup_stride rows from last column)
    let mut setup_names = Vec::with_capacity(setup_stride);
    for i in 0..setup_stride {
        let val = last_col.get(i)?;
        setup_names.push(val.to_string());
    }

    // Get data names (data_stride rows starting from data_row_start)
    let mut data_names = Vec::with_capacity(data_stride);
    for i in 0..data_stride {
        let val = last_col.get(data_row_start + i)?;
        data_names.push(val.to_string());
    }

    // Reshape by strides
    let setup_slices = flatten_by_stride(&setup_slice, setup_stride)?;
    let data_slices = flatten_by_stride(&data_slice, data_stride)?;

    // Flatten each slice into a single vector (row-major order)
    let mut all_columns: Vec<(String, Vec<AnyValue>)> = Vec::new();

    // Process setup slices
    for (i, slice) in setup_slices.iter().enumerate() {
        let flattened: Vec<AnyValue> = slice.iter().flat_map(|row| row.iter().cloned()).collect();
        all_columns.push((setup_names[i].clone(), flattened));
    }

    // Process data slices
    for (i, slice) in data_slices.iter().enumerate() {
        let flattened: Vec<AnyValue> = slice.iter().flat_map(|row| row.iter().cloned()).collect();
        all_columns.push((data_names[i].clone(), flattened));
    }

    // Add well_id column
    let well_ids = generate_default_well_ids();

    // Build DataFrame from columns
    let mut series_vec = Vec::new();

    for (name, values) in all_columns {
        // Convert AnyValue vector to Series
        let s = Series::from_any_values(PlSmallStr::from_str(&name), &values, true)?;
        series_vec.push(s.into_column());
    }

    // Add well_id as a series
    let well_id_series = Series::new("well_id".into(), well_ids);
    series_vec.push(well_id_series.into_column());

    let result_df = DataFrame::new_infer_height(series_vec)?;

    // Drop rows with fewer than 2 non-null values (equivalent to dropna(thresh=2))
    let result_df = drop_mostly_empty_rows(&result_df, 2)?;
    Ok((result_df, setup_names, data_names))
}

fn run(
    input_path: &PathBuf,
    concentration_column: &str,
    threshold_source: Option<&str>,
    threshold_value: Option<f64>,
    sheet_name: Option<&str>,
    skiprows: usize,
) -> PolarsResult<DataFrame> {
    let df = io::read_excel(input_path, sheet_name, skiprows, Some(false))?;

    // Prepare data with default strides (2 for setup, 3 for data)
    let (mut res, setup_names, data_names) = prepare_data(&df, 2, None, 2, 3)?;

    // Build adjusted/clean concentration columns if source exists
    let col_names = res.get_column_names();

    let has_concentration: String = col_names
        .iter()
        .filter(|col_name| io::strip_quotes(col_name.as_str()) == concentration_column)
        .map(|&name| name.as_str())
        .collect();

    let has_concentration = if has_concentration.is_empty() {
        None
    } else {
        Some(has_concentration)
    };
    //.any(|&name| strip_quotes(name.as_str()) == concentration_column);
    // Build column order: well_id, setup_names, data_names, optional adjusted/clean
    let mut column_order = vec!["well_id".to_string()];
    column_order.extend(setup_names.clone());
    column_order.extend(data_names.clone());
    if let Some(conc_col_use) = has_concentration {
        let adjusted_col_name = format!("{concentration_column} Adjusted");
        let clean_col_name = format!("{concentration_column} Clean");

        // Create adjusted column
        let concentration_series = res.column(&conc_col_use)?;

        let adjusted_series = clean_concentration(concentration_series)?;

        res = res
            .with_column(
                adjusted_series
                    .with_name(PlSmallStr::from_str(&adjusted_col_name))
                    .into_column(),
            )?
            .to_owned();

        // Determine which column provides the threshold
        let col_for_threshold = threshold_source.or_else(|| {
            if data_names.is_empty() {
                None
            } else {
                Some(data_names[0].as_str())
            }
        });

        // Create clean column
        let clean_series = if let (Some(threshold_col), Some(threshold_val)) =
            (col_for_threshold, threshold_value)
        {
            if res
                .get_column_names()
                .contains(&&PlSmallStr::from_str(threshold_col))
            {
                // Mask adjusted values where threshold column > threshold_value
                let threshold_series = res.column(threshold_col)?;
                let adjusted_col = res.column(&adjusted_col_name)?;

                // Create mask: true where threshold > threshold_value
                let mask = threshold_series.gt(&Column::new(
                    PlSmallStr::from_str("threshold"),
                    Series::new("threshold".into(), [threshold_val]),
                ))?;

                // Apply mask: set to null where mask is true
                let clean = adjusted_col.zip_with(
                    &mask,
                    &Series::new_null("null".into(), adjusted_col.len()).into_column(),
                )?;
                clean.with_name(PlSmallStr::from_str(&clean_col_name))
            } else {
                // Threshold column not found, use adjusted as clean
                res.column(&adjusted_col_name)?
                    .clone()
                    .with_name(PlSmallStr::from_str(&clean_col_name))
            }
        } else {
            // No threshold, use adjusted as clean
            res.column(&adjusted_col_name)?
                .clone()
                .with_name(PlSmallStr::from_str(&clean_col_name))
        };

        res = res.with_column(clean_series.into_column())?.to_owned();

        column_order.push(format!("{concentration_column} Adjusted"));
        column_order.push(format!("{concentration_column} Clean"));
    }

    // Filter to only columns that exist in the result
    let existing_columns: Vec<PlSmallStr> = column_order
        .iter()
        .filter(|c| res.get_column_names().contains(&&PlSmallStr::from_str(c)))
        .map(|c| PlSmallStr::from_str(c))
        .collect();

    // Sort by data columns if present
    let sort_cols: Vec<&str> = data_names
        .iter()
        .filter(|c| res.get_column_names().contains(&&PlSmallStr::from_str(c)))
        .map(String::as_str)
        .collect();

    let mut out = if sort_cols.is_empty() {
        res.select(existing_columns)?
    } else {
        res.sort(sort_cols, SortMultipleOptions::default())?
            .select(existing_columns)?
    };

    let binding = out.clone();
    let new_names: Vec<String> = binding
        .columns()
        .iter()
        .map(|c| io::strip_quotes(c.name().as_str()))
        .collect();

    out.set_column_names(&new_names.clone())?;
    // Select only the columns in the desired order
    out.select(&new_names)
}

fn reformat_plate_reader_data(input_file: &PathBuf, output_dir: &Path) -> PolarsResult<()> {
    let workbook: Xlsx<_> = open_workbook(input_file)
        .map_err(|e| PolarsError::ComputeError(format!("Failed to open Excel file: {e}").into()))?;
    let sheets = workbook.sheet_names().clone();
    for sheet in sheets {
        let df = run(input_file, "[Concentration]", None, None, Some(&sheet), 23)?;
        let mut output_file = output_dir.to_path_buf();
        let sheet_sanitized = sheet.replace([' ', '/'], "_");
        output_file.push(format!("{sheet_sanitized}.csv"));
        df.write_to_csv_or_stdout(output_file);
    }

    Ok(())
}

pub fn handle_command(cmd: Commands) {
    match cmd {
        Commands::ReformatPlateReaderData {
            input_file,
            output_path,
        } => {
            reformat_plate_reader_data(&input_file, &output_path).unwrap();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::read_excel;
    use polars_testing::asserts::{DataFrameEqualOptions, assert_dataframe_equal};
    use pretty_assertions::assert_eq;

    #[test]
    fn test_stride_of_one() {
        let df = df![
            "A" => [1, 3, 5],
            "B" => [2, 4, 6],
        ]
        .unwrap();
        let stride = 1;
        let result = flatten_by_stride(&df, stride).unwrap();

        let expected = vec![vec![
            vec![AnyValue::Int32(1), AnyValue::Int32(2)],
            vec![AnyValue::Int32(3), AnyValue::Int32(4)],
            vec![AnyValue::Int32(5), AnyValue::Int32(6)],
        ]];
        assert_eq!(expected, result);
    }

    #[test]
    fn test_stride_greater_than_one() {
        let df = df![
            "A" => [1, 2, 5, 7],
            "B" => [2, 4, 6, 8],
        ]
        .unwrap();
        let stride = 2;
        let result = flatten_by_stride(&df, stride).unwrap();

        let expected = vec![
            vec![
                vec![AnyValue::Int32(1), AnyValue::Int32(2)],
                vec![AnyValue::Int32(5), AnyValue::Int32(6)],
            ],
            vec![
                vec![AnyValue::Int32(2), AnyValue::Int32(4)],
                vec![AnyValue::Int32(7), AnyValue::Int32(8)],
            ],
        ];
        assert_eq!(expected, result);
    }

    #[test]
    fn test_stride_exceed_rows() {
        let df = df![
            "A" => [1, 3],
            "B" => [2, 4]
        ]
        .unwrap();
        let stride = 5;
        let result = flatten_by_stride(&df, stride).unwrap();

        let expected = vec![vec![
            vec![AnyValue::Int32(1), AnyValue::Int32(2)],
            vec![AnyValue::Int32(3), AnyValue::Int32(4)],
            vec![],
            vec![],
            vec![],
        ]];
        assert_eq!(expected, result);
    }

    #[test]
    fn test_run() {
        let input_path = PathBuf::from("test/plate_reader_data.xlsx");
        let df = run(
            &input_path,
            "[Concentration]",
            None,
            None,
            Some("1 to 2"),
            23,
        );
        let expected = read_excel(
            &PathBuf::from("test/plate_reader_data_expected.xlsx"),
            Some("1 to 2"),
            0,
            Some(true),
        );

        assert!(df.is_ok());
        assert!(expected.is_ok());
        let df = df.unwrap().sort(["well_id"], Default::default()).unwrap();
        let expected = expected
            .unwrap()
            .lazy()
            .with_column(col("[Concentration] Adjusted").cast(DataType::Float64))
            .with_column(col("[Concentration] Clean").cast(DataType::Float64))
            .collect()
            .unwrap()
            .sort(["well_id"], Default::default())
            .unwrap();

        assert_eq!(df.get_column_names(), expected.get_column_names());

        assert_dataframe_equal(
            &df,
            &expected.sort(["well_id"], Default::default()).unwrap(),
            DataFrameEqualOptions::new(),
        )
        .unwrap();
    }
}
