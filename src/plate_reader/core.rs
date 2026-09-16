use crate::io;
use crate::io::WriteToCsvOrStdout;
use crate::plate_reader::dataframe;
use anyhow::{Context, Result};
use calamine::{Reader, Xlsx, open_workbook};
use polars::datatypes::{DataType, PlSmallStr};
use polars::frame::DataFrame;
use polars::prelude::*;
use std::path::{Path, PathBuf};

/// Clean '[Concentration]' style string column:
///   - values containing '<' become 0.0
///   - numeric parts are parsed as float
///   - values containing '>' become NaN (unknown upper bound)
fn clean_concentration(col: &Column) -> Result<Series> {
    // Convert to string type if not already
    let str_series = col.cast(&DataType::String)?;
    let str_chunked = str_series.str()?;

    let len = str_chunked.len();
    let mut result = Vec::with_capacity(len);

    for opt_str in str_chunked.iter() {
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
        } else if ch == '.' && !has_decimal {
            num_str.push(ch);
            has_decimal = true;
        } else if ch == '.' && has_decimal {
            // parse two decimal points as invalid
            has_digit = false;
            break;
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

fn run(
    input_path: &PathBuf,
    concentration_column: &str,
    threshold_source: Option<&str>,
    threshold_value: Option<f64>,
    sheet_name: Option<&str>,
    skiprows: usize,
) -> Result<DataFrame> {
    let df = io::read_excel(input_path, sheet_name, skiprows, Some(false))?;

    // Prepare data with default strides (2 for setup, 3 for data)
    let (mut res, setup_names, data_names) = Context::context(
        dataframe::prepare_data(&df, 2, None, 2, 3),
        "Could not parse nested tables",
    )?;

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
    Context::context(out.select(&new_names), "Could not select columns")
}

pub fn reformat_plate_reader_data(input_file: &PathBuf, output_dir: &Path) -> Result<()> {
    let workbook: Xlsx<_> = Context::with_context(open_workbook(input_file), || {
        format!("Failed to open Excel file: {}", input_file.display())
    })?;
    let sheets = workbook.sheet_names().clone();
    for sheet in sheets {
        let df = run(input_file, "[Concentration]", None, None, Some(&sheet), 23)?;
        let mut output_file = output_dir.to_path_buf();
        let sheet_sanitized = sheet.replace([' ', '/'], "_");
        output_file.push(format!("{sheet_sanitized}.csv"));
        df.write_to_csv_or_stdout(output_file)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::read_excel;
    use crate::plate_reader::dataframe::flatten_by_stride;
    use polars_testing::asserts::{DataFrameEqualOptions, assert_dataframe_equal};
    use pretty_assertions::assert_eq;
    use rstest::rstest;

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
    fn test_stride_equal_to_rows() {
        let df = df![
            "A" => [1, 2, 3],
            "B" => [4, 5, 6],
        ]
        .unwrap();

        // stride == height
        let result = flatten_by_stride(&df, 3).unwrap();

        let expected = vec![
            vec![vec![AnyValue::Int32(1), AnyValue::Int32(4)]],
            vec![vec![AnyValue::Int32(2), AnyValue::Int32(5)]],
            vec![vec![AnyValue::Int32(3), AnyValue::Int32(6)]],
        ];

        assert_eq!(result, expected);
    }

    #[test]
    fn test_run() {
        let input_path = PathBuf::from("tests/data/plate_reader_data.xlsx");
        let df = run(
            &input_path,
            "[Concentration]",
            None,
            None,
            Some("1 to 2"),
            23,
        );
        let expected = read_excel(
            &PathBuf::from("tests/data/plate_reader_data_expected.xlsx"),
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

    #[rstest]
    #[case("123.45", Some(123.45))]
    #[case("0", Some(0.0))]
    #[case("1", Some(1.0))]
    #[case("abc", None)]
    #[case("123abc", Some(123.0))]
    #[case("123.45abc", Some(123.45))]
    #[case("123.45.67", None)]
    #[case("1..2", None)]
    #[case(">0.1", Some(0.1))]
    fn test_extract_numeric(#[case] input: &str, #[case] expected: Option<f64>) {
        let result = extract_numeric(input);
        assert_eq!(result, expected);
    }
}
