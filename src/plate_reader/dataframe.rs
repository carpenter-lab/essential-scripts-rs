use crate::io;
use polars::datatypes::{AnyValue, PlSmallStr, UInt32Chunked};
use polars::error::{PolarsError, PolarsResult};
use polars::frame::DataFrame;
use polars::prelude::{IntoColumn, NamedFrom, Series};

const LETTERS: [char; 8] = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'];

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

/// Split a block of rows into 'stride' slices (every nth row starting with offset i),
/// convert each to vector-of-vectors for flattening later.
pub fn flatten_by_stride(
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

pub fn prepare_data(
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
        let sliced = col.slice(i64::try_from(data_row_start).unwrap(), data_stride * 8);
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
