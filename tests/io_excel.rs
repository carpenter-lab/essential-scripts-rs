use essential_scripts_rs::io::{read_excel, read_from_file};
use std::path::PathBuf;

fn fixture(path: &str) -> PathBuf {
    PathBuf::from(path)
}

#[test]
fn read_excel_default_sheet_succeeds() {
    let path = fixture("test/plate_reader_data_expected.xlsx");
    let df = read_from_file(path, None)
        .collect()
        .expect("expected default-sheet read to succeed");

    assert!(df.width() > 0, "expected at least one column");
    assert!(df.height() > 0, "expected at least one row");
}

#[test]
fn read_excel_header_false_generates_column_names() {
    let path = fixture("test/plate_reader_data.xlsx");
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
fn read_excel_fails_with_no_header() {
    let path = fixture("test/plate_reader_data.xlsx");

    let df = read_excel(&path, None, 0, None);
    assert!(
        df.is_err(),
        "expected reading to fail with no header in file"
    );
}

#[test]
fn read_excel_skiprows_reduces_height() {
    let path = fixture("test/plate_reader_data_expected.xlsx");

    let df_no_skip = read_excel(&path, None, 0, None).expect("read without skiprows failed");
    let df_skip_1 = read_excel(&path, None, 1, None).expect("read with skiprows=1 failed");

    assert!(
        df_skip_1.height() <= df_no_skip.height(),
        "skiprows should not increase row count"
    );
}

#[test]
fn read_excel_nonexistent_sheet_returns_error() {
    let path = fixture("test/plate_reader_data.xlsx");

    let err = read_excel(&path, Some("this_sheet_does_not_exist"), 0, None)
        .expect_err("expected error for nonexistent sheet");

    let msg = format!("{err}");
    assert!(
        msg.contains("Failed to get sheet"),
        "expected sheet lookup error, got: {msg}"
    );
}

#[test]
fn read_excel_nonexistent_file_returns_error() {
    let path = fixture("test/this_file_does_not_exist.xlsx");

    let err = read_excel(&path, None, 0, None).expect_err("expected file-open error");
    let msg = format!("{err}");

    assert!(
        msg.contains("Failed to open Excel file"),
        "expected open-file error, got: {msg}"
    );
}

/// Optional: enable when you add a fixture with duplicate header names.
/// e.g. `test/duplicate_headers.xlsx` first row: `A,A,B`
#[test]
#[ignore = "requires test/duplicate_headers.xlsx fixture"]
fn read_excel_duplicate_headers_returns_error() {
    let path = fixture("test/duplicate_headers.xlsx");

    let err =
        read_excel(&path, None, 0, None).expect_err("expected duplicate header validation error");
    let msg = format!("{err}");

    assert!(
        msg.contains("Duplicate column names found in Excel sheet"),
        "expected duplicate-header error, got: {msg}"
    );
}
