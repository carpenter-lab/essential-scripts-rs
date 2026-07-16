use essential_scripts_rs::io::{WriteToCsvOrStdout, read_from_file};
use polars::prelude::*;
use pretty_assertions::assert_eq;
use rstest::rstest;
use std::fs;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

fn unique_path(ext: Option<&str>) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock before unix epoch")
        .as_nanos();
    match ext {
        None => std::env::temp_dir().join(format!("io_test_{nanos}")),
        Some(ext) => std::env::temp_dir().join(format!("io_test_{nanos}.{ext}")),
    }
}

fn write_text_file(path: &PathBuf, contents: &str) {
    fs::write(path, contents).expect("failed to write test input file");
}

fn collect_df(lf: LazyFrame) -> DataFrame {
    lf.collect().expect("failed to collect lazyframe")
}

#[rstest]
#[case::parse_csv(Some(b','), "a,b\n1,2\n", Some("csv"))]
#[case::parse_tsv(Some(b'\t'), "a\tb\n1\t2\n", Some("tsv"))]
#[case::parse_tsv_default_sep(None, "a\tb\n1\t2\n", Some("tsv"))]
#[should_panic(expected = "Unsupported separator")]
#[case::fail_unsupported_sep(Some(b';'), "a\tb\n1\t2\n", Some("csv"))]
#[case::parse_tsv_no_ext(None, "a\tb\n1\t2\n", None)]
fn read_from_file_with_sep(#[case] sep: Option<u8>, #[case] text: &str, #[case] ext: Option<&str>) {
    let path = unique_path(ext);
    write_text_file(&path, text);

    let df = collect_df(read_from_file(path.clone(), sep));
    fs::remove_file(path).ok();

    assert_eq!(df.shape(), (1, 2));
    assert_eq!(df.get_column_names()[0].as_str(), "a");
    assert_eq!(df.get_column_names()[1].as_str(), "b");
}

#[test]
fn write_to_flat_or_stdout_separator_override_beats_extension() {
    let out = unique_path(Some("tsv")); // extension says TSV
    let df = df![
        "a" => &[1],
        "b" => &[2]
    ]
    .expect("failed to build dataframe");

    // Explicit CSV separator should override .tsv extension.
    df.write_to_flat_or_stdout(out.clone(), Some(b','));

    let written = fs::read_to_string(&out).expect("failed to read output file");
    fs::remove_file(out).ok();

    assert!(written.contains("1,2"), "expected csv row, got: {written}");
}

#[test]
fn write_to_flat_or_stdout_uses_extension_when_separator_none() {
    let out = unique_path(Some("tsv"));
    let df = df![
        "a" => &[1],
        "b" => &[2]
    ]
    .expect("failed to build dataframe");

    df.write_to_flat_or_stdout(out.clone(), None);

    let written = fs::read_to_string(&out).expect("failed to read output file");
    fs::remove_file(out).ok();

    assert!(written.contains("1\t2"), "expected tsv row, got: {written}");
}

#[test]
#[should_panic(expected = "Unsupported separator")]
fn write_to_flat_or_stdout_unsupported_separator_panics() {
    let out = unique_path(Some("csv"));
    let df = df!["a" => &[1]].expect("failed to build dataframe");

    df.write_to_flat_or_stdout(out, Some(b';'));
}

#[test]
fn tsv_parsing_regression_tab_delimited_file_produces_two_columns() {
    let path = unique_path(Some("tsv"));
    write_text_file(&path, "col1\tcol2\nv1\tv2\n");

    let df = collect_df(read_from_file(path.clone(), Some(b'\t')));
    fs::remove_file(path).ok();

    // This test is intentionally strict: if TSV is parsed as CSV, you'll get one column.
    assert_eq!(
        df.width(),
        2,
        "TSV should parse into 2 columns; got {} columns with names {:?}",
        df.width(),
        df.get_column_names()
    );
}
