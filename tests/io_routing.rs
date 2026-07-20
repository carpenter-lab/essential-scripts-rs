use essential_scripts_rs::io::{WriteToCsvOrStdout, read_from_file};
use polars::prelude::*;
use pretty_assertions::assert_eq;
use rstest::rstest;
use rstest_reuse::{self, *};
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

#[template]
#[rstest]
#[case(Some("tsv"), Some(b','), "1,2")]
#[case(Some("csv"), Some(b','), "1,2")]
#[case(Some("tsv"), Some(b'\t'), "1\t2")]
#[case(Some("tsv"), None, "1\t2")]
#[should_panic(expected = "Unsupported separator")]
#[case(Some("csv"), Some(b';'), "1;2")]
fn write_to_file_template(
    #[case] path_ext: Option<&str>,
    #[case] file_sep: Option<u8>,
    #[case] exp: &str,
) {
}

#[apply(write_to_file_template)]
fn write_to_flat(#[case] path_ext: Option<&str>, #[case] file_sep: Option<u8>, #[case] exp: &str) {
    let out = unique_path(path_ext);
    let df = df![
        "a" => &[1],
        "b" => &[2]
    ]
    .expect("failed to build dataframe");

    df.write_to_flat_or_stdout(out.clone(), file_sep);
    let written = fs::read_to_string(&out).expect("failed to read output file");
    assert!(written.contains(exp), "expected row {exp}, got: {written}");
}

#[apply(write_to_file_template)]
fn write_to_flat_lazy(
    #[case] path_ext: Option<&str>,
    #[case] file_sep: Option<u8>,
    #[case] exp: &str,
) {
    let out = unique_path(path_ext);
    let df = df![
        "a" => &[1],
        "b" => &[2]
    ]
    .expect("failed to build dataframe")
    .lazy();

    df.write_to_flat_or_stdout(out.clone(), file_sep);
    let written = fs::read_to_string(&out).expect("failed to read output file");
    assert!(written.contains(exp), "expected row {exp}, got: {written}");
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
