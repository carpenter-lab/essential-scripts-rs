use assert_cmd::cargo::cargo_bin_cmd;
use rstest::rstest;
use rstest_reuse::{apply, template};
use std::fs;
use tempfile::tempdir;

#[rstest]
#[case("--version")]
#[case("--help")]
#[case("help")]
#[should_panic]
#[case("does-not-exist")]
fn base_command_works(#[case] arg: &str) {
    cargo_bin_cmd!("essential-scripts-rs").arg(arg).unwrap();
}

#[template]
#[rstest]
#[case::aggregate_cell_ranger_tcr("aggregate-cell-ranger-tcr")]
#[case::split_sample_id("split-sample-id")]
#[case::split_cdr3_seq("split-cdr3-seq")]
#[case::reformat_plate_reader_data("reformat-plate-reader-data")]
#[case::copy_cell_ranger_outs("copy-cell-ranger-outs")]
#[case::score_tcr_alignments("score-tcr-alignments")]
#[case::geo_fastq("geo-fastq")]
#[case::run_enrichr("run-enrichr")]
fn command_cases(#[case] arg: &str) {}

#[apply(command_cases)]
fn subcommand_help_works(#[case] arg: &str) {
    cargo_bin_cmd!("essential-scripts-rs")
        .args([arg, "--help"])
        .assert()
        .success();
    cargo_bin_cmd!("essential-scripts-rs")
        .args(["help", arg])
        .assert()
        .success();
}

#[apply(command_cases)]
fn subcommand_fails_missing_input(#[case] arg: &str) {
    cargo_bin_cmd!("essential-scripts-rs")
        .args([arg])
        .assert()
        .failure();
}

#[rstest]
#[cfg(not(feature = "tcr"))]
#[case::score_tcr_alignments("score-tcr-alignments")]
#[cfg(not(feature = "enrichment"))]
#[case::run_enrichr("run-enrichr")]
fn subcommand_help_prints_installation_help(#[case] arg: &str) {
    let output = cargo_bin_cmd!("essential-scripts-rs")
        .args([arg, "--help"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();

    let stdout = String::from_utf8_lossy(&output);
    assert!(stdout.contains("Requires re-installing with --features"));
}

#[test]
fn reformat_plate_reader_data_subcommand_fails_for_missing_input_file() {
    let tmp = tempdir().unwrap();
    let missing_input = tmp.path().join("missing.xlsx");
    let output_dir = tmp.path().join("output");

    cargo_bin_cmd!("essential-scripts-rs")
        .args([
            "reformat-plate-reader-data",
            missing_input.to_str().unwrap(),
            output_dir.to_str().unwrap(),
        ])
        .assert()
        .failure();
}

#[test]
fn split_sample_id_cli_runs() {
    let tmp = tempdir().unwrap();
    let input = tmp.path().join("input.csv");
    let output = tmp.path().join("output.csv");

    fs::write(&input, "subject:condition,other\nsubj1:cond1,x\n").unwrap();

    cargo_bin_cmd!("essential-scripts-rs")
        .args([
            "split-sample-id",
            input.to_str().unwrap(),
            output.to_str().unwrap(),
            "--column-name",
            "subject:condition",
        ])
        .assert()
        .success();

    assert!(output.is_file());
}

#[test]
fn copy_cellranger_outs_check_mode_runs() {
    let tmp = tempdir().unwrap();
    let base = tmp.path().join("base");
    fs::create_dir_all(&base).unwrap();

    let pipestance = base.join("p1");
    fs::create_dir_all(&pipestance).unwrap();

    // check mode only needs a directory tree to scan
    cargo_bin_cmd!("essential-scripts-rs")
        .args([
            "copy-cell-ranger-outs",
            "--base-path",
            base.to_str().unwrap(),
            "--check",
        ])
        .assert()
        .failure();

    fs::write(pipestance.join("p1.mri.tgz"), "some content").expect("Failed to create p1.mri.tgz");
    cargo_bin_cmd!("essential-scripts-rs")
        .args([
            "copy-cell-ranger-outs",
            "--base-path",
            base.to_str().unwrap(),
            "--check",
        ])
        .assert()
        .success();
}
