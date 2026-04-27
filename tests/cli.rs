use assert_cmd::cargo::cargo_bin_cmd;
use rstest::rstest;
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

#[rstest]
#[case::aggregate_cell_ranger_tcr("aggregate-cell-ranger-tcr")]
#[case::split_sample_id("split-sample-id")]
#[case::split_cdr3_seq("split-cdr3-seq")]
#[case::reformat_plate_reader_data("reformat-plate-reader-data")]
#[case::copy_cell_ranger_outs("copy-cell-ranger-outs")]
#[case::score_tcr_alignments("score-tcr-alignments")]
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
        .success();
}
