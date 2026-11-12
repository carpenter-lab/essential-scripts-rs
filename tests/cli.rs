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
#[case::run_gliph("run-gliph")]
#[case::tcr_match("tcr-match")]
fn command_cases(#[case] arg: &str) {}

#[apply(command_cases)]
fn subcommand_help_works(#[case] arg: &str) {
    let output = cargo_bin_cmd!("essential-scripts-rs")
        .args([arg, "--help"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    cargo_bin_cmd!("essential-scripts-rs")
        .args(["help", arg])
        .assert()
        .success();
    if arg == "run-gliph" {
        assert!(
            String::from_utf8_lossy(&output).contains("--hla-file"),
            "run-gliph help should document optional HLA input"
        );
    }
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
fn split_sample_id_cli_reads_stdin() {
    let tmp = tempdir().unwrap();
    let output = tmp.path().join("output.csv");

    cargo_bin_cmd!("essential-scripts-rs")
        .args([
            "split-sample-id",
            "-",
            output.to_str().unwrap(),
            "--column-name",
            "subject:condition",
        ])
        .write_stdin("subject:condition;other\nsubj1:cond1;x\n")
        .assert()
        .success();

    let written = fs::read_to_string(&output).unwrap();
    assert!(written.contains("subj1"), "got: {written}");
    assert!(written.contains("cond1"), "got: {written}");
}

#[test]
fn tcr_match_cli_runs() {
    let tmp = tempdir().unwrap();
    let input_path = tmp.path().join("input.txt");

    fs::write(&input_path, "CASSLGF\n").unwrap();

    let output = cargo_bin_cmd!("essential-scripts-rs")
        .args([
            "tcr-match",
            input_path.to_str().unwrap(),
            "-d",
            "iedb",
            "-s",
            "0.9",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();

    let stdout = String::from_utf8_lossy(&output);
    assert!(stdout.contains(
        "trimmed_input_sequence\tmatch_sequence\tscore\treceptor_group\tepitope\tantigen\torganism"
    ));
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

/// Copies the trimmed immGLIPH reference fixtures into a fresh cache directory so
/// the clustering commands run offline.
#[cfg(feature = "base_cmd")]
fn reference_cache() -> tempfile::TempDir {
    let cache = tempdir().unwrap();
    let fixtures =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/gliph_reference");
    for name in ["reference_list.RData", "ref_cluster_sizes.RData"] {
        fs::copy(fixtures.join(name), cache.path().join(name)).unwrap();
    }
    cache
}

/// Synthetic chains: three families of near-identical CDR3βs shared across subjects,
/// plus unrelated singletons.
#[cfg(feature = "base_cmd")]
fn write_chains(path: &std::path::Path) {
    let families = [
        ("CASSLGQAYEQYF", "TRBV7-9", "CAVRDSNYQLIW", "TRAV12-2"),
        ("CASSPDRGNTEAFF", "TRBV6-5", "CAMREGGSQGNLIF", "TRAV14/DV4"),
        ("CSARDRTGNGYTF", "TRBV20-1", "CAASGGSYIPTF", "TRAV29/DV5"),
    ];
    let substitutions = ['A', 'G', 'S', 'T', 'N', 'D'];
    let mut rows = String::from("CDR3b\tTRBV\tTRBJ\tCDR3a\tTRAV\tTRAJ\tsubject:condition\tcount\n");
    for (family, (cdr3b, trbv, cdr3a, trav)) in families.iter().enumerate() {
        for (i, aa) in substitutions.iter().enumerate() {
            let mut seq: Vec<char> = cdr3b.chars().collect();
            let position = 5 + (i % 3);
            seq[position] = *aa;
            let seq: String = seq.into_iter().collect();
            let subject = format!("S{:02}", (family * 7 + i) % 5 + 1);
            rows.push_str(&format!(
                "{seq}\t{trbv}\tTRBJ2-7\t{cdr3a}\t{trav}\tTRAJ33\t{subject}:D0\t{}\n",
                i + 1
            ));
        }
    }
    for (i, seq) in [
        "CASSQETQYF",
        "CASRWGGEQFF",
        "CSVEGTSGRADTQYF",
        "CASSYSTDTQYF",
    ]
    .iter()
    .enumerate()
    {
        rows.push_str(&format!(
            "{seq}\tTRBV5-1\tTRBJ1-1\t\t\t\tS0{}:D7\t1\n",
            i + 1
        ));
    }
    fs::write(path, rows).unwrap();
}

#[cfg(feature = "base_cmd")]
fn write_hla(path: &std::path::Path) {
    let mut rows = String::from(
        "PID\tHLA-DPA1(L)\tHLA-DPA1(R)\tHLA-DPB1(L)\tHLA-DPB1(R)\tHLA-DQA1(L)\tHLA-DQA1(R)\tHLA-DQB1(L)\tHLA-DQB1(R)\tHLA-DRB1(L)\tHLA-DRB1(R)\tHLA-DRB3/4/5(L)\tHLA-DRB3/4/5(R)\n",
    );
    for subject in 1..=5 {
        let drb1 = if subject % 2 == 0 {
            "DRB1*04:01"
        } else {
            "DRB1*07:01"
        };
        rows.push_str(&format!(
            "S{subject:02}\tDPA1*01:03\tDPA1*02:01\tDPB1*04:01\tNA\tDQA1*03:01\tDQA1*05:01\tDQB1*03:02\t\t{drb1}\tDRB1*15:01\tDRB4*01:03\tDRB5*01:01\n"
        ));
    }
    fs::write(path, rows).unwrap();
}

#[cfg(feature = "base_cmd")]
#[rstest]
#[case::gliph2(&["run-gliph"][..])]
#[case::gliph1(&["run-gliph", "--method", "gliph1"][..])]
#[case::gliph2_blosum(&["run-gliph", "--sub-matrix", "blosum62"][..])]
#[case::tcrdist_inferred(&["run-tcr-dist"][..])]
#[case::tcrdist_beta_single(&["run-tcr-dist", "--distance-mode", "beta", "--cluster-method", "single-linkage"][..]
)]
#[case::tcrdist_graph(&["run-tcr-dist", "--distance-mode", "alpha-beta", "--cluster-method", "connected-components"][..]
)]
fn clustering_cli_runs_offline(#[case] command: &[&str]) {
    let tmp = tempdir().unwrap();
    let cache = reference_cache();
    let chains = tmp.path().join("chains.tsv");
    let hla = tmp.path().join("hla.tsv");
    let output = tmp.path().join("out.tsv");
    write_chains(&chains);
    write_hla(&hla);

    cargo_bin_cmd!("essential-scripts-rs")
        .args(command)
        .args([
            chains.to_str().unwrap(),
            output.to_str().unwrap(),
            "--hla-file",
            hla.to_str().unwrap(),
            "--reference-cache-dir",
            cache.path().to_str().unwrap(),
            "--sim-depth",
            "20",
            "--cores",
            "2",
        ])
        .assert()
        .success();

    let written = fs::read_to_string(&output).unwrap();
    let header = written.lines().next().expect("output has a header");
    assert!(header.contains("TcRb"), "unexpected header: {header}");
    let columns: Vec<&str> = header.split('\t').collect();
    let alpha_pattern = columns
        .iter()
        .position(|c| *c == "pattern_alpha")
        .unwrap_or_else(|| panic!("no pattern_alpha column: {header}"));
    assert!(
        written.lines().skip(1).any(|line| !line
            .split('\t')
            .nth(alpha_pattern)
            .unwrap_or("")
            .is_empty()),
        "no alpha motif written:\n{written}"
    );
    assert!(
        written.lines().count() > 1,
        "no clusters written:\n{written}"
    );
}

#[cfg(feature = "base_cmd")]
#[rstest]
#[case::zero_sim_depth(&["--sim-depth", "0"][..])]
#[case::bad_hla_threshold(&["--hla-threshold", "2"][..])]
#[case::other_reference_missing(&["--reference", "human-v1-cd4"][..])]
fn clustering_cli_rejects_invalid_input(#[case] extra: &[&str]) {
    let tmp = tempdir().unwrap();
    let cache = reference_cache();
    let chains = tmp.path().join("chains.tsv");
    write_chains(&chains);

    cargo_bin_cmd!("essential-scripts-rs")
        .args(["run-gliph", chains.to_str().unwrap(), "-"])
        .args(["--reference-cache-dir", cache.path().to_str().unwrap()])
        .args(extra)
        .assert()
        .failure();
}

#[cfg(feature = "base_cmd")]
#[test]
fn clustering_cli_rejects_out_of_range_neighbor_overlap() {
    cargo_bin_cmd!("essential-scripts-rs")
        .args([
            "run-tcr-dist",
            "in.tsv",
            "out.tsv",
            "--min-neighbor-overlap",
            "1.5",
        ])
        .assert()
        .failure();
}
