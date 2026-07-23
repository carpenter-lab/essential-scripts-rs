use crate::geo_submission::build_reports::prepare_paths_report;
use crate::geo_submission::helper::make_progress_bar;
use crate::geo_submission::traits::*;
use crate::geo_submission::{Progress, build_reports};
use indicatif::ProgressBar;
use regex::Regex;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::{fs, io};
use walkdir::WalkDir;

#[derive(Debug, Clone)]
pub(super) struct FastqFile {
    path: PathBuf,
    sample_id: String,
    lane: String,
    read_number: String,
    md5: String,
}

impl Md5Record for FastqFile {
    fn filename(&self) -> String {
        prepare_paths_report(self)
    }
    fn md5_str(&self) -> &str {
        &self.md5
    }
}

impl HasPath for FastqFile {
    fn path(&self) -> &Path {
        &self.path
    }
}

impl HasPath for &FastqFile {
    fn path(&self) -> &Path {
        &self.path
    }
}

impl FromPathWithMd5 for FastqFile {
    fn from_path_with_md5(
        path: PathBuf,
        pb: Option<&Arc<ProgressBar>>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let filename = path
            .file_name()
            .and_then(|s| s.to_str())
            .ok_or("bad filename")?;
        let (sample_id, lane, read_number) =
            parse_fastq_filename(filename).ok_or("unparsable FastQ name")?;
        let md5 = crate::geo_submission::md5::compute_md5_with_progress(&path, pb)?;
        Ok(FastqFile {
            path,
            sample_id,
            lane,
            read_number,
            md5,
        })
    }
}

/// Parse bcl2fastq filename format: {sample_id}_S{sample_number}_{lane}_{read}_001.fastq.gz
/// Where read can be I1, I2, R1, or R2
fn parse_fastq_filename(filename: &str) -> Option<(String, String, String)> {
    let re = Regex::new(r"^(.+?)_S\d+_(L\d{3})_([IR][12])_001\.fastq(?:\.gz)?$").ok()?;
    let caps = re.captures(filename)?;
    Some((
        caps[1].to_string(),
        caps[2].to_string(),
        caps[3].to_string(),
    ))
}

/// Generate paired files report (`sample_lane` -> R1, R2 pairs)
fn generate_paired_report(groups: &BTreeMap<(String, String), Vec<&FastqFile>>) -> String {
    let mut output = String::new();
    output.push_str("# Paired FastQ Files by Sample and Lane\n");
    output.push_str("# Format: sample_lane: R1=/path/to/R1.fastq.gz, R2=/path/to/R2.fastq.gz\n\n");

    for ((_sample_id, _lane), files) in groups {
        //output.push_str(&format!("{}_{}", sample_id, lane));

        let mut r1_path = String::new();
        let mut r2_path = String::new();
        let mut i1_path = String::new();
        let mut i2_path = String::new();

        for file in files {
            let path_str = prepare_paths_report(file);
            if file.read_number == "R1" {
                r1_path = path_str;
            } else if file.read_number == "R2" {
                r2_path = path_str;
            } else if file.read_number == "I1" {
                i1_path = path_str;
            } else if file.read_number == "I2" {
                i2_path = path_str;
            }
        }

        if !r1_path.is_empty() {
            output.push_str(&format!("\t{}", r1_path));
        }
        if !r2_path.is_empty() {
            output.push_str(&format!("\t{}", r2_path));
        }
        if !i1_path.is_empty() {
            output.push_str(&format!("\t{}", i1_path));
        }
        if !i2_path.is_empty() {
            output.push_str(&format!("\t{}", i2_path));
        }
        output.push('\n');
    }

    output
}

/// Group files by (`sample_id`, `lane`)
fn group_by_lane(files: &[FastqFile]) -> BTreeMap<(String, String), Vec<&FastqFile>> {
    let mut groups = BTreeMap::new();
    for file in files {
        let key = (file.sample_id.clone(), file.lane.clone());
        groups.entry(key).or_insert_with(Vec::new).push(file);
    }
    groups
}

/// Group files by `sample_id` only
fn group_by_sample(files: &[FastqFile]) -> BTreeMap<String, Vec<&FastqFile>> {
    let mut groups = BTreeMap::new();
    for file in files {
        groups
            .entry(file.sample_id.clone())
            .or_insert_with(Vec::new)
            .push(file);
    }
    groups
}

fn scan_fastq_directories(
    dirs: &[PathBuf],
    parallel: &bool,
    jobs: &usize,
    progress: Progress,
) -> Result<Vec<FastqFile>, Box<dyn std::error::Error>> {
    let mut fastq_paths = Vec::new();

    for d in dirs {
        for entry in WalkDir::new(d)
            .into_iter()
            .filter_map(Result::ok)
            .filter(|e| {
                e.file_name().to_str().is_some_and(|s| {
                    Path::new(s)
                        .extension()
                        .is_some_and(|ext| ext.eq_ignore_ascii_case("fastq"))
                        || s.ends_with(".fastq.gz")
                })
            })
        {
            let path = entry.path();
            fastq_paths.push(path.to_path_buf());
        }
    }
    let total_bytes: u64 = fastq_paths
        .iter()
        .map(|path| fs::metadata(path).map_or(0, |m| m.len()))
        .sum();

    let pb = make_progress_bar(total_bytes, progress)?;
    let mut fastq_files: Vec<FastqFile> =
        build_reports::build_records_from_paths(fastq_paths, Some(&pb), parallel, jobs);
    fastq_files.sort_by(|a, b| {
        a.sample_id
            .cmp(&b.sample_id)
            .then_with(|| a.lane.cmp(&b.lane))
            .then_with(|| a.read_number.cmp(&b.read_number))
    });
    Ok(fastq_files)
}

/// Generate per-sample report (all files for each sample on one line)
fn generate_sample_report(groups: &BTreeMap<String, Vec<&FastqFile>>) -> String {
    let mut output = String::new();
    output.push_str("# All FastQ Files per Sample\n");
    output.push_str("# Format: sample_id: /path/to/file1.fastq.gz /path/to/file2.fastq.gz ...\n\n");

    for (sample_id, files) in groups {
        output.push_str(&format!("{sample_id}\t"));
        let paths: Vec<String> = files.iter().map(prepare_paths_report).collect();
        output.push_str(&paths.join("\t"));
        output.push('\n');
    }

    output
}

pub(super) fn match_fastq(
    input_directories: &[PathBuf],
    paired_output: Option<&PathBuf>,
    sample_output: Option<&PathBuf>,
    md5_output: Option<&PathBuf>,
    parallel_md5: &bool,
    jobs: &usize,
    progress: Option<Progress>,
) -> Result<(), Box<dyn std::error::Error>> {
    let progress = progress.unwrap_or(Progress::Progress);

    match scan_fastq_directories(input_directories, parallel_md5, jobs, progress) {
        Ok(fastq_files) => {
            if fastq_files.is_empty() {
                let err_msg = format!(
                    "No FastQ files found in directories: {}",
                    input_directories
                        .iter()
                        .map(|p| p.display().to_string())
                        .collect::<Vec<_>>()
                        .join(", ")
                );
                return Err(io::Error::new(io::ErrorKind::NotFound, err_msg.as_str()).into());
            }

            println!(
                "Found {} FastQ files in {}",
                fastq_files.len(),
                input_directories
                    .iter()
                    .map(|p| p.display().to_string())
                    .collect::<Vec<_>>()
                    .join(", ")
            );

            // Group by (sample, lane)
            let lane_groups = group_by_lane(&fastq_files);
            let paired_report = generate_paired_report(&lane_groups);
            build_reports::write_output(&paired_report, paired_output)?;

            // Group by sample only
            let sample_groups = group_by_sample(&fastq_files);
            let sample_report = generate_sample_report(&sample_groups);
            build_reports::write_output(&sample_report, sample_output)?;

            // MD5 manifest
            let md5_report = build_reports::generate_md5_report(&fastq_files);
            build_reports::write_output(&md5_report, md5_output)?;
            Ok(())
        }
        Err(e) => {
            eprintln!("Error scanning directory: {e}");
            Err(e)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::rstest;
    use tempfile::tempdir;

    fn write_file(path: &Path, content: &[u8]) {
        fs::write(path, content).expect("failed to write test file");
    }

    #[rstest]
    #[case::with_gz("sample01_S1_L001_R1_001.fastq.gz", Some(("sample01".to_string(), "L001".to_string(), "R1".to_string())
    ))]
    #[case::without_gz("sample01_S1_L001_R1_001.fastq", Some(("sample01".to_string(), "L001".to_string(), "R1".to_string())
    ))]
    #[case::with_underscore_in_sample("sample_01_with_underscores_S2_L002_R1_001.fastq", Some(("sample_01_with_underscores".to_string(), "L002".to_string(), "R1".to_string())
    ))]
    #[case::index_i1("sample01_S1_L001_I1_001.fastq", Some(("sample01".to_string(), "L001".to_string(), "I1".to_string())
    ))]
    #[case::index_i2("sample01_S1_L001_I2_001.fastq", Some(("sample01".to_string(), "L001".to_string(), "I2".to_string())
    ))]
    #[case::read_r2("sample01_S1_L001_R2_001.fastq", Some(("sample01".to_string(), "L001".to_string(), "R2".to_string())
    ))]
    #[case::read_r1("sample01_S1_L001_R1_001.fastq", Some(("sample01".to_string(), "L001".to_string(), "R1".to_string())
    ))]
    #[case::invalid_lane_digits("sample01_S1_L01_R1_001.fastq", None)]
    #[case::invalid_read_number("sample01_S1_L001_R3_001.fastq", None)]
    #[case::invalid_suffix_counter("sample01_S1_L001_R1_002.fastq", None)]
    #[case::invalid_extension("sample01_S1_L001_R1_001.fq.gz", None)]
    #[case::invalid_filename("invalid_filename.fastq.gz", None)]
    fn test_parse_fastq_filename(
        #[case] filename: &str,
        #[case] exp: Option<(String, String, String)>,
    ) {
        let result = parse_fastq_filename(filename);
        assert_eq!(result, exp);
    }
    #[test]
    fn test_group_by_lane() {
        let files = vec![
            FastqFile {
                path: PathBuf::from("/path/file1.fastq.gz"),
                sample_id: "sample01".to_string(),
                lane: "L001".to_string(),
                read_number: "R1".to_string(),
                md5: "abc123".to_string(),
            },
            FastqFile {
                path: PathBuf::from("/path/file2.fastq.gz"),
                sample_id: "sample01".to_string(),
                lane: "L001".to_string(),
                read_number: "R2".to_string(),
                md5: "def456".to_string(),
            },
            FastqFile {
                path: PathBuf::from("/path/file3.fastq.gz"),
                sample_id: "sample01".to_string(),
                lane: "L002".to_string(),
                read_number: "R1".to_string(),
                md5: "ghi789".to_string(),
            },
        ];

        let groups = group_by_lane(&files);
        assert_eq!(groups.len(), 2);
        assert!(groups.contains_key(&("sample01".to_string(), "L001".to_string())));
        assert!(groups.contains_key(&("sample01".to_string(), "L002".to_string())));
    }

    #[test]
    fn test_group_by_sample() {
        let files = vec![
            FastqFile {
                path: PathBuf::from("/path/file1.fastq.gz"),
                sample_id: "sample01".to_string(),
                lane: "L001".to_string(),
                read_number: "R1".to_string(),
                md5: "abc123".to_string(),
            },
            FastqFile {
                path: PathBuf::from("/path/file2.fastq.gz"),
                sample_id: "sample02".to_string(),
                lane: "L001".to_string(),
                read_number: "R1".to_string(),
                md5: "def456".to_string(),
            },
        ];

        let groups = group_by_sample(&files);
        assert_eq!(groups.len(), 2);
        assert!(groups.contains_key("sample01"));
        assert!(groups.contains_key("sample02"));
    }

    #[rstest]
    fn test_from_path_with_md5_success() {
        let dir = tempdir().expect("failed to create tempdir");
        let path = dir.path().join("sample01_S1_L001_R1_001.fastq");
        write_file(&path, b"ACGT\n");

        let record =
            FastqFile::from_path_with_md5(path.clone(), None).expect("record creation failed");
        assert_eq!(record.path, path);
        assert_eq!(record.sample_id, "sample01");
        assert_eq!(record.lane, "L001");
        assert_eq!(record.read_number, "R1");
        assert_eq!(record.md5.len(), 32);
    }

    #[test]
    fn test_from_path_with_md5_bad_filename_error() {
        let path = PathBuf::from("/");
        let err =
            FastqFile::from_path_with_md5(path, None).expect_err("expected bad filename error");
        assert_eq!(err.to_string(), "bad filename");
    }

    #[test]
    fn test_from_path_with_md5_unparsable_name_error() {
        let dir = tempdir().expect("failed to create tempdir");
        let path = dir.path().join("not_a_fastq_name.fastq");
        write_file(&path, b"ACGT\n");

        let err =
            FastqFile::from_path_with_md5(path, None).expect_err("expected unparsable name error");
        assert_eq!(err.to_string(), "unparsable FastQ name");
    }

    #[test]
    fn test_generate_paired_report_includes_reads() {
        let files = vec![
            FastqFile {
                path: PathBuf::from("/tmp/sample01_S1_L001_R1_001.fastq.gz"),
                sample_id: "sample01".to_string(),
                lane: "L001".to_string(),
                read_number: "R1".to_string(),
                md5: "a".to_string(),
            },
            FastqFile {
                path: PathBuf::from("/tmp/sample01_S1_L001_R2_001.fastq.gz"),
                sample_id: "sample01".to_string(),
                lane: "L001".to_string(),
                read_number: "R2".to_string(),
                md5: "b".to_string(),
            },
            FastqFile {
                path: PathBuf::from("/tmp/sample01_S1_L001_I1_001.fastq.gz"),
                sample_id: "sample01".to_string(),
                lane: "L001".to_string(),
                read_number: "I1".to_string(),
                md5: "c".to_string(),
            },
            FastqFile {
                path: PathBuf::from("/tmp/sample01_S1_L001_I2_001.fastq.gz"),
                sample_id: "sample01".to_string(),
                lane: "L001".to_string(),
                read_number: "I2".to_string(),
                md5: "d".to_string(),
            },
            FastqFile {
                path: PathBuf::from("/tmp/sample02_S1_L001_R1_001.fastq.gz"),
                sample_id: "sample02".to_string(),
                lane: "L001".to_string(),
                read_number: "R1".to_string(),
                md5: "e".to_string(),
            },
        ];

        let lane_groups = group_by_lane(&files);
        let report = generate_paired_report(&lane_groups);
        let lines: Vec<&str> = report.lines().collect();

        assert!(lines.contains(&"\tsample01_S1_L001_R1_001.fastq.gz\tsample01_S1_L001_R2_001.fastq.gz\tsample01_S1_L001_I1_001.fastq.gz\tsample01_S1_L001_I2_001.fastq.gz"));
        assert!(lines.contains(&"\tsample02_S1_L001_R1_001.fastq.gz"));
    }

    #[test]
    fn test_generate_sample_report_groups_paths_per_sample() {
        let files = vec![
            FastqFile {
                path: PathBuf::from("/tmp/sample01_S1_L001_R1_001.fastq.gz"),
                sample_id: "sample01".to_string(),
                lane: "L001".to_string(),
                read_number: "R1".to_string(),
                md5: "a".to_string(),
            },
            FastqFile {
                path: PathBuf::from("/tmp/sample01_S1_L001_R2_001.fastq.gz"),
                sample_id: "sample01".to_string(),
                lane: "L001".to_string(),
                read_number: "R2".to_string(),
                md5: "b".to_string(),
            },
            FastqFile {
                path: PathBuf::from("/tmp/sample02_S1_L001_R1_001.fastq.gz"),
                sample_id: "sample02".to_string(),
                lane: "L001".to_string(),
                read_number: "R1".to_string(),
                md5: "c".to_string(),
            },
        ];

        let sample_groups = group_by_sample(&files);
        let report = generate_sample_report(&sample_groups);
        assert!(report.contains(
            "sample01\tsample01_S1_L001_R1_001.fastq.gz\tsample01_S1_L001_R2_001.fastq.gz"
        ));
        assert!(report.contains("sample02\tsample02_S1_L001_R1_001.fastq.gz"));
    }

    #[test]
    fn test_scan_fastq_directories_finds_and_sorts_fastq_files() {
        let dir = tempdir().expect("failed to create tempdir");
        let nested = dir.path().join("nested");
        fs::create_dir(&nested).expect("failed to create nested dir");

        let paths = vec![
            dir.path().join("sample02_S1_L001_R2_001.fastq"),
            dir.path().join("sample01_S1_L001_R2_001.fastq.gz"),
            nested.join("sample01_S1_L001_R1_001.fastq"),
        ];
        for path in &paths {
            write_file(path, b"ACGT\n");
        }
        write_file(&dir.path().join("ignore.txt"), b"not fastq");

        let records = scan_fastq_directories(
            &[dir.path().to_path_buf()],
            &false,
            &1,
            Progress::NoProgress,
        )
        .expect("scan should succeed");

        assert_eq!(records.len(), 3);
        let names: Vec<String> = records.iter().map(|f| prepare_paths_report(f)).collect();
        assert_eq!(
            names,
            vec![
                "sample01_S1_L001_R1_001.fastq".to_string(),
                "sample01_S1_L001_R2_001.fastq.gz".to_string(),
                "sample02_S1_L001_R2_001.fastq".to_string()
            ]
        );
    }

    #[test]
    fn test_match_fastq_returns_not_found_for_empty_input() {
        let dir = tempdir().expect("failed to create tempdir");
        let err = match_fastq(
            &[dir.path().to_path_buf()],
            None,
            None,
            None,
            &false,
            &1,
            Some(Progress::NoProgress),
        )
        .expect_err("expected not found error");

        assert!(
            err.to_string()
                .contains("No FastQ files found in directories")
        );
    }

    #[test]
    fn test_match_fastq_writes_all_reports() {
        let dir = tempdir().expect("failed to create tempdir");
        let input = dir.path().join("input");
        fs::create_dir(&input).expect("failed to create input dir");

        let r1 = input.join("sample01_S1_L001_R1_001.fastq.gz");
        let r2 = input.join("sample01_S1_L001_R2_001.fastq.gz");
        let i1 = input.join("sample01_S1_L001_I1_001.fastq.gz");
        write_file(&r1, b"AAA\n");
        write_file(&r2, b"CCC\n");
        write_file(&i1, b"GGG\n");

        let paired_output = dir.path().join("paired.tsv");
        let sample_output = dir.path().join("sample.tsv");
        let md5_output = dir.path().join("md5.tsv");

        match_fastq(
            &[input],
            Some(&paired_output),
            Some(&sample_output),
            Some(&md5_output),
            &false,
            &1,
            Some(Progress::NoProgress),
        )
        .expect("match_fastq should succeed");

        let paired = fs::read_to_string(&paired_output).expect("failed to read paired output");
        let sample = fs::read_to_string(&sample_output).expect("failed to read sample output");
        let md5 = fs::read_to_string(&md5_output).expect("failed to read md5 output");

        assert!(paired.contains("sample01_S1_L001_R1_001.fastq.gz"));
        assert!(paired.contains("sample01_S1_L001_R2_001.fastq.gz"));
        assert!(paired.contains("sample01_S1_L001_I1_001.fastq.gz"));
        assert!(sample.contains("sample01\tsample01_S1_L001_I1_001.fastq.gz"));
        assert!(sample.contains("sample01_S1_L001_R1_001.fastq.gz"));
        assert!(sample.contains("sample01_S1_L001_R2_001.fastq.gz"));

        let md5_lines_count = md5
            .lines()
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .count();
        assert_eq!(md5_lines_count, 3);
    }
}
