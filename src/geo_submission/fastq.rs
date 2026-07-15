use crate::geo_submission::build_reports::prepare_paths_report;
use crate::geo_submission::traits::*;
use crate::geo_submission::{Progress, build_reports, make_progress_bar};
use indicatif::ProgressBar;
use regex::Regex;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Arc;
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
) {
    let progress = progress.unwrap_or(Progress::Progress);

    match scan_fastq_directories(input_directories, parallel_md5, jobs, progress) {
        Ok(fastq_files) => {
            if fastq_files.is_empty() {
                eprintln!(
                    "No FastQ files found in directories: {}",
                    input_directories
                        .iter()
                        .map(|p| p.display().to_string())
                        .collect::<Vec<_>>()
                        .join(", ")
                );
                return;
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
            if let Err(e) = build_reports::write_output(&paired_report, paired_output) {
                eprintln!("Error writing paired output: {e}");
            }

            // Group by sample only
            let sample_groups = group_by_sample(&fastq_files);
            let sample_report = generate_sample_report(&sample_groups);
            if let Err(e) = build_reports::write_output(&sample_report, sample_output) {
                eprintln!("Error writing sample output: {e}");
            }

            // MD5 manifest
            let md5_report = build_reports::generate_md5_report(&fastq_files);
            if let Err(e) = build_reports::write_output(&md5_report, md5_output) {
                eprintln!("Error writing MD5 output: {e}");
            }
        }
        Err(e) => {
            eprintln!("Error scanning directory: {e}");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

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
}
