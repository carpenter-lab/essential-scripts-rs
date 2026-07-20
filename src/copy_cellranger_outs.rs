use clap::{Args, Subcommand};
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

const PATH_TO_OUTS: &str = "outs/per_sample_outs";
const PATH_TO_COUNT_H5: &str = "count/sample_filtered_feature_bc_matrix.h5";
const PATH_TO_COUNT_MEX: &str = "count/sample_filtered_feature_bc_matrix";
const PATH_TO_COUNT_MEX_TAR: &str = "count/sample_filtered_feature_bc_matrix.tar.gz";
const PATH_TO_VDJ_ANN: &str = "vdj_t/filtered_contig_annotations.csv";

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct CopyStats {
    pub pipestances_seen: usize,
    pub samples_seen: usize,
    pub h5_copied: usize,
    pub mex_dirs_created: usize,
    pub mex_files_copied: usize,
    pub mex_dirs_cleaned_up: usize,
    pub vdj_copied: usize,
}

#[derive(Args)]
#[group(multiple = true)]
pub struct PipestanceResults {
    #[arg(
        long,
        default_value_t = false,
        help = "Copy filtered H5 matrix as <sample>.h5"
    )]
    h5: bool,

    #[arg(
        long,
        default_value_t = false,
        help = "Copy filtered MEX directory into <sample>/"
    )]
    mex: bool,

    #[arg(
        long,
        default_value_t = false,
        help = "Copy VDJ annotations as <sample>.csv"
    )]
    vdj: bool,
}

#[derive(Subcommand)]
pub enum Commands {
    #[command(about = "Copy selected outputs from Cell Ranger pipestances")]
    CopyCellRangerOuts {
        #[arg(
            short,
            long,
            help = "Base directory containing Cell Ranger pipestances"
        )]
        base_path: PathBuf,

        #[arg(
            short,
            long,
            required_unless_present_any = ["check"],
            help = "Destination directory to copy to"
        )]
        dest: Option<PathBuf>,

        #[command(flatten)]
        pipestance_results: PipestanceResults,

        #[arg(
            long,
            default_value_t = false,
            help = "Check each pipestance for presence of a *.mri.tgz marker and print results. No copying is performed."
        )]
        check: bool,
    },
}

pub fn handle_command(cmd: Commands) {
    match cmd {
        Commands::CopyCellRangerOuts {
            base_path,
            dest,
            pipestance_results,
            check,
        } => {
            if let Err(e) = copy_cellranger_outs_main(
                &base_path,
                Option::from(&dest),
                &pipestance_results,
                check,
            ) {
                eprintln!("{e}");
            }
        }
    }
}

pub fn copy_cellranger_outs_main(
    base_path: &Path,
    dest: Option<&PathBuf>,
    outs: &PipestanceResults,
    check: bool,
) -> io::Result<CopyStats> {
    let mut stats = CopyStats::default();
    // If in check-only mode, validate pipestances and exit
    if check {
        let _ = check_pipestances(base_path)?;
        return Ok(stats);
    }
    let PipestanceResults { h5, mex, vdj } = outs;

    if !h5 && !mex && !vdj {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Nothing to do: specify at least one of --h5, --mex, --vdj",
        ));
    }

    let Some(dest) = dest else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Missing destination path",
        ));
    };

    // Ensure destination exists
    if let Err(e) = fs::create_dir_all(dest) {
        return Err(io::Error::new(
            e.kind(),
            format!(
                "Failed to create destination directory {}: {}",
                dest.display(),
                e
            ),
        ));
    }
    let entries = fs::read_dir(base_path).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!(
                "Base path does not exist or cannot be read: {}",
                base_path.display()
            ),
        )
    })?;
    for entry in entries.flatten() {
        let dir_path = entry.path();
        if !dir_path.is_dir() {
            continue;
        }
        stats.pipestances_seen += 1;

        match copy_outputs_for_pipestance(&dir_path, dest, *h5, *mex, *vdj) {
            Ok(s) => {
                stats.samples_seen += s.samples_seen;
                stats.h5_copied += s.h5_copied;
                stats.mex_dirs_created += s.mex_dirs_created;
                stats.mex_files_copied += s.mex_files_copied;
                stats.mex_dirs_cleaned_up += s.mex_dirs_cleaned_up;
                stats.vdj_copied += s.vdj_copied;
            }
            Err(e) => eprintln!("{}: {}", dir_path.display(), e),
        }
    }

    Ok(stats)
}

fn has_mri_tgz(pipestance_dir: &Path) -> io::Result<bool> {
    let mut ok = false;
    for entry in fs::read_dir(pipestance_dir)? {
        let entry = entry?;
        let p = entry.path();
        if p.is_file()
            && let Some(name) = p.file_name().and_then(|s| s.to_str())
            && name.ends_with(".mri.tgz")
        {
            ok = true;
            break;
        }
    }

    Ok(ok)
}

pub(crate) fn check_pipestances(base_path: &Path) -> io::Result<usize> {
    let entries = fs::read_dir(base_path)?;
    let mut total = 0usize;
    let mut ok_count = 0usize;
    for entry in entries.flatten() {
        let dir_path = entry.path();
        if !dir_path.is_dir() {
            continue;
        }
        total += 1;
        match has_mri_tgz(&dir_path) {
            Ok(true) => {
                ok_count += 1;
                println!("{}: OK", dir_path.display());
            }
            Ok(false) => {
                eprintln!("{}: MISSING .mri.tgz", dir_path.display());
            }
            Err(e) => {
                eprintln!("{}: ERROR reading directory: {}", dir_path.display(), e);
            }
        }
    }
    println!(
        "Checked {} pipestance(s): {} OK, {} missing",
        total,
        ok_count,
        total - ok_count
    );
    if total == 0 {
        Err(io::Error::new(
            io::ErrorKind::NotFound,
            "No pipestances found",
        ))
    } else if total != ok_count {
        Err(io::Error::other(format!(
            "Some pipestances failed validation: {} OK, {} missing",
            ok_count,
            total - ok_count
        )))
    } else {
        Ok(ok_count)
    }
}

pub fn list_samples(pipestance_dir: &Path) -> io::Result<Vec<PathBuf>> {
    let base = pipestance_dir.join(PATH_TO_OUTS);
    let mut samples = Vec::new();
    if base.is_dir() {
        for entry in fs::read_dir(base)? {
            let entry = entry?;
            let p = entry.path();
            if p.is_dir() {
                samples.push(p);
            }
        }
    }
    Ok(samples)
}

fn copy_file_counting(src: &Path, dst: &Path, counter: &mut usize) {
    if let Err(e) = fs::copy(src, dst) {
        eprintln!(
            "Failed to copy {} to {}: {}",
            src.display(),
            dst.display(),
            e
        );
    } else {
        *counter += 1;
    }
}

fn copy_outputs_for_pipestance(
    pipestance_dir: &Path,
    dest: &Path,
    h5: bool,
    mex: bool,
    vdj: bool,
) -> io::Result<CopyStats> {
    let mut stats = CopyStats::default();

    for sample_dir in list_samples(pipestance_dir)? {
        stats.samples_seen += 1;
        let sample_name = sample_dir
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("sample");

        if h5 {
            let src = sample_dir.join(PATH_TO_COUNT_H5);
            if src.is_file() {
                let dst = dest.join(format!("{sample_name}.h5"));
                copy_file_counting(&src, &dst, &mut stats.h5_copied);
            }
        }

        if mex {
            let src_dir = sample_dir.join(PATH_TO_COUNT_MEX);
            let dst_dir = dest.join(sample_name);
            // Use create_dir (not create_dir_all) to atomically determine whether
            // we created the directory in this run, avoiding a TOCTOU race condition.
            let dst_dir_newly_created = match fs::create_dir(&dst_dir) {
                Ok(()) => {
                    stats.mex_dirs_created += 1;
                    true
                }
                Err(e) if e.kind() == io::ErrorKind::AlreadyExists => false,
                Err(e) => {
                    eprintln!("Failed to create {}: {}", dst_dir.display(), e);
                    continue;
                }
            };
            if src_dir.is_dir() {
                // Copy non-recursively, mirroring the Python behavior
                match fs::read_dir(&src_dir) {
                    Ok(files) => {
                        for f in files.flatten() {
                            let from = f.path();
                            if from.is_file() {
                                let to = dst_dir.join(f.file_name());
                                copy_file_counting(&from, &to, &mut stats.mex_files_copied);
                            }
                        }
                    }
                    Err(e) => {
                        eprintln!("Failed to read {}: {}", src_dir.display(), e);
                        // Only remove the destination directory if it was created in this run
                        cleanup_mex_dir_if_new(&mut stats, &dst_dir, dst_dir_newly_created);
                    }
                }
            } else {
                let src = sample_dir.join(PATH_TO_COUNT_MEX_TAR);
                if src.is_file() {
                    let dst = dest.join(format!("{sample_name}.tar.gz"));
                    copy_file_counting(&src, &dst, &mut stats.mex_files_copied);
                }
                // Source MEX missing; remove created empty dir only if new
                cleanup_mex_dir_if_new(&mut stats, &dst_dir, dst_dir_newly_created);
            }
        }

        if vdj {
            let src = sample_dir.join(PATH_TO_VDJ_ANN);
            if src.is_file() {
                let dst = dest.join(format!("{sample_name}.csv"));
                copy_file_counting(&src, &dst, &mut stats.vdj_copied);
            }
        }
    }
    Ok(stats)
}

fn cleanup_mex_dir_if_new(stats: &mut CopyStats, dst_dir: &PathBuf, dst_dir_newly_created: bool) {
    if dst_dir_newly_created {
        if let Err(e) = fs::remove_dir(dst_dir) {
            eprintln!("Failed to remove {}: {}", dst_dir.display(), e);
        } else {
            stats.mex_dirs_cleaned_up += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::{fixture, rstest};
    use std::fs;
    use tempfile::{TempDir, tempdir};

    #[fixture]
    fn temp_dir() -> TempDir {
        tempdir().unwrap()
    }

    fn sample_dir(temp_dir: TempDir, write_mex: bool) -> (PathBuf, PathBuf, PathBuf, TempDir) {
        let base = temp_dir.path().join("base");
        let dest = temp_dir.path().join("dest");
        let sample = Path::new("p1/outs/per_sample_outs/S1");
        let sample = base.join(sample);
        let count_dir = sample.join("count");

        let mex_src = count_dir.join("sample_filtered_feature_bc_matrix");
        let vdj_dir = sample.join("vdj_t");
        fs::create_dir_all(&count_dir).unwrap();

        fs::create_dir_all(&vdj_dir).unwrap();
        fs::create_dir_all(&dest).unwrap();

        fs::write(
            count_dir.join("sample_filtered_feature_bc_matrix.h5"),
            b"abc",
        )
        .unwrap();

        fs::write(vdj_dir.join("filtered_contig_annotations.csv"), b"contigs").unwrap();

        if write_mex {
            fs::create_dir_all(&mex_src).unwrap();
            fs::write(mex_src.join("barcodes.tsv.gz"), b"barcodes").unwrap();
            fs::write(mex_src.join("features.tsv.gz"), b"features").unwrap();
            fs::write(mex_src.join("matrix.mtx.gz"), b"matrix").unwrap();
            // Nested file should NOT be copied (non-recursive copy behavior)
            let nested = mex_src.join("nested");
            fs::create_dir_all(&nested).unwrap();
            fs::write(nested.join("ignored.txt"), b"ignore me").unwrap();
        }

        (sample, base, dest, temp_dir)
    }

    #[fixture]
    fn sample_dir_mex(temp_dir: TempDir) -> (PathBuf, PathBuf, PathBuf, TempDir) {
        sample_dir(temp_dir, true)
    }

    #[fixture]
    fn sample_dir_no_mex(temp_dir: TempDir) -> (PathBuf, PathBuf, PathBuf, TempDir) {
        sample_dir(temp_dir, false)
    }

    #[rstest]
    #[case(PipestanceResults{ h5: true, mex: false, vdj: false })]
    #[case(PipestanceResults{ h5: false, mex: true, vdj: false })]
    #[case(PipestanceResults{ h5: false, mex: false, vdj: true })]
    #[case(PipestanceResults{ h5: true, mex: true, vdj: true })]
    #[should_panic(expected = "Nothing to do: specify at least one of --h5, --mex, --vdj")]
    #[case(PipestanceResults{ h5: false, mex: false, vdj: false })]
    fn copies_output(
        #[from(sample_dir_mex)] sample_dir: (PathBuf, PathBuf, PathBuf, TempDir),
        #[case] pipestance_results: PipestanceResults,
    ) {
        let (_sample, base, dest, _temp_dir) = sample_dir;
        let stats =
            copy_cellranger_outs_main(base.as_path(), Some(&dest), &pipestance_results, false)
                .unwrap();

        assert_eq!(stats.pipestances_seen, 1);
        assert_eq!(stats.samples_seen, 1);

        if pipestance_results.h5 {
            let out = dest.join("S1.h5");
            assert!(out.exists());
            assert!(out.is_file());
            assert_eq!(fs::read(out).unwrap(), b"abc");
        }

        if pipestance_results.mex {
            let dst_dir = dest.join("S1");
            assert!(dst_dir.is_dir());
            assert!(dst_dir.join("barcodes.tsv.gz").is_file());
            assert!(dst_dir.join("features.tsv.gz").is_file());
            assert!(dst_dir.join("matrix.mtx.gz").is_file());
            assert!(!dst_dir.join("ignored.txt").exists());
            assert!(!dst_dir.join("nested").exists());

            assert_eq!(stats.mex_dirs_created, 1);
            assert_eq!(stats.mex_files_copied, 3);
            assert_eq!(stats.mex_dirs_cleaned_up, 0);
        }
        if pipestance_results.vdj {
            let out = dest.join("S1.csv");
            assert!(out.is_file());
            assert_eq!(fs::read(out).unwrap(), b"contigs");

            assert_eq!(stats.vdj_copied, 1);
        }
    }

    #[rstest]
    fn copies_output_no_mex(temp_dir: TempDir) {
        let base = temp_dir.path().join("base");
        let dest = temp_dir.path().join("dest");

        // Sample exists, but MEX source dir is intentionally missing
        let sample_root = base.join("p1/outs/per_sample_outs/S1");
        fs::create_dir_all(&sample_root).unwrap();

        let stats = copy_cellranger_outs_main(
            &base,
            Some(&dest.clone()),
            &PipestanceResults {
                h5: false,
                mex: true,
                vdj: false,
            },
            false,
        )
        .unwrap();

        // Destination S1 dir should have been created then cleaned up
        assert!(!dest.join("S1").exists());
        assert!(dest.is_dir());

        assert_eq!(stats.pipestances_seen, 1);
        assert_eq!(stats.samples_seen, 1);
        assert_eq!(stats.mex_dirs_created, 1);
        assert_eq!(stats.mex_files_copied, 0);
        assert_eq!(stats.mex_dirs_cleaned_up, 1);
    }

    #[rstest]
    fn does_not_remove_preexisting_mex_dest_dir_when_source_missing(temp_dir: TempDir) {
        let base = temp_dir.path().join("base");
        let dest = temp_dir.path().join("dest");

        let sample_root = base.join("p1/outs/per_sample_outs/S1");
        fs::create_dir_all(&sample_root).unwrap();

        let preexisting = dest.join("S1");
        fs::create_dir_all(&preexisting).unwrap();
        fs::write(preexisting.join("keep.txt"), b"keep").unwrap();

        let stats = copy_cellranger_outs_main(
            &base,
            Some(&dest),
            &PipestanceResults {
                h5: false,
                mex: true,
                vdj: false,
            },
            false,
        )
        .unwrap();

        assert!(preexisting.is_dir());
        assert!(preexisting.join("keep.txt").is_file());

        assert_eq!(stats.mex_dirs_created, 0);
        assert_eq!(stats.mex_dirs_cleaned_up, 0);
    }
    #[rstest]
    fn test_mri_tgz(temp_dir: TempDir) {
        let base = temp_dir.path().join("base");
        fs::create_dir_all(&base).unwrap();
        fs::write(base.join("test.mri.tgz"), b"").expect("Failed to write test.mri.tgz");
        assert!(has_mri_tgz(&base).is_ok())
    }
    #[rstest]
    fn test_handle_command(temp_dir: TempDir) {
        let base = temp_dir.path().join("base");
        let dest = temp_dir.path().join("dest");
        let cmd = Commands::CopyCellRangerOuts {
            base_path: base,
            dest: Some(dest),
            pipestance_results: PipestanceResults {
                h5: false,
                mex: true,
                vdj: false,
            },
            check: false,
        };
        handle_command(cmd);
    }
}
