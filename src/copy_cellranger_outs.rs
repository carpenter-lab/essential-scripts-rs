use clap::Subcommand;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

const PATH_TO_OUTS: &str = "outs/per_sample_outs";
const PATH_TO_COUNT_H5: &str = "count/sample_filtered_feature_bc_matrix.h5";
const PATH_TO_COUNT_MEX: &str = "count/sample_filtered_feature_bc_matrix";
const PATH_TO_VDJ_ANN: &str = "vdj_t/filtered_contig_annotations.csv";

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
        dest: PathBuf,

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

        #[arg(
            long,
            default_value_t = false,
            help = "Check each pipestance for presence of a *.mri.tgz marker and print results. No copying is performed."
        )]
        check: bool,
    },
}

pub fn handle_command(cmd: Commands) -> () {
    match cmd {
        Commands::CopyCellRangerOuts {
            base_path,
            dest,
            h5,
            mex,
            vdj,
            check,
        } => copy_cellranger_outs_main(&base_path, &dest, h5, mex, vdj, check),
    }
}

fn copy_cellranger_outs_main(
    base_path: &Path,
    dest: &Path,
    h5: bool,
    mex: bool,
    vdj: bool,
    check: bool,
) {
    // If in check-only mode, validate pipestances and exit
    if check {
        if let Err(e) = check_pipestances(base_path) {
            eprintln!("{}", e);
        }
        return;
    }

    if !h5 && !mex && !vdj {
        eprintln!("Nothing to do: specify at least one of --h5, --mex, --vdj");
        return;
    }

    // Ensure destination exists
    if let Err(e) = fs::create_dir_all(dest) {
        eprintln!(
            "Failed to create destination directory {}: {}",
            dest.display(),
            e
        );
        return;
    }

    let Ok(entries) = fs::read_dir(base_path) else {
        eprintln!(
            "Base path does not exist or cannot be read: {}",
            base_path.display()
        );
        return;
    };

    for entry in entries.flatten() {
        let dir_path = entry.path();
        if !dir_path.is_dir() {
            continue;
        }

        if let Err(e) = copy_outputs_for_pipestance(&dir_path, dest, h5, mex, vdj) {
            eprintln!("{}: {}", dir_path.display(), e);
        }
    }
}

fn has_mri_tgz(pipestance_dir: &Path) -> io::Result<bool> {
    let mut ok = false;
    for entry in fs::read_dir(pipestance_dir)? {
        let entry = entry?;
        let p = entry.path();
        if p.is_file() {
            if let Some(name) = p.file_name().and_then(|s| s.to_str()) {
                if name.ends_with(".mri.tgz") {
                    ok = true;
                    break;
                }
            }
        }
    }
    Ok(ok)
}

fn check_pipestances(base_path: &Path) -> io::Result<()> {
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
                println!("{}: MISSING .mri.tgz", dir_path.display());
            }
            Err(e) => {
                println!("{}: ERROR reading directory: {}", dir_path.display(), e);
            }
        }
    }
    println!(
        "Checked {} pipestance(s): {} OK, {} missing",
        total,
        ok_count,
        total - ok_count
    );
    Ok(())
}

fn list_samples(pipestance_dir: &Path) -> io::Result<Vec<PathBuf>> {
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

fn copy_outputs_for_pipestance(
    pipestance_dir: &Path,
    dest: &Path,
    h5: bool,
    mex: bool,
    vdj: bool,
) -> io::Result<()> {
    for sample_dir in list_samples(pipestance_dir)? {
        let sample_name = sample_dir
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("sample");

        if h5 {
            let src = sample_dir.join(PATH_TO_COUNT_H5);
            if src.is_file() {
                let dst = dest.join(format!("{}.h5", sample_name));
                let _ = fs::copy(&src, &dst);
            }
        }

        if mex {
            let src_dir = sample_dir.join(PATH_TO_COUNT_MEX);
            let dst_dir = dest.join(sample_name);
            if let Err(e) = fs::create_dir_all(&dst_dir) {
                eprintln!("Failed to create {}: {}", dst_dir.display(), e);
            } else if src_dir.is_dir() {
                // Copy non-recursively, mirroring the Python behavior
                match fs::read_dir(&src_dir) {
                    Ok(files) => {
                        for f in files.flatten() {
                            let from = f.path();
                            if from.is_file() {
                                let to = dst_dir.join(f.file_name());
                                if let Err(e) = fs::copy(&from, &to) {
                                    eprintln!(
                                        "Failed to copy {} to {}: {}",
                                        from.display(),
                                        to.display(),
                                        e
                                    );
                                }
                            }
                        }
                    }
                    Err(_) => {
                        // Remove created directory if we cannot read the source
                        if let Err(e) = fs::remove_dir(&dst_dir) {
                            eprintln!("Failed to remove {}: {}", dst_dir.display(), e);
                        }
                    }
                }
            } else {
                // Source MEX missing; remove created empty dir
                if let Err(e) = fs::remove_dir(&dst_dir) {
                    eprintln!("Failed to remove {}: {}", dst_dir.display(), e);
                }
            }
        }

        if vdj {
            let src = sample_dir.join(PATH_TO_VDJ_ANN);
            if src.is_file() {
                let dst = dest.join(format!("{}.csv", sample_name));
                if let Err(e) = fs::copy(&src, &dst) {
                    eprintln!(
                        "Failed to copy {} to {}: {}",
                        src.display(),
                        dst.display(),
                        e
                    );
                }
            }
        }
    }
    Ok(())
}
