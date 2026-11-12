use crate::cache;
use clap::ValueEnum;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fmt::Display;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
pub const ALPHABET: &str = "ARNDCQEGHILKMFPSTWYV";
pub const P_KMIN: usize = 1;
pub const P_KMAX: usize = 30;
pub const P_BETA: f32 = 0.11387;

pub const BLM_QIJ: [[f32; 20]; 20] = [
    [
        0.0215, 0.0023, 0.0019, 0.0022, 0.0016, 0.0019, 0.003, 0.0058, 0.0011, 0.0032, 0.0044,
        0.0033, 0.0013, 0.0016, 0.0022, 0.0063, 0.0037, 0.0004, 0.0013, 0.0051,
    ],
    [
        0.0023, 0.0178, 0.002, 0.0016, 0.0004, 0.0025, 0.0027, 0.0017, 0.0012, 0.0012, 0.0024,
        0.0062, 0.0008, 0.0009, 0.001, 0.0023, 0.0018, 0.0003, 0.0009, 0.0016,
    ],
    [
        0.0019, 0.002, 0.0141, 0.0037, 0.0004, 0.0015, 0.0022, 0.0029, 0.0014, 0.001, 0.0014,
        0.0024, 0.0005, 0.0008, 0.0009, 0.0031, 0.0022, 0.0002, 0.0007, 0.0012,
    ],
    [
        0.0022, 0.0016, 0.0037, 0.0213, 0.0004, 0.0016, 0.0049, 0.0025, 0.001, 0.0012, 0.0015,
        0.0024, 0.0005, 0.0008, 0.0012, 0.0028, 0.0019, 0.0002, 0.0006, 0.0013,
    ],
    [
        0.0016, 0.0004, 0.0004, 0.0004, 0.0119, 0.0003, 0.0004, 0.0008, 0.0002, 0.0011, 0.0016,
        0.0005, 0.0004, 0.0005, 0.0004, 0.001, 0.0009, 0.0001, 0.0003, 0.0014,
    ],
    [
        0.0019, 0.0025, 0.0015, 0.0016, 0.0003, 0.0073, 0.0035, 0.0014, 0.001, 0.0009, 0.0016,
        0.0031, 0.0007, 0.0005, 0.0008, 0.0019, 0.0014, 0.0002, 0.0007, 0.0012,
    ],
    [
        0.003, 0.0027, 0.0022, 0.0049, 0.0004, 0.0035, 0.0161, 0.0019, 0.0014, 0.0012, 0.002,
        0.0041, 0.0007, 0.0009, 0.0014, 0.003, 0.002, 0.0003, 0.0009, 0.0017,
    ],
    [
        0.0058, 0.0017, 0.0029, 0.0025, 0.0008, 0.0014, 0.0019, 0.0378, 0.001, 0.0014, 0.0021,
        0.0025, 0.0007, 0.0012, 0.0014, 0.0038, 0.0022, 0.0004, 0.0008, 0.0018,
    ],
    [
        0.0011, 0.0012, 0.0014, 0.001, 0.0002, 0.001, 0.0014, 0.001, 0.0093, 0.0006, 0.001, 0.0012,
        0.0004, 0.0008, 0.0005, 0.0011, 0.0007, 0.0002, 0.0015, 0.0006,
    ],
    [
        0.0032, 0.0012, 0.001, 0.0012, 0.0011, 0.0009, 0.0012, 0.0014, 0.0006, 0.0184, 0.0114,
        0.0016, 0.0025, 0.003, 0.001, 0.0017, 0.0027, 0.0004, 0.0014, 0.012,
    ],
    [
        0.0044, 0.0024, 0.0014, 0.0015, 0.0016, 0.0016, 0.002, 0.0021, 0.001, 0.0114, 0.0371,
        0.0025, 0.0049, 0.0054, 0.0014, 0.0024, 0.0033, 0.0007, 0.0022, 0.0095,
    ],
    [
        0.0033, 0.0062, 0.0024, 0.0024, 0.0005, 0.0031, 0.0041, 0.0025, 0.0012, 0.0016, 0.0025,
        0.0161, 0.0009, 0.0009, 0.0016, 0.0031, 0.0023, 0.0003, 0.001, 0.0019,
    ],
    [
        0.0013, 0.0008, 0.0005, 0.0005, 0.0004, 0.0007, 0.0007, 0.0007, 0.0004, 0.0025, 0.0049,
        0.0009, 0.004, 0.0012, 0.0004, 0.0009, 0.001, 0.0002, 0.0006, 0.0023,
    ],
    [
        0.0016, 0.0009, 0.0008, 0.0008, 0.0005, 0.0005, 0.0009, 0.0012, 0.0008, 0.003, 0.0054,
        0.0009, 0.0012, 0.0183, 0.0005, 0.0012, 0.0012, 0.0008, 0.0042, 0.0026,
    ],
    [
        0.0022, 0.001, 0.0009, 0.0012, 0.0004, 0.0008, 0.0014, 0.0014, 0.0005, 0.001, 0.0014,
        0.0016, 0.0004, 0.0005, 0.0191, 0.0017, 0.0014, 0.0001, 0.0005, 0.0012,
    ],
    [
        0.0063, 0.0023, 0.0031, 0.0028, 0.001, 0.0019, 0.003, 0.0038, 0.0011, 0.0017, 0.0024,
        0.0031, 0.0009, 0.0012, 0.0017, 0.0126, 0.0047, 0.0003, 0.001, 0.0024,
    ],
    [
        0.0037, 0.0018, 0.0022, 0.0019, 0.0009, 0.0014, 0.002, 0.0022, 0.0007, 0.0027, 0.0033,
        0.0023, 0.001, 0.0012, 0.0014, 0.0047, 0.0125, 0.0003, 0.0009, 0.0036,
    ],
    [
        0.0004, 0.0003, 0.0002, 0.0002, 0.0001, 0.0002, 0.0003, 0.0004, 0.0002, 0.0004, 0.0007,
        0.0003, 0.0002, 0.0008, 0.0001, 0.0003, 0.0003, 0.0065, 0.0009, 0.0004,
    ],
    [
        0.0013, 0.0009, 0.0007, 0.0006, 0.0003, 0.0007, 0.0009, 0.0008, 0.0015, 0.0014, 0.0022,
        0.001, 0.0006, 0.0042, 0.0005, 0.001, 0.0009, 0.0009, 0.0102, 0.0015,
    ],
    [
        0.0051, 0.0016, 0.0012, 0.0013, 0.0014, 0.0012, 0.0017, 0.0018, 0.0006, 0.012, 0.0095,
        0.0019, 0.0023, 0.0026, 0.0012, 0.0024, 0.0036, 0.0004, 0.0015, 0.0196,
    ],
];

#[derive(Clone, Debug)]
pub struct Peptide {
    pub seq: String,
    pub len: usize,
    pub aff: f32,
    pub indices: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IedbDataRow {
    pub amino_acid_sequence: String,
    pub original_sequence: String,
    pub receptor_group: String,
    pub epitope: String,
    pub organism: String,
    pub antigen: String,
}

pub fn fmatrix_k1() -> [[f32; 20]; 20] {
    let mut marg = [0.0f32; 20];
    for j in 0..20 {
        let mut sum = 0.0f32;
        for item in BLM_QIJ[j] {
            sum += item;
        }
        marg[j] = sum;
    }
    let mut k1 = [[0.0f32; 20]; 20];
    for j in 0..20 {
        for k in 0..20 {
            let val = BLM_QIJ[j][k] / (marg[j] * marg[k]);
            k1[j][k] = val.powf(P_BETA);
        }
    }
    k1
}

pub fn trim(line: &str) -> &str {
    if line.len() >= 2 && line.starts_with('C') && (line.ends_with('F') || line.ends_with('W')) {
        &line[1..line.len() - 1]
    } else {
        line
    }
}

pub fn is_tcr_gene(line: &str) -> bool {
    let mut fields = line.split_whitespace();
    for i in 0..=6 {
        if let Some(field) = fields.next() {
            if (i == 4 || i == 5 || i == 6) && field.starts_with("TRB") {
                return true;
            }
        } else {
            break;
        }
    }
    false
}

pub fn get_nth_field(line: &str, n: usize) -> &str {
    line.split_whitespace().nth(n).unwrap_or("")
}

pub fn read_iedb_data<P: AsRef<Path>>(iedb_data_file: P) -> anyhow::Result<Vec<String>> {
    let file = File::open(iedb_data_file)?;
    let reader = BufReader::new(file);
    let mut iedb_data = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if let Some(sequence) = line.split_whitespace().next()
            && sequence != "trimmed_seq"
        {
            iedb_data.push(sequence.to_string());
        }
    }

    Ok(iedb_data)
}

pub fn create_iedb_map<P: AsRef<Path>>(
    iedb_data_file: &P,
) -> anyhow::Result<HashMap<String, Vec<IedbDataRow>>> {
    let file = File::open(iedb_data_file)?;
    let reader = BufReader::new(file);
    let mut iedb_map: HashMap<String, Vec<IedbDataRow>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let values: Vec<&str> = line.split('\t').collect();
        if values.is_empty() {
            continue;
        }
        let seq = values[0].to_string();
        let row = IedbDataRow {
            amino_acid_sequence: values.first().copied().unwrap_or("").to_string(),
            original_sequence: values.get(1).copied().unwrap_or("").to_string(),
            receptor_group: values.get(2).copied().unwrap_or("").to_string(),
            epitope: values.get(3).copied().unwrap_or("").to_string(),
            organism: values.get(4).copied().unwrap_or("").to_string(),
            antigen: values.get(5).copied().unwrap_or("").to_string(),
        };

        iedb_map.entry(seq).or_default().push(row);
    }

    Ok(iedb_map)
}

pub fn read_airr_data<P: AsRef<Path>>(airr_data_file: P) -> anyhow::Result<Vec<String>> {
    let file = File::open(&airr_data_file)?;
    let mut reader = BufReader::new(file);
    let mut first_line = String::new();
    if reader.read_line(&mut first_line)? == 0 {
        return Ok(Vec::new());
    }

    let mut column = None;
    for (idx, temp) in first_line
        .trim_end_matches(['\r', '\n'])
        .split('\t')
        .enumerate()
    {
        if temp == "cdr3_aa" {
            column = Some(idx);
            break;
        }
    }

    let column_idx = match column {
        Some(c) => c,
        None => {
            eprintln!("Cannot find cdr3_aa column in AIRR TSV input file");
            return Ok(Vec::new());
        }
    };

    let mut airr_data = Vec::new();
    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        let values: Vec<&str> = trimmed.split('\t').collect();
        if let Some(&val) = values.get(column_idx) {
            airr_data.push(val.to_string());
        }
        line.clear();
    }

    Ok(airr_data)
}

#[inline]
pub fn k3_sum(k1: &[[f32; 20]; 20], pep1: &Peptide, pep2: &Peptide) -> f32 {
    let len1 = pep1.len;
    let len2 = pep2.len;
    let max_k = P_KMAX.min(len1).min(len2);
    if len1 == 0 || len2 == 0 || max_k < P_KMIN {
        return 0.0;
    }

    let mut k3 = 0.0f32;
    let buf_len = len1 * len2;

    if buf_len <= 1024 {
        let mut prev = [0.0f32; 1024];
        let mut curr = [0.0f32; 1024];

        for k in P_KMIN..=max_k {
            let limit1 = len1 - k;
            let limit2 = len2 - k;
            let k_offset = k - 1;

            for start1 in 0..=limit1 {
                let j1 = pep1.indices[start1 + k_offset];
                let row_idx = start1 * len2;
                let k1_row = &k1[j1];

                for start2 in 0..=limit2 {
                    let j2 = pep2.indices[start2 + k_offset];
                    let term = k1_row[j2];
                    let idx = row_idx + start2;
                    let k2 = if k == 1 { term } else { prev[idx] * term };
                    curr[idx] = k2;
                    k3 += k2;
                }
            }
            prev[..buf_len].copy_from_slice(&curr[..buf_len]);
        }
    } else {
        let mut prev = vec![0.0f32; buf_len];
        let mut curr = vec![0.0f32; buf_len];

        for k in P_KMIN..=max_k {
            let limit1 = len1 - k;
            let limit2 = len2 - k;
            let k_offset = k - 1;

            for start1 in 0..=limit1 {
                let j1 = pep1.indices[start1 + k_offset];
                let row_idx = start1 * len2;
                let k1_row = &k1[j1];

                for start2 in 0..=limit2 {
                    let j2 = pep2.indices[start2 + k_offset];
                    let term = k1_row[j2];
                    let idx = row_idx + start2;
                    let k2 = if k == 1 { term } else { prev[idx] * term };
                    curr[idx] = k2;
                    k3 += k2;
                }
            }
            std::mem::swap(&mut prev, &mut curr);
        }
    }

    k3
}

pub fn populate_peptide_affinities(k1: &[[f32; 20]; 20], peplist: &mut [Peptide]) {
    peplist.par_iter_mut().for_each(|pep| {
        pep.indices = pep
            .seq
            .bytes()
            .map(|b| ALPHABET.bytes().position(|x| x == b).unwrap())
            .collect();
        pep.aff = k3_sum(k1, pep, pep);
    });
}

pub struct SequenceMatch<'a> {
    pub seq1: &'a str,
    pub seq2: &'a str,
    pub score: f32,
    pub target_index: usize,
}

pub fn find_matches<'a>(
    k1: &[[f32; 20]; 20],
    peplist1: &'a [Peptide],
    peplist2: &'a [Peptide],
    threshold: f32,
    memory_limit: usize,
) -> anyhow::Result<Vec<SequenceMatch<'a>>> {
    let current_memory_usage = AtomicUsize::new(0);
    let approx_item_size = size_of::<SequenceMatch>() + 64;

    let raw_matches: anyhow::Result<Vec<Vec<SequenceMatch<'_>>>> = peplist1
        .par_iter()
        .map(|pep1| {
            peplist2
                .iter()
                .enumerate()
                .filter_map(|(idx2, pep2)| {
                    let score = k3_sum(k1, pep1, pep2) / (pep1.aff * pep2.aff).sqrt();

                    if score <= threshold {
                        return None;
                    }

                    // Check and update memory usage safely
                    let prev_mem = current_memory_usage.fetch_add(approx_item_size, Ordering::Relaxed);
                    if prev_mem + approx_item_size > memory_limit {
                        let limit_gb = memory_limit / (1024 * 1024 * 1024);
                        return Some(Err(anyhow::anyhow!(
                        "TCRMatch exceeded the {limit_gb}GB limit for RAM usage.\n\
                        Consider partitioning your input or increasing the memory limit using the -m flag.\n\
                        Aborting."
                    )));
                    }

                    Some(Ok(SequenceMatch {
                        seq1: pep1.seq.as_str(),
                        seq2: pep2.seq.as_str(),
                        score,
                        target_index: idx2,
                    }))
                })
                .collect::<anyhow::Result<Vec<_>>>()
        })
        .collect();

    let thread_matches = raw_matches?;
    let mut seen = HashSet::new();
    let mut unique_vec = Vec::new();

    for local in thread_matches {
        for item in local {
            let key = (
                item.seq1,
                item.seq2,
                item.score.to_bits(),
                item.target_index,
            );
            if seen.insert(key) {
                unique_vec.push(item);
            }
        }
    }

    Ok(unique_vec)
}

pub fn multi_calc_k3(
    k1: &[[f32; 20]; 20],
    peplist1: &[Peptide],
    peplist2: &[Peptide],
    threshold: f32,
    iedb_map: &HashMap<String, Vec<IedbDataRow>>,
    memory_limit: usize,
) -> anyhow::Result<Vec<String>> {
    let unique_vec = find_matches(k1, peplist1, peplist2, threshold, memory_limit)?;
    let mut tcrmatch_output = Vec::new();

    for tuple in unique_vec {
        let match_peptide = tuple.seq1;
        if let Some(match_row_vec) = iedb_map.get(match_peptide) {
            for match_row in match_row_vec {
                let output_string = format!(
                    "{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\n",
                    tuple.seq1,
                    tuple.seq2,
                    tuple.score,
                    match_row.receptor_group,
                    match_row.epitope,
                    match_row.antigen,
                    match_row.organism
                );
                tcrmatch_output.push(output_string);
            }
        }
    }

    Ok(tcrmatch_output)
}

#[derive(Clone, Debug, ValueEnum)]
pub enum Database {
    Iedb,
    Cedar,
}

impl Display for Database {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Database::Iedb => write!(f, "IEDB"),
            Database::Cedar => write!(f, "CEDAR"),
        }
    }
}

pub trait GetData {
    fn get_data(&self) -> anyhow::Result<PathBuf>;
}

impl GetData for Database {
    fn get_data(&self) -> anyhow::Result<PathBuf> {
        match self {
            Database::Iedb => get_iedb_data(),
            Database::Cedar => get_cedar_data(),
        }
    }
}

fn get_iedb_data() -> anyhow::Result<PathBuf> {
    let url =
        "https://raw.githubusercontent.com/IEDB/TCRMatch/refs/heads/master/data/IEDB_data.tsv";
    cache::cached_file(
        None,
        "IEDB_data.tsv",
        url,
        false,
        false,
        41200,
        "e0fa6002f24eef484be629eb8fc7d893",
    )
}

fn get_cedar_data() -> anyhow::Result<PathBuf> {
    let url =
        "https://raw.githubusercontent.com/IEDB/TCRMatch/refs/heads/master/data/CEDAR_data.tsv";
    cache::cached_file(
        None,
        "CEDAR_data.tsv",
        url,
        false,
        false,
        3792,
        "372ab8e287eb12d3b2cd3039784c52ff",
    )
}

#[allow(clippy::too_many_arguments)]
pub fn run_tcrmatch<W: Write, T: GetData>(
    input_file: &Path,
    database: &T,
    threads: usize,
    threshold: f32,
    airr: bool,
    keep_flanking: bool,
    trust4: bool,
    max_ram_gb: Option<usize>,
    mut writer: W,
) -> anyhow::Result<()> {
    if threads < 1 {
        anyhow::bail!("Number of threads cannot be less than one.");
    }
    if !(0.0..=1.0).contains(&threshold) {
        anyhow::bail!("Threshold must be between 0 and 1");
    }

    let memory_limit = max_ram_gb.unwrap_or(2) * 1024 * 1024 * 1024;
    let trimming = !keep_flanking;
    let database = database.get_data()?;

    let iedb_map = create_iedb_map(&database)?;
    let iedb_data = read_iedb_data(database)?;

    if iedb_data.is_empty() {
        eprintln!(
            "Failed to read database file - please specify the path to database file if not ./data/IEDB_data.tsv"
        );
        anyhow::bail!("Failed to read database file");
    }

    let k1 = fmatrix_k1();
    let mut peplist1 = Vec::new();
    let mut inputlines = Vec::new();
    let mut trust4header = String::new();

    if airr {
        let airr_data = read_airr_data(input_file)?;
        if airr_data.is_empty() {
            return Ok(());
        }
        for seq in airr_data {
            let mut invalid_residue = false;
            for (i, b) in seq.bytes().enumerate() {
                if !ALPHABET.as_bytes().contains(&b) {
                    eprintln!(
                        "Invalid amino acid found in {} at position {}. Skipping.",
                        seq,
                        i + 1
                    );
                    invalid_residue = true;
                    break;
                }
            }
            if invalid_residue {
                continue;
            }
            let len = seq.len();
            peplist1.push(Peptide {
                seq,
                len,
                aff: -99.9,
                indices: Vec::new(),
            });
        }
    } else if trust4 {
        let file = File::open(input_file)?;
        let mut reader = BufReader::new(file);
        let mut header_line = String::new();
        reader.read_line(&mut header_line)?;
        trust4header = header_line.trim_end_matches(['\r', '\n']).to_string();

        let valid_header = "#count\tfrequency\tCDR3nt\tCDR3aa\tV\tD\tJ\tC\tcid\tcid_full_length";
        if trust4header != valid_header {
            eprintln!("ERROR: Invalid TRUST4 header. Please check input file format.");
            eprintln!("Input file: {}", input_file.display());
            anyhow::bail!("Invalid TRUST4 header");
        }

        let mut line = String::new();
        while reader.read_line(&mut line)? > 0 {
            let trimmed_line = line.trim_end_matches(['\r', '\n']);
            if is_tcr_gene(trimmed_line) {
                inputlines.push(trimmed_line.to_string());

                let seq_raw = get_nth_field(trimmed_line, 3);
                let mut invalid_char_found = false;
                for b in seq_raw.bytes() {
                    if !ALPHABET.as_bytes().contains(&b) {
                        invalid_char_found = true;
                        break;
                    }
                }

                if invalid_char_found {
                    eprintln!("Invalid amino acid found in {seq_raw}");
                    eprintln!("Skipping");
                    line.clear();
                    continue;
                }

                let seq = if trimming {
                    trim(seq_raw).to_string()
                } else {
                    seq_raw.to_string()
                };

                let len = seq.len();
                peplist1.push(Peptide {
                    seq,
                    len,
                    aff: -99.9,
                    indices: Vec::new(),
                });
            }
            line.clear();
        }
    } else {
        let file = File::open(input_file)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line?;
            let trimmed_line = line.trim_end_matches(['\r', '\n']);
            let mut invalid_residue = false;
            for (i, b) in trimmed_line.bytes().enumerate() {
                if !ALPHABET.as_bytes().contains(&b) {
                    eprintln!(
                        "Invalid amino acid found in {} at position {}. Skipping.",
                        trimmed_line,
                        i + 1
                    );
                    invalid_residue = true;
                    break;
                }
            }
            if invalid_residue {
                continue;
            }

            let seq = if trimming {
                trim(trimmed_line).to_string()
            } else {
                trimmed_line.to_string()
            };

            let len = seq.len();
            peplist1.push(Peptide {
                seq,
                len,
                aff: -99.9,
                indices: Vec::new(),
            });
        }
    }

    let mut peplist2 = Vec::with_capacity(iedb_data.len());
    for seq in iedb_data {
        for (i, b) in seq.bytes().enumerate() {
            if !ALPHABET.as_bytes().contains(&b) {
                eprintln!("Invalid amino acid found in {} at position {}", seq, i + 1);
                anyhow::bail!("Invalid amino acid in database");
            }
        }
        let len = seq.len();
        peplist2.push(Peptide {
            seq,
            len,
            aff: -99.9,
            indices: Vec::new(),
        });
    }

    {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()?;
        pool.install(|| {
            populate_peptide_affinities(&k1, &mut peplist1);
            populate_peptide_affinities(&k1, &mut peplist2);
        });
    }

    let tcrmatch_output = {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()?;
        pool.install(|| {
            multi_calc_k3(
                &k1,
                &peplist1,
                &peplist2,
                threshold,
                &iedb_map,
                memory_limit,
            )
        })?
    };

    let was_trimmed = if trimming { "trimmed_" } else { "" };
    let tcrmatch_output_header = format!(
        "{}input_sequence\tmatch_sequence\tscore\treceptor_group\tepitope\tantigen\torganism",
        was_trimmed
    );

    if trust4 {
        writeln!(writer, "{trust4header}\t{tcrmatch_output_header}")?;
        for input_line in &inputlines {
            let aa_seq_trust = get_nth_field(input_line, 3);
            let trimmed_trust = trim(aa_seq_trust);
            for out in &tcrmatch_output {
                let seq_tcrmatch = get_nth_field(out, 0);
                if trimmed_trust == seq_tcrmatch {
                    write!(writer, "{input_line}\t{out}")?;
                }
            }
        }
    } else {
        writeln!(writer, "{tcrmatch_output_header}")?;
        for out in &tcrmatch_output {
            write!(writer, "{out}")?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trim() {
        assert_eq!(trim("CASSLGF"), "ASSLG");
        assert_eq!(trim("CASW"), "AS");
        assert_eq!(trim("ASSLGF"), "ASSLGF");
        assert_eq!(trim("CASSLGA"), "CASSLGA");
        assert_eq!(trim("C"), "C");
        assert_eq!(trim(""), "");
    }

    #[test]
    fn test_is_tcr_gene() {
        let line1 = "10 0.05 AAAG CASSLGF TRBV2*01 TRBD1*01 TRBJ1*01 TRBC1 cid1 full";
        assert!(is_tcr_gene(line1));

        let line2 = "10 0.05 AAAG CASSLGF IGHV1*01 IGHD1*01 IGHJ1*01 IGHC1 cid1 full";
        assert!(!is_tcr_gene(line2));
    }

    #[test]
    fn test_fmatrix_k1_properties() {
        let k1 = fmatrix_k1();
        for i in 0..20 {
            for j in 0..20 {
                assert!(k1[i][j] > 0.0, "k1[{}][{}] was not positive", i, j);
            }
        }
    }

    #[test]
    fn test_k3_sum_and_score() {
        let k1 = fmatrix_k1();
        let mut pep1 = Peptide {
            seq: "CASSLGF".to_string(),
            len: 7,
            aff: -99.9,
            indices: Vec::new(),
        };
        let mut pep2 = Peptide {
            seq: "CASSLGF".to_string(),
            len: 7,
            aff: -99.9,
            indices: Vec::new(),
        };
        let mut list = vec![pep1.clone(), pep2.clone()];
        populate_peptide_affinities(&k1, &mut list);
        pep1 = list[0].clone();
        pep2 = list[1].clone();

        assert!(pep1.aff > 0.0);
        let score = k3_sum(&k1, &pep1, &pep2) / (pep1.aff * pep2.aff).sqrt();
        assert!(
            (score - 1.0).abs() < 1e-4,
            "Self score should be 1.0, got {}",
            score
        );
    }

    struct TestDatabase {
        data: String,
        path: PathBuf,
    }

    impl GetData for TestDatabase {
        fn get_data(&self) -> anyhow::Result<PathBuf> {
            std::fs::write(&self.path, &self.data)?;
            Ok(PathBuf::from(&self.path))
        }
    }

    #[test]
    fn test_run_tcrmatch_plain_text() {
        let temp_dir = tempfile::tempdir().unwrap();
        let input_file = temp_dir.path().join("input.txt");

        let test_db = TestDatabase {
            data: String::from(
                "trimmed_seq\toriginal_seq\treceptor_group\tepitope\torganism\tantigen\n\
             ASSLG\tCASSLGF\tgroup1\tepi1\torg1\tant1\n",
            ),
            path: temp_dir.path().join("db.tsv"),
        };

        std::fs::write(&input_file, "CASSLGF\n").unwrap();

        let mut output = Vec::new();
        run_tcrmatch(
            &input_file,
            &test_db,
            1,
            0.9,
            false,
            false,
            false,
            None,
            &mut output,
        )
        .unwrap();

        let out_str = String::from_utf8(output).unwrap();
        assert!(out_str.contains("trimmed_input_sequence\tmatch_sequence\tscore\treceptor_group\tepitope\tantigen\torganism"));
        assert!(out_str.contains("ASSLG\tASSLG\t1.000000\tgroup1\tepi1\tant1\torg1"));
    }

    #[test]
    fn test_run_tcrmatch_keep_flanking() {
        let temp_dir = tempfile::tempdir().unwrap();
        let input_file = temp_dir.path().join("input.txt");

        let test_db = TestDatabase {
            data: String::from(
                "trimmed_seq\toriginal_seq\treceptor_group\tepitope\torganism\tantigen\n\
             CASSLGF\tCASSLGF\tgroup1\tepi1\torg1\tant1\n",
            ),
            path: temp_dir.path().join("db.tsv"),
        };

        std::fs::write(&input_file, "CASSLGF\n").unwrap();

        let mut output = Vec::new();
        run_tcrmatch(
            &input_file,
            &test_db,
            1,
            0.9,
            false,
            true,
            false,
            None,
            &mut output,
        )
        .unwrap();

        let out_str = String::from_utf8(output).unwrap();
        assert!(out_str.contains(
            "input_sequence\tmatch_sequence\tscore\treceptor_group\tepitope\tantigen\torganism"
        ));
        assert!(out_str.contains("CASSLGF\tCASSLGF\t1.000000\tgroup1\tepi1\tant1\torg1"));
    }

    #[test]
    fn test_run_tcrmatch_airr() {
        let temp_dir = tempfile::tempdir().unwrap();
        let db_file = temp_dir.path().join("db.tsv");
        let input_file = temp_dir.path().join("input.tsv");

        let test_db = TestDatabase {
            data: String::from(
                "trimmed_seq\toriginal_seq\treceptor_group\tepitope\torganism\tantigen\n\
             CASSLGF\tCASSLGF\tgroup1\tepi1\torg1\tant1\n",
            ),
            path: db_file,
        };

        std::fs::write(
            &input_file,
            "sequence_id\tcdr3_aa\tv_call\n\
             seq1\tCASSLGF\tTRBV2\n",
        )
        .unwrap();

        let mut output = Vec::new();
        run_tcrmatch(
            &input_file,
            &test_db,
            1,
            0.9,
            true,
            true,
            false,
            None,
            &mut output,
        )
        .unwrap();

        let out_str = String::from_utf8(output).unwrap();
        assert!(out_str.contains("CASSLGF\tCASSLGF\t1.000000\tgroup1\tepi1\tant1\torg1"));
    }

    #[test]
    fn test_run_tcrmatch_trust4() {
        let temp_dir = tempfile::tempdir().unwrap();
        let db_file = temp_dir.path().join("db.tsv");
        let input_file = temp_dir.path().join("input.txt");

        let test_db = TestDatabase {
            data: String::from(
                "trimmed_seq\toriginal_seq\treceptor_group\tepitope\torganism\tantigen\n\
             ASSLG\tCASSLGF\tgroup1\tepi1\torg1\tant1\n",
            ),
            path: db_file,
        };

        std::fs::write(
            &input_file,
            "#count\tfrequency\tCDR3nt\tCDR3aa\tV\tD\tJ\tC\tcid\tcid_full_length\n\
             10\t0.05\tAAAG\tCASSLGF\tTRBV2*01\tTRBD1*01\tTRBJ1*01\tTRBC1\tcid1\tfull\n",
        )
        .unwrap();

        let mut output = Vec::new();
        run_tcrmatch(
            &input_file,
            &test_db,
            1,
            0.9,
            false,
            false,
            true,
            None,
            &mut output,
        )
        .unwrap();

        let out_str = String::from_utf8(output).unwrap();
        assert!(out_str.contains("#count\tfrequency\tCDR3nt\tCDR3aa\tV\tD\tJ\tC\tcid\tcid_full_length\ttrimmed_input_sequence\tmatch_sequence\tscore\treceptor_group\tepitope\tantigen\torganism"));
        assert!(out_str.contains("10\t0.05\tAAAG\tCASSLGF\tTRBV2*01\tTRBD1*01\tTRBJ1*01\tTRBC1\tcid1\tfull\tASSLG\tASSLG\t1.000000\tgroup1\tepi1\tant1\torg1"));
    }

    #[test]
    fn test_run_tcrmatch_invalid_params() {
        let temp_dir = tempfile::tempdir().unwrap();
        let db_file = Database::Iedb;
        let input_file = temp_dir.path().join("input.txt");

        let mut output = Vec::new();
        assert!(
            run_tcrmatch(
                &input_file,
                &db_file,
                0,
                0.9,
                false,
                false,
                false,
                None,
                &mut output,
            )
            .is_err()
        );

        assert!(
            run_tcrmatch(
                &input_file,
                &db_file,
                1,
                1.5,
                false,
                false,
                false,
                None,
                &mut output,
            )
            .is_err()
        );
    }
}
