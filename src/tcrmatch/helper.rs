use crate::cache;
use crate::progress::{Progress, ProgressSettings, SetDrawTarget, progress_bar_with_message};
use clap::ValueEnum;
use indicatif::ProgressBar;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fmt::Display;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
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

/// Margin below the threshold within which a pair is still scored exactly.
///
/// [`QueryProfile::k3_sum_fast`] adds terms in a different order from [`k3_sum`], so the two
/// round differently. This margin is far larger than the combined `f32`
/// rounding error of both, so the prefilter never drops a pair that `k3_sum`
/// would accept.
pub const PREFILTER_MARGIN: f64 = 1e-2;

// `k3_sum_fast` sums runs of every length from 1, so it assumes P_KMIN == 1.
const _: () = assert!(P_KMIN == 1);

/// Precomputed `k1` values for one input peptide, used by the fast prefilter.
///
/// `values[c * len + i]` holds `k1[a_i][c]` for input residue `a_i` and database
/// residue `c`. Walking the grid one database residue at a time then reads one
/// contiguous row per step, with no indexed lookups inside the inner loop. The
/// profile is built once per input and reused against the whole database.
///
/// It also owns the scratch rows for the recurrence and the per-column sums,
/// so scoring allocates nothing. Use one profile per thread.
pub struct QueryProfile {
    len: usize,
    values: Vec<f32>,
    prev: Vec<f32>,
    curr: Vec<f32>,
    /// `col_sums[i]` adds up every run ending in input column `i`. One sum per
    /// column, rather than a single running total, means no addition waits on
    /// the previous one, so the inner loop can be vectorized.
    col_sums: Vec<f32>,
}

impl QueryProfile {
    pub fn new(k1: &[[f32; 20]; 20], pep: &Peptide) -> Self {
        let len = pep.indices.len();
        let mut values = vec![0.0f32; 20 * len];
        for (i, &ai) in pep.indices.iter().enumerate() {
            for (c, &v) in k1[ai].iter().enumerate() {
                values[c * len + i] = v;
            }
        }
        Self {
            len,
            values,
            prev: vec![0.0; len + 1],
            curr: vec![0.0; len + 1],
            col_sums: vec![0.0; len],
        }
    }

    /// Approximates [`k3_sum`] of this input against `pep2` in
    /// O(len1 * len2) instead of O(len1 * len2 * k).
    ///
    /// `k3_sum` adds, for every diagonal of the `len1 x len2` grid, the product
    /// of `k1` terms over every contiguous run of length `P_KMIN..=max_k`.
    /// Summing the runs that end at each cell gives the same total in one pass:
    /// `s[j][i] = k1[a_i][b_j] * (1 + s[j-1][i-1])`.
    ///
    /// Works in `f32` for speed and adds terms in a different order from
    /// `k3_sum`, so the two differ by rounding (well under 1e-3 relative,
    /// against a [`PREFILTER_MARGIN`] of 1e-2). Use it to discard pairs, never
    /// to report a score.
    ///
    /// Returns `None` when this shortcut does not apply: an empty sequence, or
    /// both sequences longer than `P_KMAX`, where `k3_sum` caps the run length.
    pub fn k3_sum_fast(&mut self, pep2: &Peptide) -> Option<f64> {
        self.k3_sum_fast_residues(&pep2.indices)
    }

    /// [`Self::k3_sum_fast`] on raw residue codes (indices into [`ALPHABET`]),
    /// so callers can pass a packed `u8` database without building a `Peptide`.
    pub fn k3_sum_fast_residues<T: Copy>(&mut self, b: &[T]) -> Option<f64>
    where
        usize: From<T>,
    {
        let n = self.len;
        if n == 0 || b.is_empty() || n.min(b.len()) > P_KMAX {
            return None;
        }

        // curr[0] is reset every row and curr[1..] is fully overwritten, so
        // only prev and the column sums need clearing between pairs.
        self.prev.fill(0.0);
        self.col_sums.fill(0.0);
        let (mut prev, mut curr) = (&mut self.prev, &mut self.curr);
        for &bj in b {
            let bj = usize::from(bj);
            let row = &self.values[bj * n..(bj + 1) * n];
            curr[0] = 0.0;
            update_row(&mut curr[1..], &mut self.col_sums, row, &prev[..n]);
            std::mem::swap(&mut prev, &mut curr);
        }
        Some(f64::from(self.col_sums.iter().sum::<f32>()))
    }
}

/// One row of the [`QueryProfile::k3_sum_fast`] recurrence:
/// `curr[i] = row[i] * (1 + prev[i])` and `sums[i] += curr[i]`.
///
/// Takes separate slices, resliced to one length, so the compiler knows the
/// outputs cannot overlap the inputs and needs no per-row length checks.
///
/// Note: as of Rust 1.99 / LLVM this still compiles to scalar code (one `f32`
/// per instruction), with or without `-C target-cpu=native`, even when written
/// in fixed blocks of 4. Vectorizing it would take explicit SIMD.
#[inline(always)]
fn update_row(curr: &mut [f32], sums: &mut [f32], row: &[f32], prev: &[f32]) {
    let n = curr.len();
    let (sums, row, prev) = (&mut sums[..n], &row[..n], &prev[..n]);
    for (((c, s), &k), &p) in curr.iter_mut().zip(sums).zip(row).zip(prev) {
        let run_sum = k * (1.0 + p);
        *c = run_sum;
        *s += run_sum;
    }
}

/// Approximates [`k3_sum`]; see [`QueryProfile::k3_sum_fast`]. Builds a
/// profile on every call, so loops over many pairs should reuse a
/// [`QueryProfile`] instead.
pub fn k3_sum_fast(k1: &[[f32; 20]; 20], pep1: &Peptide, pep2: &Peptide) -> Option<f64> {
    QueryProfile::new(k1, pep1).k3_sum_fast(pep2)
}

/// Residue codes of many peptides packed into one contiguous `u8` buffer.
///
/// Each `Peptide` keeps its codes in its own `Vec<usize>`: 8 bytes per residue
/// in a separate allocation. The prefilter reads the whole database once per
/// input, so packing it (1 byte per residue, one block) keeps that pass in cache.
pub struct PackedResidues {
    codes: Vec<u8>,
    /// Peptide `i` occupies `codes[offsets[i]..offsets[i + 1]]`.
    offsets: Vec<usize>,
}

impl PackedResidues {
    pub fn new(peps: &[Peptide]) -> Self {
        let mut codes = Vec::with_capacity(peps.iter().map(|pep| pep.indices.len()).sum());
        let mut offsets = Vec::with_capacity(peps.len() + 1);
        offsets.push(0);
        for pep in peps {
            // Codes index ALPHABET, so they are always below 20.
            codes.extend(pep.indices.iter().map(|&code| code as u8));
            offsets.push(codes.len());
        }
        Self { codes, offsets }
    }

    #[inline]
    pub fn get(&self, i: usize) -> &[u8] {
        &self.codes[self.offsets[i]..self.offsets[i + 1]]
    }
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
    progress: &ProgressBar,
) -> anyhow::Result<Vec<SequenceMatch<'a>>> {
    let current_memory_usage = AtomicUsize::new(0);
    let approx_item_size = size_of::<SequenceMatch>() + 64;

    // The prefilter normalizes every pair; precompute 1/sqrt(aff) once per
    // peptide so the hot loop only multiplies. The exact score below keeps
    // the original formula so reported scores round exactly as before.
    let inv_sqrt_aff = |pep: &Peptide| 1.0 / f64::from(pep.aff).sqrt();
    let db_inv_sqrt_aff: Vec<f64> = peplist2.iter().map(inv_sqrt_aff).collect();
    let prefilter_cutoff = f64::from(threshold) - PREFILTER_MARGIN;
    let db_residues = PackedResidues::new(peplist2);

    let raw_matches: anyhow::Result<Vec<Vec<SequenceMatch<'_>>>> = peplist1
        .par_iter()
        .map(|pep1| {
            let pep1_inv_sqrt_aff = inv_sqrt_aff(pep1);
            // Built once per input, then reused against every database sequence.
            let mut profile = QueryProfile::new(k1, pep1);
            let pep1_matches = peplist2
                .iter()
                .zip(&db_inv_sqrt_aff)
                .enumerate()
                .filter_map(|(idx2, (pep2, &pep2_inv_sqrt_aff))| {
                    // Cheap prefilter: skip pairs clearly below the threshold.
                    // Anything near or above it is rescored exactly below.
                    if let Some(fast) = profile.k3_sum_fast_residues(db_residues.get(idx2)) {
                        let approx = fast * pep1_inv_sqrt_aff * pep2_inv_sqrt_aff;
                        if approx <= prefilter_cutoff {
                            return None;
                        }
                    }

                    let score = k3_sum(k1, pep1, pep2) / (pep1.aff * pep2.aff).sqrt();

                    // NaN (an empty sequence) must never count as a match.
                    if score.is_nan() || score <= threshold {
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
                .collect::<anyhow::Result<Vec<_>>>();

            // Advance by residues, not sequences: an input's cost grows with its
            // length, so this keeps the rate and ETA steady in any input order.
            progress.inc(pep1.len as u64);
            pep1_matches
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

/// Writes one output row for each IEDB annotation of the match's database
/// sequence. When `prefix` is given, it and a tab start every row.
pub fn write_match_rows<W: Write>(
    writer: &mut W,
    prefix: Option<&str>,
    m: &SequenceMatch<'_>,
    iedb_map: &HashMap<String, Vec<IedbDataRow>>,
) -> std::io::Result<()> {
    // Annotations belong to the database sequence (seq2), not the input (seq1).
    let Some(rows) = iedb_map.get(m.seq2) else {
        return Ok(());
    };
    for row in rows {
        if let Some(prefix) = prefix {
            write!(writer, "{prefix}\t")?;
        }
        writeln!(
            writer,
            "{}\t{}\t{:.6}\t{}\t{}\t{}\t{}",
            m.seq1, m.seq2, m.score, row.receptor_group, row.epitope, row.antigen, row.organism
        )?;
    }
    Ok(())
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
    progress: Progress,
    writer: W,
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

    // Empty sequences (blank lines, empty AIRR cdr3_aa, or CDR3s such as "CF"
    // that trimming empties) score 0/0 = NaN against everything.
    let n_before = peplist1.len();
    peplist1.retain(|pep| !pep.seq.is_empty());
    let n_empty = n_before - peplist1.len();
    if n_empty > 0 {
        eprintln!("Skipping {n_empty} empty input sequence(s).");
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

    // The bar draws to stderr, so it never mixes with results on stdout.
    // It counts input residues; find_matches advances it by each input's length.
    let total_residues: usize = peplist1.iter().map(|pep| pep.len).sum();
    let progress_bar = progress_bar_with_message(
        total_residues,
        "Matching input residues to database",
        ProgressSettings::PositionSpeed,
    )?;
    progress.apply(&progress_bar);

    let matches = {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()?;
        pool.install(|| {
            find_matches(
                &k1,
                &peplist1,
                &peplist2,
                threshold,
                memory_limit,
                &progress_bar,
            )
        })?
    };
    progress_bar.finish();

    let was_trimmed = if trimming { "trimmed_" } else { "" };
    let tcrmatch_output_header = format!(
        "{}input_sequence\tmatch_sequence\tscore\treceptor_group\tepitope\tantigen\torganism",
        was_trimmed
    );

    // A locked stdout flushes on every newline; buffer the rows instead.
    let mut writer = BufWriter::new(writer);

    if trust4 {
        writeln!(writer, "{trust4header}\t{tcrmatch_output_header}")?;

        // Group matches by input sequence, keeping their order within each group.
        let mut matches_by_input: HashMap<&str, Vec<&SequenceMatch<'_>>> = HashMap::new();
        for m in &matches {
            matches_by_input.entry(m.seq1).or_default().push(m);
        }

        for input_line in &inputlines {
            let trimmed_trust = trim(get_nth_field(input_line, 3));
            for m in matches_by_input.get(trimmed_trust).into_iter().flatten() {
                write_match_rows(&mut writer, Some(input_line), m, &iedb_map)?;
            }
        }
    } else {
        writeln!(writer, "{tcrmatch_output_header}")?;
        for m in &matches {
            write_match_rows(&mut writer, None, m, &iedb_map)?;
        }
    }

    writer.flush()?;
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

    #[test]
    fn test_write_match_rows_joins_on_database_sequence() {
        let k1 = fmatrix_k1();
        let pep = |s: &str| Peptide {
            seq: s.to_string(),
            len: s.len(),
            aff: -99.9,
            indices: Vec::new(),
        };
        // The input sequence is not in the database, so annotations must come
        // from the matched database sequence.
        let mut input = vec![pep("ASSLGQ")];
        let mut db = vec![pep("ASSLGE")];
        populate_peptide_affinities(&k1, &mut input);
        populate_peptide_affinities(&k1, &mut db);

        let mut iedb_map = HashMap::new();
        iedb_map.insert(
            "ASSLGE".to_string(),
            vec![IedbDataRow {
                amino_acid_sequence: "ASSLGE".to_string(),
                original_sequence: "CASSLGEF".to_string(),
                receptor_group: "group1".to_string(),
                epitope: "epi1".to_string(),
                organism: "org1".to_string(),
                antigen: "ant1".to_string(),
            }],
        );

        let matches = find_matches(&k1, &input, &db, 0.0, usize::MAX, &ProgressBar::hidden()).unwrap();
        assert_eq!(matches.len(), 1);

        let mut buf = Vec::new();
        write_match_rows(&mut buf, None, &matches[0], &iedb_map).unwrap();
        let out = String::from_utf8(buf).unwrap();

        let rows: Vec<&str> = out.lines().collect();
        assert_eq!(rows.len(), 1, "expected one annotated row, got {out:?}");
        let fields: Vec<&str> = rows[0].split('\t').collect();
        assert_eq!(fields[0], "ASSLGQ");
        assert_eq!(fields[1], "ASSLGE");
        assert_eq!(&fields[3..], ["group1", "epi1", "ant1", "org1"]);

        let mut buf = Vec::new();
        write_match_rows(&mut buf, Some("line\tdata"), &matches[0], &iedb_map).unwrap();
        let out = String::from_utf8(buf).unwrap();
        assert!(out.starts_with("line\tdata\tASSLGQ\tASSLGE\t"), "got {out:?}");
        assert!(out.ends_with('\n'));
    }

    /// Deterministic random peptides (xorshift), so tests need no RNG crate.
    fn random_peptides(seed: u64, count: usize, min_len: usize, max_len: usize) -> Vec<Peptide> {
        let mut state = seed.max(1);
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let alphabet = ALPHABET.as_bytes();
        let mut peps: Vec<Peptide> = (0..count)
            .map(|_| {
                let len = min_len + (next() as usize) % (max_len - min_len + 1);
                let seq: String = (0..len)
                    .map(|_| alphabet[(next() as usize) % alphabet.len()] as char)
                    .collect();
                Peptide {
                    len: seq.len(),
                    seq,
                    aff: -99.9,
                    indices: Vec::new(),
                }
            })
            .collect();
        populate_peptide_affinities(&fmatrix_k1(), &mut peps);
        peps
    }

    #[test]
    fn test_k3_sum_fast_agrees_with_k3_sum() {
        let k1 = fmatrix_k1();
        let a = random_peptides(1, 60, 1, 30);
        let b = random_peptides(2, 60, 1, 45);
        let packed = PackedResidues::new(&b);
        let (mut max_rel, mut max_norm_err) = (0.0f64, 0.0f64);
        for p1 in &a {
            // Reuse one profile, as find_matches does, to check the scratch
            // rows are reset between pairs.
            let mut profile = QueryProfile::new(&k1, p1);
            for (j, p2) in b.iter().enumerate() {
                let exact = f64::from(k3_sum(&k1, p1, p2));
                let fast = profile.k3_sum_fast(p2).expect("shortcut applies");
                let fresh = k3_sum_fast(&k1, p1, p2).expect("shortcut applies");
                assert_eq!(fast.to_bits(), fresh.to_bits(), "reused profile drifted");
                let packed_fast = profile
                    .k3_sum_fast_residues(packed.get(j))
                    .expect("shortcut applies");
                assert_eq!(packed_fast.to_bits(), fast.to_bits(), "packed codes differ");
                let rel = (fast - exact).abs() / exact.abs().max(1e-30);
                assert!(rel < 1e-3, "{} vs {}: fast {fast}, exact {exact}", p1.seq, p2.seq);
                max_rel = max_rel.max(rel);

                let norm = |x: f64| x / (f64::from(p1.aff) * f64::from(p2.aff)).sqrt();
                let norm_err = (norm(fast) - norm(exact)).abs();
                assert!(
                    norm_err < PREFILTER_MARGIN / 10.0,
                    "normalized score error too close to the prefilter margin"
                );
                max_norm_err = max_norm_err.max(norm_err);
            }
        }
        eprintln!("max relative error {max_rel:e}, max normalized error {max_norm_err:e}");
    }

    #[test]
    fn test_k3_sum_fast_falls_back_when_shortcut_does_not_apply() {
        let k1 = fmatrix_k1();
        let empty = random_peptides(3, 1, 0, 0).remove(0);
        let short = random_peptides(4, 1, 10, 10).remove(0);
        let long_a = random_peptides(5, 1, P_KMAX + 1, P_KMAX + 5).remove(0);
        let long_b = random_peptides(6, 1, P_KMAX + 1, P_KMAX + 5).remove(0);

        assert!(k3_sum_fast(&k1, &empty, &short).is_none());
        assert!(k3_sum_fast(&k1, &short, &empty).is_none());
        assert!(k3_sum_fast(&k1, &long_a, &long_b).is_none());
        assert!(k3_sum_fast(&k1, &short, &long_a).is_some());
    }

    #[test]
    fn test_find_matches_identical_to_exact_scan() {
        let k1 = fmatrix_k1();
        // Mutate copies of the database so many pairs land near each threshold.
        let db = random_peptides(7, 80, 8, 18);
        let mut input = random_peptides(8, 40, 8, 18);
        for (i, pep) in input.iter_mut().enumerate().take(30) {
            let mut seq = db[i].seq.clone().into_bytes();
            let pos = i % seq.len();
            seq[pos] = ALPHABET.as_bytes()[(i * 7) % 20];
            pep.seq = String::from_utf8(seq).unwrap();
            pep.len = pep.seq.len();
        }
        populate_peptide_affinities(&k1, &mut input);

        for threshold in [0.5f32, 0.8, 0.9, 0.97] {
            let mut expected = Vec::new();
            for p1 in &input {
                for (j, p2) in db.iter().enumerate() {
                    let score = k3_sum(&k1, p1, p2) / (p1.aff * p2.aff).sqrt();
                    if score > threshold {
                        expected.push((p1.seq.clone(), j, score.to_bits()));
                    }
                }
            }

            let mut got: Vec<_> = find_matches(&k1, &input, &db, threshold, usize::MAX, &ProgressBar::hidden())
                .unwrap()
                .into_iter()
                .map(|m| (m.seq1.to_string(), m.target_index, m.score.to_bits()))
                .collect();
            got.sort();
            expected.sort();
            expected.dedup();

            assert!(!expected.is_empty(), "test data has no matches at {threshold}");
            assert_eq!(got, expected, "mismatch at threshold {threshold}");
        }
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
            Progress::NoProgress,
            &mut output,
        )
        .unwrap();

        let out_str = String::from_utf8(output).unwrap();
        assert!(out_str.contains("trimmed_input_sequence\tmatch_sequence\tscore\treceptor_group\tepitope\tantigen\torganism"));
        assert!(out_str.contains("ASSLG\tASSLG\t1.000000\tgroup1\tepi1\tant1\torg1"));
    }

    #[test]
    fn test_packed_residues_round_trip() {
        let peps = random_peptides(13, 30, 0, 25);
        let packed = PackedResidues::new(&peps);
        for (i, pep) in peps.iter().enumerate() {
            let codes: Vec<usize> = packed.get(i).iter().map(|&c| usize::from(c)).collect();
            assert_eq!(codes, pep.indices, "peptide {i} ({})", pep.seq);
        }
    }

    #[test]
    fn test_find_matches_advances_progress_by_residues() {
        let k1 = fmatrix_k1();
        let input = random_peptides(11, 25, 5, 20);
        let db = random_peptides(12, 10, 8, 15);
        let total_residues: u64 = input.iter().map(|pep| pep.len as u64).sum();

        let pb = ProgressBar::hidden();
        find_matches(&k1, &input, &db, 0.97, usize::MAX, &pb).unwrap();
        assert_eq!(pb.position(), total_residues);
    }

    #[test]
    fn test_find_matches_rejects_empty_sequence() {
        let k1 = fmatrix_k1();
        let input = random_peptides(9, 1, 0, 0);
        let db = random_peptides(10, 5, 8, 15);

        let matches = find_matches(&k1, &input, &db, 0.0, usize::MAX, &ProgressBar::hidden()).unwrap();
        assert!(matches.is_empty(), "empty sequence matched {} entries", matches.len());
    }

    #[test]
    fn test_run_tcrmatch_skips_empty_input_sequences() {
        let temp_dir = tempfile::tempdir().unwrap();
        let input_file = temp_dir.path().join("input.txt");

        let test_db = TestDatabase {
            data: String::from(
                "trimmed_seq\toriginal_seq\treceptor_group\tepitope\torganism\tantigen\n\
             ASSLG\tCASSLGF\tgroup1\tepi1\torg1\tant1\n\
             ASSQDRDTQY\tASSQDRDTQY\tgroup2\tepi2\torg2\tant2\n",
            ),
            path: temp_dir.path().join("db.tsv"),
        };

        // A blank line, and "CF", which trimming reduces to an empty sequence.
        std::fs::write(&input_file, "CASSLGF\n\nCF\n").unwrap();

        let mut output = Vec::new();
        run_tcrmatch(
            &input_file, &test_db, 1, 0.0, false, false, false, None, Progress::NoProgress, &mut output,
        )
        .unwrap();

        let out_str = String::from_utf8(output).unwrap();
        let rows: Vec<&str> = out_str.lines().skip(1).collect();
        assert!(!rows.is_empty(), "expected matches for CASSLGF");
        assert!(
            rows.iter().all(|row| row.starts_with("ASSLG\t")),
            "an empty input sequence produced output rows:\n{out_str}"
        );
    }

    #[test]
    fn test_run_tcrmatch_trust4_rows_follow_input_lines() {
        let temp_dir = tempfile::tempdir().unwrap();
        let input_file = temp_dir.path().join("input.txt");

        let test_db = TestDatabase {
            data: String::from(
                "trimmed_seq\toriginal_seq\treceptor_group\tepitope\torganism\tantigen\n\
             ASSLG\tCASSLGF\tgroup1\tepi1\torg1\tant1\n\
             ASSQDRDTQY\tASSQDRDTQY\tgroup2\tepi2\torg2\tant2\n",
            ),
            path: temp_dir.path().join("db.tsv"),
        };

        // Lines 1 and 3 share a CDR3; every line matches both entries at threshold 0.
        let lines = [
            "10\t0.05\tAAAG\tCASSLGF\tTRBV2*01\tTRBD1*01\tTRBJ1*01\tTRBC1\tcid1\tfull",
            "7\t0.03\tAAAC\tCASSQDRDTQYF\tTRBV5*01\tTRBD1*01\tTRBJ2*01\tTRBC1\tcid2\tfull",
            "3\t0.01\tAAAT\tCASSLGF\tTRBV2*01\tTRBD2*01\tTRBJ1*01\tTRBC1\tcid3\tfull",
        ];
        std::fs::write(
            &input_file,
            format!(
                "#count\tfrequency\tCDR3nt\tCDR3aa\tV\tD\tJ\tC\tcid\tcid_full_length\n{}\n",
                lines.join("\n")
            ),
        )
        .unwrap();

        let mut output = Vec::new();
        run_tcrmatch(
            &input_file, &test_db, 1, 0.0, false, false, true, None, Progress::NoProgress, &mut output,
        )
        .unwrap();

        let out_str = String::from_utf8(output).unwrap();
        let rows: Vec<&str> = out_str.lines().skip(1).collect();
        assert_eq!(rows.len(), 6, "expected two rows per input line:\n{out_str}");
        for (i, line) in lines.iter().enumerate() {
            for row in &rows[2 * i..2 * i + 2] {
                assert!(row.starts_with(&format!("{line}\t")), "row out of order: {row}");
            }
        }
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
            Progress::NoProgress,
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
            Progress::NoProgress,
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
            Progress::NoProgress,
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
                Progress::NoProgress,
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
                Progress::NoProgress,
                &mut output,
            )
            .is_err()
        );
    }
}
