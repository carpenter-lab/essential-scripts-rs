// Exact match to old distance.rs

use super::lookup_cdr_seqs::{GAP, VGeneCdrs};
use crate::progress::{ProgressSettings, progress_bar_with_message};
use anyhow::Context;
use rayon::prelude::*;
use rust_igraph::Graph;
use std::cell::RefCell;
use std::collections::HashMap;
use std::fmt::Display;
use std::ops::{Add, AddAssign, Index, Mul, MulAssign};
use std::sync::atomic::{AtomicUsize, Ordering};

use std::array;

/// Block width for the SIMD-friendly kernels below. Fixed-size blocks give LLVM a known
/// trip count with no bounds checks, so the inner loops auto-vectorize on stable Rust
/// (NEON on aarch64, SSE/AVX on x86_64 with `target-cpu=native`).
/// Should default to 16 if AVX2 is not available on x86_64
#[cfg(target_feature = "avx512bw")]
const LANES: usize = 64;
#[cfg(any(
    all(not(target_feature = "avx512bw"), target_feature = "avx2"),
    target_arch = "aarch64"
))]
const LANES: usize = 32;
#[cfg(not(any(
    target_feature = "avx512bw",
    target_feature = "avx2",
    target_arch = "aarch64"
)))]
const LANES: usize = 16;

/// Elements per rayon task for element-wise matrix operations.
const ELEMENTWISE_CHUNK: usize = 1 << 16;

pub struct TcrDistanceMatrix {
    distances: Vec<u8>,
    pub(crate) dim: usize,
    names: Vec<String>,
}

impl TcrDistanceMatrix {
    pub fn from_flat(distances: Vec<u8>, names: Vec<String>, n: usize) -> Self {
        TcrDistanceMatrix {
            distances,
            dim: n,
            names,
        }
    }

    pub fn new(rows: Vec<Vec<u8>>, names: Vec<String>, n: usize) -> Self {
        TcrDistanceMatrix {
            distances: rows.into_iter().flatten().collect(),
            dim: n,
            names,
        }
    }
    pub fn dist(&self) -> Vec<u8> {
        self.distances.clone()
    }
    pub fn names(&self) -> Vec<String> {
        self.names.clone()
    }

    pub fn to_igraph(&self, threshold: u8) -> anyhow::Result<Graph> {
        let n = u32::try_from(self.dim).context("Too many sequences for an igraph graph")?;

        // Undirected graph on all `dim` vertices, so sequences without any edge stay present.
        let mut graph = Graph::new(n, false)?;

        // Stream the lower triangle (i > j) straight into the graph: each pair is visited once,
        // so there are no duplicate edges, and no intermediate edge list is materialized.
        // This must stay a single `add_edges` call, since igraph re-sorts all edges on every call.
        // Edge ids follow this insertion order; an edge's weight can be recovered later as
        // `self[(to, from)]` instead of storing a weight vector.
        let edges = (1..self.dim).flat_map(move |i| {
            (0..i).filter_map(move |j| (self[(i, j)] <= threshold).then_some((i as u32, j as u32)))
        });
        graph.add_edges(edges)?;

        Ok(graph)
    }
    pub fn max(&self) -> u8 {
        self.distances.iter().max().copied().unwrap_or(0)
    }
}

/// `lhs[i] = lhs[i].saturating_add(rhs[i])`; compiles to `uqadd` (NEON) / `paddusb` (x86).
#[inline]
fn saturating_add_slice(lhs: &mut [u8], rhs: &[u8]) {
    debug_assert_eq!(lhs.len(), rhs.len());
    let (l_blocks, l_tail) = lhs.as_chunks_mut::<LANES>();
    let (r_blocks, r_tail) = rhs.as_chunks::<LANES>();
    for (a, b) in l_blocks.iter_mut().zip(r_blocks) {
        for i in 0..LANES {
            a[i] = a[i].saturating_add(b[i]);
        }
    }
    for (a, &b) in l_tail.iter_mut().zip(r_tail) {
        *a = a.saturating_add(b);
    }
}

/// `x = x.saturating_mul(rhs)` for every element. `u8::saturating_mul` has no SIMD
/// equivalent, so this widens to u16 (u8 * u8 <= 65_025 never overflows), clamps to 255
/// and narrows back, which vectorizes to widening multiply + min/saturating narrow.
#[inline]
fn saturating_scale_slice(data: &mut [u8], rhs: u8) {
    let r = u16::from(rhs);
    let scale = |x: u8| (u16::from(x) * r).min(255) as u8;
    let (blocks, tail) = data.as_chunks_mut::<LANES>();
    for block in blocks {
        for x in block.iter_mut() {
            *x = scale(*x);
        }
    }
    for x in tail {
        *x = scale(*x);
    }
}

impl AddAssign<&TcrDistanceMatrix> for TcrDistanceMatrix {
    /// In-place, parallel saturating sum; no new matrix is allocated.
    fn add_assign(&mut self, rhs: &TcrDistanceMatrix) {
        assert_eq!(
            self.dim, rhs.dim,
            "Cannot add matrices with differing dimensions ({} vs {})",
            self.dim, rhs.dim
        );
        self.distances
            .par_chunks_mut(ELEMENTWISE_CHUNK)
            .zip(rhs.distances.par_chunks(ELEMENTWISE_CHUNK))
            .for_each(|(lhs, rhs)| saturating_add_slice(lhs, rhs));
    }
}

impl AddAssign for TcrDistanceMatrix {
    fn add_assign(&mut self, rhs: Self) {
        *self += &rhs;
    }
}

impl<'b> Add<&'b TcrDistanceMatrix> for &TcrDistanceMatrix {
    type Output = TcrDistanceMatrix;

    fn add(self, rhs: &'b TcrDistanceMatrix) -> Self::Output {
        let mut out = TcrDistanceMatrix {
            distances: self.distances.clone(),
            dim: self.dim,
            names: self.names.clone(),
        };
        out += rhs;
        out
    }
}

impl Add for TcrDistanceMatrix {
    type Output = Self;

    /// Reuses `self`'s buffer.
    fn add(mut self, rhs: Self) -> Self::Output {
        self += &rhs;
        self
    }
}

impl MulAssign<u8> for TcrDistanceMatrix {
    /// In-place, parallel saturating scale.
    fn mul_assign(&mut self, rhs: u8) {
        match rhs {
            1 => {}
            0 => self.distances.fill(0),
            _ => self
                .distances
                .par_chunks_mut(ELEMENTWISE_CHUNK)
                .for_each(|chunk| saturating_scale_slice(chunk, rhs)),
        }
    }
}

impl Mul<u8> for &TcrDistanceMatrix {
    type Output = TcrDistanceMatrix;

    fn mul(self, rhs: u8) -> Self::Output {
        let mut out = TcrDistanceMatrix {
            distances: self.distances.clone(),
            dim: self.dim,
            names: self.names.clone(),
        };
        out *= rhs;
        out
    }
}

impl Mul<u8> for TcrDistanceMatrix {
    type Output = Self;

    /// Reuses `self`'s buffer.
    fn mul(mut self, rhs: u8) -> Self::Output {
        self *= rhs;
        self
    }
}

impl Display for TcrDistanceMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.dim {
            for j in 0..self.dim {
                let dist = self.index((i, j));
                write!(f, "{}\t", dist)?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

impl Index<(usize, usize)> for TcrDistanceMatrix {
    type Output = u8;

    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        if i == j {
            return &0;
        }

        // Ensure `i` is always the larger index to look up inside the lower triangle
        let (row, col) = if i > j { (i, j) } else { (j, i) };

        // Calculate the index in the flattened lower-triangular vector
        let index = row * (row - 1) / 2 + col;
        &self.distances[index]
    }
}

pub type FastDistMatrixU8 = [[u8; 256]; 256];

fn build_fast_dist_matrix_u8(
    dist_matrix: &HashMap<(u8, u8), i32>,
    default_penalty: i32,
) -> Box<FastDistMatrixU8> {
    let mut fast_matrix = Box::new([[default_penalty.min(255) as u8; 256]; 256]);
    for (&(a, b), &cost) in dist_matrix {
        fast_matrix[a as usize][b as usize] = cost.min(255) as u8;
    }
    fast_matrix
}

#[inline(always)]
pub fn needleman_wunsch_simd<const LANES: usize>(
    seq1: &[u8],
    seqs2: [&[u8]; LANES],
    fast_matrix: &FastDistMatrixU8,
    gap_penalty: u8,
    max_cap: u8,
    scratch: &mut Vec<[u8; LANES]>, // Pass as a mutable Vec so we can resize safely!
) -> [u8; LANES] {
    let m = seq1.len();
    let n = seqs2[0].len();
    let fail_score = max_cap.saturating_add(1);

    // O(1) Length Pruning across all lanes at once
    let len_diff = (m.abs_diff(n) as u8).saturating_mul(gap_penalty);
    if len_diff >= max_cap {
        return [fail_score; LANES];
    }

    let width = n + 1;

    // SAFE DYNAMIC RESIZING:
    // Ensure the scratch buffer is at least `width` elements long.
    if scratch.len() < width {
        scratch.resize(width, [0; LANES]);
    }
    let dp_slice = &mut scratch[..width];

    // Base cases: Row 0 initialization
    for (j, cell) in dp_slice.iter_mut().enumerate() {
        *cell = [(j as u8).saturating_mul(gap_penalty); LANES];
    }

    let gap_vec = [gap_penalty; LANES];

    // Main DP loop
    for i in 1..=m {
        let char1 = seq1[i - 1] as usize;
        let cost_row = &fast_matrix[char1];

        let mut prev_diag = dp_slice[0];
        let mut left = [(i as u8).saturating_mul(gap_penalty); LANES];
        dp_slice[0] = left;

        let mut row_min = left;

        for j in 1..=n {
            let char_dists: [u8; LANES] =
                array::from_fn(|lane| cost_row[seqs2[lane][j - 1] as usize]);

            let up = dp_slice[j];
            let mut best = [0u8; LANES];

            for lane in 0..LANES {
                let match_mismatch = prev_diag[lane].saturating_add(char_dists[lane]);
                let gap_seq2 = up[lane].saturating_add(gap_vec[lane]);
                let gap_seq1 = left[lane].saturating_add(gap_vec[lane]);

                let b = match_mismatch.min(gap_seq2).min(gap_seq1);
                best[lane] = b;
                row_min[lane] = row_min[lane].min(b);
            }

            prev_diag = up;
            left = best;
            dp_slice[j] = best;
        }

        // SIMD Early-Exit
        if row_min.iter().all(|&min_val| min_val >= max_cap) {
            return [fail_score; LANES];
        }
    }

    dp_slice[n]
}

#[inline(always)]
fn needleman_wunsch_u8(
    seq1: &[u8],
    seq2: &[u8],
    fast_matrix: &FastDistMatrixU8,
    gap_penalty: u8,
    max_cap: u8, // Typically 255 (or your highest planned threshold + 1)
    scratch: &mut Vec<u8>,
) -> u8 {
    let m = seq1.len();
    let n = seq2.len();

    // Length pruning
    let len_diff = (m.abs_diff(n) as u8).saturating_mul(gap_penalty);
    if len_diff >= max_cap {
        return max_cap;
    }

    let width = n + 1;
    if scratch.len() < width {
        scratch.resize(width, 0);
    }
    let dp = &mut scratch[..width];

    for (j, cell) in dp.iter_mut().enumerate() {
        *cell = (j as u8).saturating_mul(gap_penalty);
    }

    for i in 1..=m {
        let char1 = seq1[i - 1] as usize;
        let cost_row = &fast_matrix[char1];

        let mut prev_diag = dp[0];
        let mut left = (i as u8).saturating_mul(gap_penalty);
        dp[0] = left;
        let mut row_min = left;

        for j in 1..=n {
            let char2 = seq2[j - 1] as usize;
            let char_dist = cost_row[char2];

            let up = dp[j];
            let match_mismatch = prev_diag.saturating_add(char_dist);
            let gap_seq2 = up.saturating_add(gap_penalty);
            let gap_seq1 = left.saturating_add(gap_penalty);

            let best = match_mismatch.min(gap_seq2).min(gap_seq1);

            prev_diag = up;
            left = best;
            dp[j] = best;
            row_min = row_min.min(best);
        }

        // Early-exit if every possible path through this row exceeds max_cap
        if row_min >= max_cap {
            return max_cap;
        }
    }

    dp[n]
}

fn compute_lower_triangle_matrix_simd_u8<const LANES: usize>(
    sequences: &[String],
    dist_matrix: &HashMap<(u8, u8), i32>,
    gap_penalty: i32,
    max_cap: u8,
    name: &str,
) -> anyhow::Result<TcrDistanceMatrix> {
    let n = sequences.len();
    if n <= 1 {
        return Ok(TcrDistanceMatrix::from_flat(Vec::new(), Vec::new(), n));
    }

    let fast_matrix = build_fast_dist_matrix_u8(dist_matrix, gap_penalty);
    let fail_score = max_cap.saturating_add(1);
    let total_pairs = n * (n - 1) / 2;

    let mut flat_distances = vec![fail_score; total_pairs];

    // Group sequence slices by byte length
    let max_len = sequences.iter().map(|s| s.len()).max().unwrap_or(0);
    let mut length_buckets: Vec<Vec<(usize, &[u8])>> = vec![Vec::new(); max_len + 1];
    for (idx, seq) in sequences.iter().enumerate() {
        length_buckets[seq.len()].push((idx, seq.as_bytes()));
    }

    let max_len_diff = if gap_penalty > 0 {
        (max_cap as usize) / (gap_penalty as usize)
    } else {
        max_len
    };

    // Carve memory into disjoint row slices
    let mut row_slices: Vec<(usize, &mut [u8])> = Vec::with_capacity(n - 1);
    let mut remainder = flat_distances.as_mut_slice();
    for i in 1..n {
        let (row_slice, rest) = remainder.split_at_mut(i);
        row_slices.push((i, row_slice));
        remainder = rest;
    }

    // Interleave small and large rows for perfect thread load balancing
    let mut balanced_rows = Vec::with_capacity(n - 1);
    let mut iter = row_slices.into_iter();
    while let Some(front) = iter.next() {
        balanced_rows.push(front);
        if let Some(back) = iter.next_back() {
            balanced_rows.push(back);
        }
    }

    // Thread-local scratch buffers adapt dynamically to LANES
    thread_local! {
        static SIMD_SCRATCH: RefCell<Vec<[u8; 0]>> = const { RefCell::new(Vec::new()) }; // Stub fallback
        static SCALAR_SCRATCH: RefCell<Vec<u8>> = RefCell::new(Vec::with_capacity(64));
    }

    // Concrete scratch type allocation per thread
    thread_local! {
        // We use dynamic raw allocations inside RefCell to safely handle generic array sizes
        static GENERIC_SIMD_SCRATCH: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
    }

    let pb = progress_bar_with_message(
        total_pairs,
        &format!("Computing distances on {} sequences", name),
        ProgressSettings::Percent,
    )?;
    let pending_progress = AtomicUsize::new(0);

    balanced_rows.into_par_iter().for_each(|(i, row_slice)| {
        let bytes_i = sequences[i].as_bytes();
        let len_i = bytes_i.len();

        let min_bucket = len_i.saturating_sub(max_len_diff);
        let max_bucket = (len_i + max_len_diff).min(max_len);

        // 1. Allocate a fresh, 100% safe local DP row buffer for this thread.
        // Since this happens once per ROW (not per alignment!), there is virtually zero overhead.
        let mut simd_scratch = Vec::with_capacity(max_len + 1);

        SCALAR_SCRATCH.with(|scalar_cell| {
            let mut scalar_scratch = scalar_cell.borrow_mut();

            for bucket in length_buckets.iter().take(max_bucket + 1).skip(min_bucket) {
                let cutoff = bucket.partition_point(|&(j, _)| j < i);
                let valid_targets = &bucket[..cutoff];

                // A. Chunking by LANES
                let (chunks, remainder_targets) = valid_targets.as_chunks::<LANES>();

                for chunk in chunks {
                    let seqs2: [&[u8]; LANES] = array::from_fn(|lane| chunk[lane].1);

                    let scores = needleman_wunsch_simd::<LANES>(
                        bytes_i,
                        seqs2,
                        &fast_matrix,
                        gap_penalty as u8,
                        max_cap,
                        &mut simd_scratch, // Pure safe array slice!
                    );

                    for lane in 0..LANES {
                        row_slice[chunk[lane].0] = scores[lane];
                    }
                }

                // B. Remainder (0..LANES-1 targets)
                for &(j, bytes_j) in remainder_targets {
                    row_slice[j] = needleman_wunsch_u8(
                        bytes_i,
                        bytes_j,
                        &fast_matrix,
                        gap_penalty as u8,
                        max_cap,
                        &mut scalar_scratch,
                    );
                }
            }
        });
        let prev = pending_progress.fetch_add(i, Ordering::Relaxed);
        if prev + i >= 100_000 {
            let accumulated = pending_progress.swap(0, Ordering::Relaxed);
            if accumulated > 0 {
                pb.inc(accumulated as u64);
            }
        }
    });

    pb.finish();

    Ok(TcrDistanceMatrix::from_flat(
        flat_distances,
        sequences.to_vec(),
        n,
    ))
}

pub fn compute_matrix_auto_detect(
    sequences: &[String],
    dist_matrix: &HashMap<(u8, u8), i32>,
    gap_penalty: i32,
    max_cap: u8,
    name: &str,
) -> anyhow::Result<TcrDistanceMatrix> {
    compute_lower_triangle_matrix_simd_u8::<LANES>(
        sequences,
        dist_matrix,
        gap_penalty,
        max_cap,
        name,
    )
}

/// Position-wise distance between two IMGT-aligned germline CDRs (tcrdist's fixed gap
/// positions): gap vs gap costs 0, residue vs gap costs `gap_penalty`, residue vs residue
/// uses `fast_matrix`. Unequal lengths are padded with gaps.
fn aligned_cdr_distance(a: &str, b: &str, fast_matrix: &FastDistMatrixU8, gap_penalty: u8) -> u32 {
    let gap = GAP as u8;
    let (a, b) = (a.as_bytes(), b.as_bytes());
    (0..a.len().max(b.len()))
        .map(|i| {
            match (
                a.get(i).copied().unwrap_or(gap),
                b.get(i).copied().unwrap_or(gap),
            ) {
                (x, y) if x == gap && y == gap => 0,
                (x, y) if x == gap || y == gap => u32::from(gap_penalty),
                (x, y) => u32::from(fast_matrix[x as usize][y as usize]),
            }
        })
        .sum()
}

/// Germline distances between distinct V alleles: every record's index into the table,
/// the number of distinct alleles `u`, and the `u × u` weighted distance table.
fn germline_allele_table(
    cdrs: &[&VGeneCdrs],
    dist_matrix: &HashMap<(u8, u8), i32>,
    gap_penalty: u8,
    weight: u8,
) -> (Vec<usize>, usize, Vec<u8>) {
    let fast_matrix = build_fast_dist_matrix_u8(dist_matrix, i32::from(gap_penalty));
    let mut unique: Vec<&VGeneCdrs> = Vec::new();
    let mut index_of: HashMap<&VGeneCdrs, usize> = HashMap::new();
    let keys: Vec<usize> = cdrs
        .iter()
        .map(|&c| {
            *index_of.entry(c).or_insert_with(|| {
                unique.push(c);
                unique.len() - 1
            })
        })
        .collect();

    let u = unique.len();
    let mut table = vec![0u8; u * u];
    for i in 0..u {
        for j in 0..i {
            let (a, b) = (unique[i], unique[j]);
            let d = aligned_cdr_distance(&a.cdr1, &b.cdr1, &fast_matrix, gap_penalty)
                + aligned_cdr_distance(&a.cdr2, &b.cdr2, &fast_matrix, gap_penalty)
                + aligned_cdr_distance(&a.pmhc, &b.pmhc, &fast_matrix, gap_penalty);
            let d = (d * u32::from(weight)).min(u32::from(u8::MAX)) as u8;
            table[i * u + j] = d;
            table[j * u + i] = d;
        }
    }
    (keys, u, table)
}

/// Adds the germline (CDR1 + CDR2 + CDR2.5) distances of one chain, multiplied by
/// `weight`, into `total` in place.
///
/// The germline CDRs are fixed by the V allele, so distances are computed once per pair of
/// distinct alleles and then looked up for every pair of records. Writing straight into
/// `total` avoids allocating (and then summing) a second full matrix.
pub fn add_germline_distances(
    total: &mut TcrDistanceMatrix,
    cdrs: &[&VGeneCdrs],
    dist_matrix: &HashMap<(u8, u8), i32>,
    gap_penalty: u8,
    weight: u8,
    name: &str,
) -> anyhow::Result<()> {
    let n = cdrs.len();
    anyhow::ensure!(
        total.dim == n,
        "{name}: {n} records but the distance matrix has {}",
        total.dim
    );
    if n <= 1 {
        return Ok(());
    }
    let (keys, u, table) = germline_allele_table(cdrs, dist_matrix, gap_penalty, weight);

    let pb = progress_bar_with_message(
        n - 1,
        &format!("Computing distances on {name}"),
        ProgressSettings::Percent,
    )?;
    let mut rows: Vec<(usize, &mut [u8])> = Vec::with_capacity(n - 1);
    let mut remainder = total.distances.as_mut_slice();
    for i in 1..n {
        let (row, rest) = remainder.split_at_mut(i);
        rows.push((i, row));
        remainder = rest;
    }
    rows.into_par_iter().with_min_len(64).for_each(|(i, row)| {
        let lookup = &table[keys[i] * u..(keys[i] + 1) * u];
        for (cell, &key_j) in row.iter_mut().zip(&keys) {
            *cell = cell.saturating_add(lookup[key_j]);
        }
        pb.inc(1);
    });
    pb.finish();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    pub fn compute_germline_matrix(
        cdrs: &[&VGeneCdrs],
        dist_matrix: &HashMap<(u8, u8), i32>,
        gap_penalty: u8,
        weight: u8,
        name: &str,
    ) -> anyhow::Result<TcrDistanceMatrix> {
        let n = cdrs.len();
        let mut matrix =
            TcrDistanceMatrix::from_flat(vec![0u8; n * n.saturating_sub(1) / 2], Vec::new(), n);
        add_germline_distances(&mut matrix, cdrs, dist_matrix, gap_penalty, weight, name)?;
        Ok(matrix)
    }
    #[test]
    fn to_igraph_keeps_all_vertices_and_only_edges_within_threshold() {
        // Lower triangle order: (1,0) (2,0) (2,1) (3,0) (3,1) (3,2)
        let names = ["A", "B", "C", "D"].map(String::from).to_vec();
        let matrix = TcrDistanceMatrix::from_flat(vec![5, 50, 10, 200, 3, 100], names, 4);

        let graph = matrix.to_igraph(10).unwrap();

        assert_eq!(graph.vcount(), 4);
        let edges: Vec<(u32, u32)> = (0..graph.ecount() as u32)
            .map(|eid| graph.edge(eid).unwrap())
            .collect();
        // Edge ids follow insertion order; undirected edges are stored as (smaller, larger).
        assert_eq!(edges, vec![(0, 1), (1, 2), (1, 3)]);
        // Weights are recoverable from the matrix by edge endpoints.
        assert!(
            edges
                .iter()
                .all(|&(u, v)| matrix[(v as usize, u as usize)] <= 10)
        );
    }

    /// Deterministic pseudo-random bytes covering the full u8 range, incl. 0/1/255.
    fn test_bytes(len: usize, seed: u32) -> Vec<u8> {
        let mut state = seed.wrapping_mul(2_654_435_761).max(1);
        (0..len)
            .map(|i| match i % 17 {
                0 => 0,
                1 => 1,
                2 => 255,
                _ => {
                    state ^= state << 13;
                    state ^= state >> 17;
                    state ^= state << 5;
                    (state >> 24) as u8
                }
            })
            .collect()
    }

    #[test]
    fn simd_kernels_match_scalar_saturating_ops() {
        // Lengths below, at and around the block width, to exercise the tail handling.
        for len in [0, 1, 63, 64, 65, 127, 128, 1000, 4099] {
            let lhs = test_bytes(len, 1);
            let rhs = test_bytes(len, 2);

            let mut got = lhs.clone();
            saturating_add_slice(&mut got, &rhs);
            let want: Vec<u8> = lhs
                .iter()
                .zip(&rhs)
                .map(|(a, b)| a.saturating_add(*b))
                .collect();
            assert_eq!(got, want, "add, len {len}");

            for scale in [0, 1, 2, 3, 16, 127, 128, 255] {
                let mut got = lhs.clone();
                saturating_scale_slice(&mut got, scale);
                let want: Vec<u8> = lhs.iter().map(|a| a.saturating_mul(scale)).collect();
                assert_eq!(got, want, "mul by {scale}, len {len}");
            }
        }
    }

    #[test]
    fn matrix_ops_span_multiple_rayon_chunks() {
        let len = ELEMENTWISE_CHUNK * 2 + 37;
        let a = test_bytes(len, 3);
        let b = test_bytes(len, 4);
        let m = |d: Vec<u8>| TcrDistanceMatrix::from_flat(d, Vec::new(), 0);

        let sum = m(a.clone()) + m(b.clone());
        let want: Vec<u8> = a
            .iter()
            .zip(&b)
            .map(|(x, y)| x.saturating_add(*y))
            .collect();
        assert_eq!(sum.distances, want);

        for scale in [0, 1, 7] {
            let scaled = &m(a.clone()) * scale;
            let want: Vec<u8> = a.iter().map(|x| x.saturating_mul(scale)).collect();
            assert_eq!(scaled.distances, want, "scale {scale}");
        }
    }

    fn cdrs(cdr1: &str, cdr2: &str, pmhc: &str) -> VGeneCdrs {
        VGeneCdrs {
            cdr1: cdr1.to_owned(),
            cdr2: cdr2.to_owned(),
            pmhc: pmhc.to_owned(),
        }
    }

    #[test]
    fn germline_matrix_scores_gaps_and_mismatches_per_allele_pair() {
        let mismatch =
            crate::tcr_cluster::blosum::AminoAcidSubMatrix::Blosum62.build_mismatch_matrix(4);
        let a = cdrs("MN.HEY", "SVG", "TDQ");
        let b = cdrs("MNAHEY", "SVG", "TDQ"); // one residue vs gap
        let c = cdrs("MN.HEY", "SVG", "TDQ"); // identical to a
        let records = [&a, &b, &c, &b];

        let matrix = compute_germline_matrix(&records, &mismatch, 4, 2, "test").unwrap();

        assert_eq!(matrix[(0, 1)], 8, "gap penalty 4 x weight 2");
        assert_eq!(matrix[(0, 2)], 0);
        assert_eq!(matrix[(1, 3)], 0);
        assert_eq!(matrix[(2, 3)], 8);
        assert_eq!(matrix[(1, 2)], matrix[(2, 1)]);
    }
}
