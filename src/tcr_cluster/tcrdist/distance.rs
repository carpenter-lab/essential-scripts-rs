// Exact match to old distance.rs

use super::lookup_cdr_seqs::{GAP, VGeneCdrs};
use crate::progress::{ProgressSettings, progress_bar_with_message};
use anyhow::Context;
use rayon::prelude::*;
use rust_igraph::Graph;
use std::collections::HashMap;
use std::fmt::Display;
use std::ops::{Add, AddAssign, Index, Mul, MulAssign};
use std::sync::atomic::{AtomicUsize, Ordering};

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
    /// The flattened lower triangle: row `i` (for `i >= 1`) holds the distances to records
    /// `0..i` and starts at offset `i * (i - 1) / 2`.
    pub(crate) fn lower_triangle(&self) -> &[u8] {
        &self.distances
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

/// CDR3 residues tcrdist3 leaves out of the comparison at the N- and C-terminus.
const CDR3_NTRIM: usize = 3;
const CDR3_CTRIM: usize = 2;

/// tcrdist3's CDR3 distance (pwseqdist `nb_tcrdist` with tcrdist3's CDR3 settings:
/// `ntrim=3, ctrim=2, fixed_gappos=False`), before the CDR3 weight is applied.
///
/// - Equal lengths: no gaps; the substitution distances at each position are summed,
///   skipping the first 3 and last 2 residues.
/// - Unequal lengths: one contiguous gap block of the length difference, costing
///   `gap_penalty` per missing residue, is placed in the longer sequence. Residues before
///   it are compared from the N-terminus and residues after it from the C-terminus, and
///   the block position giving the smallest distance is used. Positions searched are
///   `5..=short_len - 5`, shifted down symmetrically for sequences shorter than 10.
///
/// Saturates at 255.
#[inline]
pub fn tcrdist_cdr3(a: &[u8], b: &[u8], fast_matrix: &FastDistMatrixU8, gap_penalty: u8) -> u8 {
    let d = |x: u8, y: u8| u32::from(fast_matrix[x as usize][y as usize]);
    let (short, long) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    let (s_len, l_len) = (short.len(), long.len());
    let gap = (l_len - s_len) as u32 * u32::from(gap_penalty);
    if gap >= u32::from(u8::MAX) {
        return u8::MAX;
    }
    if s_len == l_len {
        let total: u32 = (CDR3_NTRIM..s_len.saturating_sub(CDR3_CTRIM))
            .map(|i| d(a[i], b[i]))
            .sum();
        return total.min(u32::from(u8::MAX)) as u8;
    }

    let (mut min_pos, mut max_pos) = (5isize, s_len as isize - 5);
    while min_pos > max_pos {
        min_pos -= 1;
        max_pos += 1;
    }
    let mut best = u32::MAX;
    for gap_pos in min_pos..=max_pos {
        let gap_pos = gap_pos as usize;
        let n_side: u32 = (CDR3_NTRIM..gap_pos).map(|i| d(short[i], long[i])).sum();
        let c_side: u32 = (CDR3_CTRIM..s_len - gap_pos)
            .map(|c| d(short[s_len - 1 - c], long[l_len - 1 - c]))
            .sum();
        best = best.min(n_side + c_side);
        if best == 0 {
            break;
        }
    }
    (best + gap).min(u32::from(u8::MAX)) as u8
}

fn compute_lower_triangle_matrix_u8(
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
    let total_pairs = n * (n - 1) / 2;
    let mut flat_distances = vec![u8::MAX; total_pairs];

    // Group sequence slices by byte length, so pairs whose length difference alone
    // exceeds `max_cap` are never visited.
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

    // Interleave small and large rows for thread load balancing
    let mut balanced_rows = Vec::with_capacity(n - 1);
    let mut iter = row_slices.into_iter();
    while let Some(front) = iter.next() {
        balanced_rows.push(front);
        if let Some(back) = iter.next_back() {
            balanced_rows.push(back);
        }
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

        for bucket in length_buckets.iter().take(max_bucket + 1).skip(min_bucket) {
            let cutoff = bucket.partition_point(|&(j, _)| j < i);
            for &(j, bytes_j) in &bucket[..cutoff] {
                row_slice[j] = tcrdist_cdr3(bytes_i, bytes_j, &fast_matrix, gap_penalty as u8);
            }
        }
        let prev = pending_progress.fetch_add(i, Ordering::Relaxed);
        if prev + i >= 100_000 {
            let accumulated = pending_progress.swap(0, Ordering::Relaxed);
            if accumulated > 0 {
                pb.inc(accumulated as u64);
            }
        }
    });

    pb.finish_and_clear();
    eprintln!(
        "Computed {} distances in {}",
        name,
        indicatif::HumanDuration(pb.elapsed())
    );

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
    compute_lower_triangle_matrix_u8(sequences, dist_matrix, gap_penalty, max_cap, name)
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
    pb.finish_and_clear();
    eprintln!(
        "Computed {} germline distances in {}",
        name,
        indicatif::HumanDuration(pb.elapsed())
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tcr_cluster::blosum::Chain;
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

    fn blosum_fast_matrix() -> Box<FastDistMatrixU8> {
        let m = crate::tcr_cluster::blosum::AminoAcidSubMatrix::Blosum62
            .build_mismatch_matrix(4, &Chain::Beta);
        build_fast_dist_matrix_u8(&m, 4)
    }

    /// Line-by-line port of pwseqdist's `nb_tcrdist` with tcrdist3's CDR3 settings
    /// (`ntrim=3, ctrim=2, fixed_gappos=False, dist_weight=1`), keeping Python's argument
    /// order and empty-range semantics.
    fn nb_tcrdist_reference(seq1: &[u8], seq2: &[u8], m: &FastDistMatrixU8, gap: u32) -> u32 {
        let d = |x: u8, y: u8| u32::from(m[x as usize][y as usize]);
        let (q_l, s_l) = (seq1.len() as isize, seq2.len() as isize);
        let (ntrim, ctrim) = (3isize, 2isize);
        if q_l == s_l {
            return (ntrim..q_l - ctrim)
                .map(|i| d(seq1[i as usize], seq2[i as usize]))
                .sum();
        }
        let short_len = q_l.min(s_l);
        let len_diff = (q_l - s_l).unsigned_abs() as u32;
        let (mut min_gappos, mut max_gappos) = (5isize, short_len - 1 - 4);
        while min_gappos > max_gappos {
            min_gappos -= 1;
            max_gappos += 1;
        }
        let mut min_dist: Option<u32> = None;
        for gappos in min_gappos..=max_gappos {
            let mut tmp = 0;
            let remainder = short_len - gappos;
            for n_i in ntrim..gappos {
                tmp += d(seq1[n_i as usize], seq2[n_i as usize]);
            }
            for c_i in ctrim..remainder {
                tmp += d(
                    seq1[(q_l - 1 - c_i) as usize],
                    seq2[(s_l - 1 - c_i) as usize],
                );
            }
            if min_dist.is_none_or(|m| tmp < m) {
                min_dist = Some(tmp);
            }
            if min_dist == Some(0) {
                break;
            }
        }
        min_dist.unwrap() + len_diff * gap
    }

    #[test]
    fn cdr3_distance_matches_hand_checked_tcrdist3_values() {
        let m = blosum_fast_matrix();
        let dist = |a: &str, b: &str| tcrdist_cdr3(a.as_bytes(), b.as_bytes(), &m, 4);
        // Differences only in the trimmed first 3 / last 2 residues are ignored.
        assert_eq!(dist("CASSLGQETQYF", "CSRSLGQETQWY"), 0);
        // One L->P substitution: min(4, 4 - (-3)) = 4.
        assert_eq!(dist("CASSLGQETQYF", "CASSPGQETQYF"), 4);
        // Inserted G inside the gap window: the gap alone.
        assert_eq!(dist("CASSLGQETQYF", "CASSLGGQETQYF"), 4);
        assert_eq!(dist("CASSLGGQETQYF", "CASSLGQETQYF"), 4);
        // Same length, shifted by one: no gaps allowed, so every shifted position counts
        // (a free alignment would pay one insertion + one deletion = 8).
        assert!(dist("CASSGLGQETQYF", "CASSLGQETQYFF") > 8);
        // An insertion next to the C-terminus can't put its gap there: the residues
        // between the gap window and the insertion are compared out of register.
        assert!(dist("CASSLGQETQYF", "CASSLGQETQAYF") > 4);
    }

    #[test]
    fn cdr3_distance_matches_pwseqdist_port() {
        let m = blosum_fast_matrix();
        let alphabet = b"ACDEFGHIKLMNPQRSTVWY";
        let mut state = 0x2545_F491u32;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            state
        };
        let mut random_seq = |len: usize| -> Vec<u8> {
            (0..len)
                .map(|_| alphabet[next() as usize % alphabet.len()])
                .collect()
        };
        for len_a in 0usize..20 {
            for len_b in len_a.saturating_sub(4)..=len_a + 4 {
                for _ in 0..20 {
                    let a = random_seq(len_a);
                    // Mostly-similar pairs so the gap search has close calls to get right.
                    let mut b = random_seq(len_b);
                    for (i, x) in b.iter_mut().enumerate().take(len_a.min(len_b)) {
                        if i % 3 != 0 {
                            *x = a[i];
                        }
                    }
                    let want = nb_tcrdist_reference(&a, &b, &m, 4).min(255);
                    assert_eq!(u32::from(tcrdist_cdr3(&a, &b, &m, 4)), want, "{a:?} {b:?}");
                    assert_eq!(u32::from(tcrdist_cdr3(&b, &a, &m, 4)), want, "symmetric");
                }
            }
        }
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
        let mismatch = crate::tcr_cluster::blosum::AminoAcidSubMatrix::Blosum62
            .build_mismatch_matrix(4, &Chain::Beta);
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
