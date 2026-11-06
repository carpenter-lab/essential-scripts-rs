//! Single-linkage clustering straight from the distance matrix: the connected components of
//! the "distance <= cutoff" graph, for several cutoffs at once, without building any graph.

use super::distance::TcrDistanceMatrix;
use rayon::prelude::*;

/// Disjoint-set forest with union by size and path halving.
#[derive(Clone)]
struct UnionFind {
    parent: Vec<u32>,
    size: Vec<u32>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n as u32).collect(),
            size: vec![1; n],
        }
    }

    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            let grandparent = self.parent[self.parent[x as usize] as usize];
            self.parent[x as usize] = grandparent;
            x = grandparent;
        }
        x
    }

    fn union(&mut self, a: u32, b: u32) {
        let (mut a, mut b) = (self.find(a), self.find(b));
        if a == b {
            return;
        }
        if self.size[a as usize] < self.size[b as usize] {
            std::mem::swap(&mut a, &mut b);
        }
        self.parent[b as usize] = a;
        self.size[a as usize] += self.size[b as usize];
    }

    /// Adds every union recorded in `other` (same vertex set). Each parent link joins two
    /// vertices of one set, so replaying the links rebuilds `other`'s partition in `self`.
    fn absorb(&mut self, other: &UnionFind) {
        for (v, &p) in other.parent.iter().enumerate() {
            if v as u32 != p {
                self.union(v as u32, p);
            }
        }
    }

    /// Component label of every vertex, numbered in order of each component's lowest
    /// vertex (the same labels igraph's `connected_components` gives).
    fn memberships(&mut self) -> Vec<u32> {
        let mut label = vec![u32::MAX; self.parent.len()];
        let mut next = 0;
        (0..self.parent.len() as u32)
            .map(|v| {
                let root = self.find(v) as usize;
                if label[root] == u32::MAX {
                    label[root] = next;
                    next += 1;
                }
                label[root]
            })
            .collect()
    }
}

/// Single-linkage component memberships for each of `cutoffs` (in the order given), where
/// records at distance <= the cutoff are linked.
///
/// One parallel pass over the matrix serves every cutoff: each pair is unioned only at the
/// smallest cutoff that admits it, and each larger cutoff then absorbs the one below it.
/// Memory is O(cutoffs × records) per worker; no edge list or graph is built.
pub fn single_linkage_memberships(matrix: &TcrDistanceMatrix, cutoffs: &[u8]) -> Vec<Vec<u32>> {
    let n = matrix.dim;
    assert!(
        u32::try_from(n).is_ok(),
        "Too many records for single linkage"
    );
    let mut levels = cutoffs.to_vec();
    levels.sort_unstable();
    levels.dedup();

    // Smallest level whose cutoff admits distance `d`; `NONE` when `d` exceeds every cutoff.
    const NONE: u8 = u8::MAX;
    assert!(levels.len() < usize::from(NONE), "Too many cutoffs");
    let mut level_of = [NONE; 256];
    for (d, level) in level_of.iter_mut().enumerate() {
        if let Some(k) = levels.iter().position(|&c| d <= usize::from(c)) {
            *level = k as u8;
        }
    }

    let tri = matrix.lower_triangle();
    let new_levels = || vec![UnionFind::new(n); levels.len()];
    // Few, large row blocks: every block gets its own forests, merged pairwise at the end.
    let min_rows = (n / (rayon::current_num_threads() * 4)).max(1);
    let mut forests = (1..n)
        .into_par_iter()
        .with_min_len(min_rows)
        .fold(new_levels, |mut forests, i| {
            let start = i * (i - 1) / 2;
            for (j, &d) in tri[start..start + i].iter().enumerate() {
                let k = level_of[usize::from(d)];
                if k != NONE {
                    forests[usize::from(k)].union(i as u32, j as u32);
                }
            }
            forests
        })
        .reduce(new_levels, |mut a, b| {
            for (x, y) in a.iter_mut().zip(&b) {
                x.absorb(y);
            }
            a
        });

    // Every pair within a smaller cutoff is also within the larger ones.
    for k in 1..forests.len() {
        let (lower, upper) = forests.split_at_mut(k);
        upper[0].absorb(&lower[k - 1]);
    }
    let per_level: Vec<Vec<u32>> = forests.iter_mut().map(UnionFind::memberships).collect();
    cutoffs
        .iter()
        .map(|c| per_level[levels.binary_search(c).expect("cutoff is a level")].clone())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tcr_cluster::tcrdist::clustering::graph_components;

    fn matrix(n: usize, flat: Vec<u8>) -> TcrDistanceMatrix {
        let names = (0..n).map(|i| i.to_string()).collect();
        TcrDistanceMatrix::from_flat(flat, names, n)
    }

    #[test]
    fn links_chains_and_respects_each_cutoff() {
        // Lower triangle order: (1,0) (2,0) (2,1) (3,0) (3,1) (3,2) (4,0) (4,1) (4,2) (4,3)
        // 0-1 at 5, 1-2 at 10 (0 and 2 are far apart, joined only by chaining), 3-4 at 20.
        let m = matrix(5, vec![5, 200, 10, 255, 255, 255, 255, 255, 255, 20]);
        let got = single_linkage_memberships(&m, &[20, 4, 5, 10]);
        assert_eq!(got[0], vec![0, 0, 0, 1, 1]);
        assert_eq!(got[1], vec![0, 1, 2, 3, 4]);
        assert_eq!(got[2], vec![0, 0, 1, 2, 3]);
        assert_eq!(got[3], vec![0, 0, 0, 1, 2]);
    }

    #[test]
    fn matches_igraph_connected_components() {
        let n = 300;
        let mut state = 0x9E37_79B9u32;
        let flat: Vec<u8> = (0..n * (n - 1) / 2)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 17;
                state ^= state << 5;
                (state % 256) as u8
            })
            .collect();
        let m = matrix(n, flat);
        let cutoffs = [0, 1, 2, 3, 5, 254];
        for (cutoff, got) in cutoffs.iter().zip(single_linkage_memberships(&m, &cutoffs)) {
            let graph = m.to_igraph(*cutoff).unwrap();
            assert_eq!(got, graph_components(&graph).unwrap(), "cutoff {cutoff}");
        }
    }

    #[test]
    fn handles_empty_and_single_record() {
        assert!(single_linkage_memberships(&matrix(0, vec![]), &[5])[0].is_empty());
        assert_eq!(
            single_linkage_memberships(&matrix(1, vec![]), &[5]),
            vec![vec![0]]
        );
    }
}
