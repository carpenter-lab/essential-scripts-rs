//! Linkage straight from the distance matrix, for several cutoffs at once, without building
//! an igraph graph: single linkage (connected components of the "distance <= cutoff" graph)
//! and shared-neighbor linkage (the same, keeping only links whose neighborhoods overlap).

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

/// Every record's neighbors within `max_cutoff`, as `(index, distance)` sorted by index
/// (the record itself excluded). Built from one parallel pass over the matrix; memory is
/// one entry per neighbor, so `max_cutoff` must keep neighborhoods small.
pub fn neighbor_lists(matrix: &TcrDistanceMatrix, max_cutoff: u8) -> Vec<Vec<(u32, u8)>> {
    let n = matrix.dim;
    assert!(u32::try_from(n).is_ok(), "Too many records for linkage");
    let tri = matrix.lower_triangle();
    // Row `i` of the lower triangle holds the neighbors `j < i`, in increasing order.
    let lower: Vec<Vec<(u32, u8)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let start = i * i.saturating_sub(1) / 2;
            tri[start..start + i]
                .iter()
                .enumerate()
                .filter(|&(_, &d)| d <= max_cutoff)
                .map(|(j, &d)| (j as u32, d))
                .collect()
        })
        .collect();
    // Append each `i` to the lists of its lower neighbors. Rows are visited in increasing
    // `i`, so every list stays sorted: its own lower entries, then the larger indices.
    let mut lists: Vec<Vec<(u32, u8)>> = Vec::with_capacity(n);
    for (i, row) in lower.into_iter().enumerate() {
        for &(j, d) in &row {
            lists[j as usize].push((i as u32, d));
        }
        lists.push(row);
    }
    lists
}

/// Shared-neighbor linkage for each of `cutoffs` (in the order given; each must be <= the
/// cutoff `lists` were built with).
///
/// Records `i` and `j` within the cutoff are linked when the Jaccard index of their closed
/// neighborhoods (each record plus its neighbors within the cutoff) is at least
/// `min_overlap`; the result is the connected components of the kept links. Members of a
/// tight group share most neighbors and stay linked, as does an isolated pair (overlap 1),
/// while a single link between two dense groups shares few and is dropped. With
/// `min_overlap` 0 this is single linkage.
pub fn shared_neighbor_memberships(
    lists: &[Vec<(u32, u8)>],
    cutoffs: &[u8],
    min_overlap: f32,
) -> Vec<Vec<u32>> {
    let n = lists.len();
    cutoffs
        .iter()
        .map(|&cutoff| {
            let degree: Vec<u32> = lists
                .par_iter()
                .map(|list| list.iter().filter(|&&(_, d)| d <= cutoff).count() as u32)
                .collect();
            let links: Vec<(u32, u32)> = (0..n)
                .into_par_iter()
                .flat_map_iter(|i| {
                    let list_i = &lists[i];
                    let degree = &degree;
                    list_i
                        .iter()
                        .filter(move |&&(_, d)| d <= cutoff)
                        .take_while(move |&&(j, _)| (j as usize) < i)
                        .filter(move |&&(j, _)| {
                            let shared = common_within(list_i, &lists[j as usize], cutoff);
                            // Closed neighbourhoods: both also contain `i` and `j`.
                            let inter = shared + 2;
                            let union = degree[i] + degree[j as usize] + 2 - inter;
                            inter as f32 >= min_overlap * union as f32
                        })
                        .map(move |&(j, _)| (i as u32, j))
                })
                .collect();
            let mut forest = UnionFind::new(n);
            for (i, j) in links {
                forest.union(i, j);
            }
            forest.memberships()
        })
        .collect()
}

/// Number of indices present in both sorted lists with a distance <= `cutoff` in each.
fn common_within(a: &[(u32, u8)], b: &[(u32, u8)], cutoff: u8) -> u32 {
    let (mut x, mut y, mut count) = (0, 0, 0);
    while x < a.len() && y < b.len() {
        let ((ia, da), (ib, db)) = (a[x], b[y]);
        match ia.cmp(&ib) {
            std::cmp::Ordering::Less => x += 1,
            std::cmp::Ordering::Greater => y += 1,
            std::cmp::Ordering::Equal => {
                if da <= cutoff && db <= cutoff {
                    count += 1;
                }
                x += 1;
                y += 1;
            }
        }
    }
    count
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tcr_cluster::tcrdist::clustering::graph_components;
    use std::collections::HashSet;

    fn matrix(n: usize, flat: Vec<u8>) -> TcrDistanceMatrix {
        let names = (0..n).map(|i| i.to_string()).collect();
        TcrDistanceMatrix::from_flat(flat, names, n)
    }

    #[test]
    fn shared_neighbors_cut_bridges_but_keep_pairs() {
        // Two triangles {0,1,2} and {3,4,5} joined by the single link 2-3, plus an isolated
        // pair {6,7}. Every listed link is at distance 1.
        let n = 8;
        let linked = [
            (1, 0),
            (2, 0),
            (2, 1),
            (3, 2),
            (4, 3),
            (5, 3),
            (5, 4),
            (7, 6),
        ];
        let mut flat = Vec::new();
        for i in 1..n {
            for j in 0..i {
                flat.push(if linked.contains(&(i, j)) { 1 } else { 99 });
            }
        }
        let m = matrix(n, flat);
        let lists = neighbor_lists(&m, 1);
        let single = &single_linkage_memberships(&m, &[1])[0];
        assert_eq!(single[0], single[5], "single linkage chains through 2-3");

        // Bridge 2-3: closed neighborhoods {0,1,2,3} and {2,3,4,5} share 2 of 6 = 0.33.
        // Triangle links share 3 of 4 or more; the isolated pair shares 2 of 2.
        let got = &shared_neighbor_memberships(&lists, &[1], 0.5)[0];
        assert_eq!(got[0], got[1]);
        assert_eq!(got[0], got[2]);
        assert_eq!(got[3], got[4]);
        assert_eq!(got[3], got[5]);
        assert_ne!(got[2], got[3], "the bridge is cut");
        assert_eq!(got[6], got[7], "an isolated pair stays linked");
        assert_eq!(
            shared_neighbor_memberships(&lists, &[1], 0.3)[0][2],
            shared_neighbor_memberships(&lists, &[1], 0.3)[0][3],
            "a lower overlap keeps the bridge"
        );
        // Cutoff 0: nothing is linked.
        let none = &shared_neighbor_memberships(&lists, &[0], 0.5)[0];
        assert_eq!(none.iter().collect::<HashSet<_>>().len(), n);
    }

    #[test]
    fn zero_overlap_is_single_linkage() {
        let n = 250;
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
        let cutoffs = [0, 1, 3, 6];
        let lists = neighbor_lists(&m, 6);
        assert_eq!(
            shared_neighbor_memberships(&lists, &cutoffs, 0.0),
            single_linkage_memberships(&m, &cutoffs)
        );
        // Lists hold exactly the neighbors within the cutoff, sorted.
        for (i, list) in lists.iter().enumerate() {
            let want: Vec<(u32, u8)> = (0..n)
                .filter(|&j| j != i && m[(i, j)] <= 6)
                .map(|j| (j as u32, m[(i, j)]))
                .collect();
            assert_eq!(list, &want);
        }
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
