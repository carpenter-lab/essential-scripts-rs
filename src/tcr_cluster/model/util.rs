pub struct CumulativeSampler {
    cumulative_weights: Vec<f64>,
    total: f64,
}

impl CumulativeSampler {
    pub(crate) fn new(weights: &[f64]) -> Self {
        let mut cumulative_weights = Vec::with_capacity(weights.len());
        let mut running_sum = 0.0;

        for &weight in weights {
            if weight > 0.0 && weight.is_finite() {
                running_sum += weight;
            }
            cumulative_weights.push(running_sum);
        }

        Self {
            cumulative_weights,
            total: running_sum,
        }
    }

    #[inline]
    pub(crate) fn sample(&self, state: &mut u64) -> usize {
        let len = self.cumulative_weights.len();
        if len == 0 {
            return 0;
        }

        *state = state.wrapping_mul(6364136223846793005).wrapping_add(1);

        if self.total <= 0.0 {
            return (*state as usize) % len;
        }

        let target = ((*state >> 11) as f64 / ((1_u64 << 53) as f64)) * self.total;

        // Binary Search O(log W) instead of linear scan O(W)
        match self.cumulative_weights.binary_search_by(|&val| {
            val.partial_cmp(&target)
                .unwrap_or(std::cmp::Ordering::Greater)
        }) {
            Ok(idx) => idx.min(len - 1),
            Err(idx) => idx.min(len - 1),
        }
    }
}

pub trait RoundDecimal {
    fn round_dp(self, places: u32) -> f64;
}

impl RoundDecimal for f64 {
    fn round_dp(self, places: u32) -> f64 {
        let scale = 10_f64.powi(places as i32);
        (self * scale).round() / scale
    }
}

#[derive(Clone, Debug)]
pub struct UnionFind {
    pub parent: Vec<usize>,
    pub rank: Vec<u8>,
}

impl UnionFind {
    pub fn new(size: usize) -> Self {
        Self {
            parent: (0..size).collect(),
            rank: vec![0; size],
        }
    }

    pub fn root(&mut self, node: usize) -> usize {
        if self.parent[node] != node {
            self.parent[node] = self.root(self.parent[node]);
        }
        self.parent[node]
    }

    pub fn join(&mut self, left: usize, right: usize) {
        let mut left = self.root(left);
        let mut right = self.root(right);
        if left == right {
            return;
        }
        if self.rank[left] < self.rank[right] {
            std::mem::swap(&mut left, &mut right);
        }
        self.parent[right] = left;
        if self.rank[left] == self.rank[right] {
            self.rank[left] += 1;
        }
    }
}

pub fn seeded_random(seed: &str) -> u64 {
    let mut value = 0xcbf29ce484222325_u64;
    for byte in seed.bytes() {
        value ^= u64::from(byte);
        value = value.wrapping_mul(0x100000001b3);
    }
    value.max(1)
}

pub fn random_uniform_index(state: &mut u64, length: usize) -> usize {
    *state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
    (*state as usize) % length.max(1)
}
