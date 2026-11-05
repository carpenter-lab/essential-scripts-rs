use super::*;

pub(super) fn hamming(left: &str, right: &str) -> Option<usize> {
    if left.len() != right.len() {
        return None;
    }
    Some(
        left.bytes()
            .zip(right.bytes())
            .filter(|(left, right)| left != right)
            .count(),
    )
}

pub(super) fn gliph1_candidates(
    mut candidates: Vec<Cluster>,
    sample_sequences: &[String],
    gccutoff: usize,
) -> Vec<Cluster> {
    let sequence_index = sample_sequences
        .iter()
        .enumerate()
        .map(|(index, sequence)| (sequence.as_str(), index))
        .collect::<HashMap<_, _>>();
    let mut union = UnionFind::new(sample_sequences.len());
    let mut edge_types: HashMap<(usize, usize), BTreeSet<ClusterKind>> = HashMap::new();
    let mut edge_fisher: HashMap<(usize, usize), f64> = HashMap::new();
    let progress = progress_bar_with_message(
        candidates.len() + sample_sequences.len(),
        "Clustering GLIPH1 network",
        true,
        true,
    )
    .expect("Failed to create progress bar for GLIPH1 network");
    for cluster in candidates.drain(..) {
        let indexes = cluster
            .member_sequences
            .iter()
            .filter_map(|sequence| sequence_index.get(sequence.as_str()).copied())
            .collect::<Vec<_>>();
        for i in 0..indexes.len() {
            for j in i + 1..indexes.len() {
                let (left, right) = (indexes[i].min(indexes[j]), indexes[i].max(indexes[j]));
                union.join(left, right);
                edge_types
                    .entry((left, right))
                    .or_default()
                    .insert(cluster.kind);
                edge_fisher
                    .entry((left, right))
                    .and_modify(|score| *score = score.min(cluster.fisher_score))
                    .or_insert(cluster.fisher_score);
            }
        }
        progress.inc(1);
    }
    for left in 0..sample_sequences.len() {
        for right in left + 1..sample_sequences.len() {
            if hamming(
                stats::motif_region(&sample_sequences[left]),
                stats::motif_region(&sample_sequences[right]),
            )
            .is_some_and(|distance| distance <= gccutoff)
            {
                union.join(left, right);
                edge_types
                    .entry((left, right))
                    .or_default()
                    .insert(ClusterKind::Global);
                edge_fisher.entry((left, right)).or_insert(1.0);
            }
            progress.inc(1);
        }
    }
    progress.finish_with_message("GLIPH1 clustering complete");
    let mut components: BTreeMap<usize, BTreeSet<String>> = BTreeMap::new();
    for (index, sequence) in sample_sequences.iter().enumerate() {
        let root = union.root(index);
        components.entry(root).or_default().insert(sequence.clone());
    }
    let mut out = Vec::new();
    for (root, members) in components
        .into_iter()
        .filter(|(_, members)| members.len() >= 2)
    {
        let kinds = edge_types
            .iter()
            .filter(|((left, right), _)| {
                members.contains(&sample_sequences[*left])
                    && members.contains(&sample_sequences[*right])
            })
            .flat_map(|(_, kinds)| kinds.iter().cloned())
            .collect::<BTreeSet<_>>();
        let kind = if kinds.len() > 1 {
            ClusterKind::Local
        } else {
            kinds.iter().next().cloned().unwrap_or(ClusterKind::Global)
        };
        out.push(Cluster::new(
            "not_implemented".into(),
            kind,
            edge_fisher
                .iter()
                .filter(|((left, right), _)| {
                    union.root(*left) == root && union.root(*right) == root
                })
                .map(|(_, score)| *score)
                .fold(1.0, f64::min),
            members.into_iter().collect(),
        ));
    }
    out
}
