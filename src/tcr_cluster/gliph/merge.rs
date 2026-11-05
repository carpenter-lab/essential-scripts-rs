use super::*;

pub(super) fn merged_gliph2_pattern(
    local_patterns: &BTreeSet<String>,
    global_patterns: &BTreeSet<String>,
    fallback_members: &BTreeSet<String>,
) -> String {
    if let Some(pattern) = local_patterns
        .iter()
        .max_by(|left, right| left.len().cmp(&right.len()).then_with(|| right.cmp(left)))
    {
        return pattern.into();
    }

    if let Some(pattern) = global_patterns
        .iter()
        .max_by(|left, right| left.len().cmp(&right.len()).then_with(|| right.cmp(left)))
    {
        return pattern.into();
    }

    fallback_members
        .iter()
        .next()
        .cloned()
        .unwrap_or_else(|| "unknown".to_owned())
}

pub(super) fn merge_gliph2_candidates(
    candidates: Vec<Cluster>,
    _sample_sequences: &[String],
) -> Vec<Cluster> {
    if candidates.is_empty() {
        return Vec::new();
    }

    let member_sets: Vec<BTreeSet<String>> = candidates
        .iter()
        .map(|cluster| cluster.member_sequences.iter().cloned().collect())
        .collect();

    let mut union = UnionFind::new(candidates.len());

    for left in 0..candidates.len() {
        for right in left + 1..candidates.len() {
            let left_kind = candidates[left].kind;
            let right_kind = candidates[right].kind;

            // Preserve local-local and global-global candidates as separate groups.
            // Only merge candidate groups when a local motif group and a global
            // similarity group substantially describe the same member set.
            if left_kind == right_kind {
                continue;
            }

            let shared = member_sets[left].intersection(&member_sets[right]).count();
            if shared < MIN_SHARED_MEMBERS
                || (shared == member_sets[left].len() && shared == member_sets[right].len())
            {
                continue;
            }

            let smaller = member_sets[left].len().min(member_sets[right].len());
            let overlap_fraction = shared as f64 / smaller.max(1) as f64;
            if overlap_fraction >= MIN_OVERLAP_FRACTION {
                union.join(left, right);
            }
        }
    }

    let mut components: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for index in 0..candidates.len() {
        let root = union.root(index);
        components.entry(root).or_default().push(index);
    }

    let mut merged = Vec::new();
    for indexes in components.into_values() {
        let mut members = BTreeSet::new();
        let mut kinds = BTreeSet::new();
        let mut local_patterns = BTreeSet::new();
        let mut global_patterns = BTreeSet::new();
        let mut fisher_score = 1.0_f64;

        for index in indexes {
            let cluster = &candidates[index];
            members.extend(cluster.member_sequences.iter().cloned());
            kinds.insert(cluster.kind);
            match cluster.kind {
                ClusterKind::Local => local_patterns.insert(cluster.pattern.clone()),
                ClusterKind::Global => global_patterns.insert(cluster.pattern.clone()),
                ClusterKind::Single => false,
            };

            fisher_score = fisher_score.min(cluster.fisher_score);
        }

        if members.len() < 2 {
            continue;
        }

        let kind = if kinds.contains(&ClusterKind::Local) && kinds.contains(&ClusterKind::Global) {
            ClusterKind::Global
        } else {
            kinds.iter().next().cloned().unwrap_or(ClusterKind::Global)
        };

        merged.push(Cluster::new(
            merged_gliph2_pattern(&local_patterns, &global_patterns, &members),
            kind,
            fisher_score,
            members.into_iter().collect(),
        ))
    }

    merged
}
