use super::*;

/// Adds a [`ClusterKind::Single`] cluster for every CDR3β carried by at least two records
/// (even when it is also in a motif cluster) and for every other CDR3β that no cluster
/// contains, so every record appears in the output at least once.
pub(super) fn add_singleton_clusters(
    mut clusters: Vec<Cluster>,
    records: &[ChainRecord],
    reference_sequences: &[String],
) -> Vec<Cluster> {
    let clustered: HashSet<&str> = clusters
        .iter()
        .flat_map(|cluster| cluster.member_sequences.iter().map(String::as_str))
        .collect();
    let mut beta_counts: BTreeMap<&str, usize> = BTreeMap::new();
    for record in records {
        *beta_counts.entry(record.beta.as_str()).or_default() += 1;
    }
    let mut reference_counts: HashMap<&str, usize> = HashMap::new();
    for sequence in reference_sequences {
        if beta_counts.contains_key(sequence.as_str()) {
            *reference_counts.entry(sequence.as_str()).or_default() += 1;
        }
    }
    let log_fact_table = LogFactorialTable::new(reference_sequences.len() + records.len());

    let singles: Vec<Cluster> = beta_counts
        .into_iter()
        .filter(|&(beta, observed)| observed >= 2 || !clustered.contains(beta))
        .map(|(beta, observed)| {
            let fisher = ContingencyTable2x2::from_reported_counts(
                observed,
                records.len(),
                reference_counts.get(beta).copied().unwrap_or(0),
                reference_sequences.len(),
            )
            .map_or(1.0, |table| log_fact_table.fishers_exact(&table).0);
            Cluster::new(
                beta.to_owned(),
                ClusterKind::Single,
                fisher,
                vec![beta.to_owned()],
            )
        })
        .collect();
    eprintln!(
        "Added {} single clusters of identical CDR3β sequences",
        singles.len()
    );

    clusters.extend(singles);
    clusters
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(beta: &str) -> ChainRecord {
        ChainRecord {
            beta: beta.to_owned(),
            ..ChainRecord::default()
        }
    }

    #[test]
    fn groups_identical_cdr3b_seen_at_least_twice() {
        let records = [
            record("CASSIRSSYEQYF"),
            record("CASSLGGYGYTF"),
            record("CASSIRSSYEQYF"),
            record("CASSIRSSYEQYF"),
            record("CASSPLAGGYEQYF"),
        ];
        let reference = vec!["CASSPRSSYEQYF".to_owned(); 10];
        let existing = vec![Cluster::new(
            "SSIR".to_owned(),
            ClusterKind::Local,
            0.01,
            vec!["CASSIRSSYEQYF".to_owned(), "CASSPLAGGYEQYF".to_owned()],
        )];

        let clusters = add_singleton_clusters(existing, &records, &reference);
        let singles: Vec<&Cluster> = clusters
            .iter()
            .filter(|c| c.kind == ClusterKind::Single)
            .collect();
        let patterns: Vec<&str> = singles.iter().map(|c| c.pattern.as_str()).collect();
        assert_eq!(
            patterns,
            ["CASSIRSSYEQYF", "CASSLGGYGYTF"],
            "shared CDR3β is a single even when clustered; unclustered seen-once CDR3β too"
        );
        assert_eq!(
            singles[0].member_sequences,
            vec!["CASSIRSSYEQYF".to_owned()]
        );
        assert!(singles[0].fisher_score < 1.0);
    }
}
