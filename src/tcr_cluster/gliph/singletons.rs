use super::*;

pub(super) fn add_singleton_clusters(
    mut clusters: Vec<Cluster>,
    sample_sequences: &[String],
    reference_sequences: &[String],
) -> Vec<Cluster> {
    let clustered_sequences = clusters
        .iter()
        .flat_map(|cluster| cluster.member_sequences.iter().cloned())
        .collect::<HashSet<_>>();
    let log_fact_table = LogFactorialTable::new(reference_sequences.len() + sample_sequences.len());
    let pb = progress_bar_with_message(
        sample_sequences.len(),
        "Adding singleton clusters",
        true,
        true,
    )
    .expect("Failed to create progress bar for adding singleton clusters");

    let mut singletons = sample_sequences
        .into_par_iter()
        .filter(|s| !clustered_sequences.contains(*s))
        .map(|sequence| {
            let observed = sample_sequences.iter().filter(|s| s.eq(&sequence)).count();
            let reference_count = reference_sequences
                .iter()
                .filter(|s| s.eq(&sequence))
                .count()
                + observed;
            let cont_table = ContingencyTable2x2::from_reported_counts(
                observed,
                sample_sequences.len(),
                reference_count,
                reference_sequences.len(),
            )
            .expect("Failed to create contingency table for global cluster");
            let (fisher, _) = log_fact_table.fishers_exact(&cont_table);
            pb.inc(1);
            Cluster::new(
                sequence.clone(),
                ClusterKind::Single,
                fisher,
                vec![sequence.clone()],
            )
        })
        .collect();

    clusters.append(&mut singletons);
    clusters
}
