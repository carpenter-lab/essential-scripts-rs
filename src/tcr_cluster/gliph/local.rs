use super::*;

pub(super) fn local_candidates(
    method: Algorithm,
    sample_sequences: &[String],
    reference_sequences: &[String],
    sim_depth: usize,
    config: &GlobalFilterConfig,
) -> Vec<Cluster> {
    let stats = stats::motif_statistics_iter(sample_sequences, reference_sequences);

    let pb = progress_bar_with_message(
        sample_sequences.len(),
        "Testing motif enrichment",
        true,
        true,
    )
    .expect("Failed to create progress bar for motif enrichment testing");
    let tic = std::time::Instant::now();
    let clusters: Vec<Cluster> = stats
        .filter_map(|(motif, stat)| {
            let (p_value, fold_change) = stats::get_p_value_and_fc(
                method,
                &motif,
                &stat,
                sim_depth,
                sample_sequences,
                reference_sequences,
            );

            let threshold = MOTIF_MINIMUM_FOLD_CHANGE[motif.len().saturating_sub(2).min(2)];
            let insignificant = if method == Algorithm::Gliph1 {
                p_value >= 0.05
            } else {
                p_value > 0.05
            };
            if stat.members.len() < config.min_members || insignificant || fold_change < threshold {
                return None;
            }
            Some(Cluster::new(
                motif,
                ClusterKind::Local,
                p_value,
                stat.members.into_iter().collect(),
            ))
        })
        .progress_with(pb)
        .collect();
    eprintln!(
        "Tested {} clusters in {}",
        clusters.len(),
        indicatif::HumanDuration(std::time::Instant::now() - tic)
    );
    clusters
}
