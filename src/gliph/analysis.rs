use crate::gliph::blosum::{AminoAcidSubMatrix, BLOSUM_ORDER};
use crate::gliph::stats::{ContingencyTable2x2, LogFactorialTable};
use crate::gliph::types::*;
use crate::gliph::{blosum, score, stats};
use crate::progress::progress_bar_with_message;
use anyhow::{Result, bail};
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

const MOTIF_MINIMUM_FOLD_CHANGE: [f64; 3] = [1000.0, 100.0, 10.0];
const MIN_SHARED_MEMBERS: usize = 4;
const MIN_OVERLAP_FRACTION: f64 = 0.80;

fn is_amino_acid_sequence(sequence: &str) -> bool {
    sequence.bytes().all(|aa| BLOSUM_ORDER.contains(&aa))
}

fn filtered_sequences(records: &[ChainRecord]) -> Result<Vec<ChainRecord>> {
    let valid: Vec<ChainRecord> = records
        .iter()
        .filter(|record| {
            let sequence = record.beta.as_str();
            sequence.len() >= 8
                && sequence.starts_with('C')
                && sequence.ends_with('F')
                && is_amino_acid_sequence(sequence)
        })
        .cloned()
        .collect();
    if valid.is_empty() {
        bail!("No valid CDR3b sequences remain after GLIPH sequence filtering");
    }
    Ok(valid)
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

fn local_candidates(
    method: Algorithm,
    sample_sequences: &[String],
    reference_sequences: &[String],
    sim_depth: usize,
    config: &GlobalFilterConfig,
) -> Vec<Cluster> {
    let stats = stats::motif_statistics_par_iter(sample_sequences, reference_sequences);

    let pb = progress_bar_with_message(sample_sequences.len(), "Testing motif enrichment")
        .expect("Failed to create progress bar for motif enrichment testing");
    let tic = std::time::Instant::now();
    let clusters: Vec<Cluster> = stats
        .progress_with(pb)
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
        .collect();
    eprintln!(
        "Tested {} clusters in {}",
        clusters.len(),
        indicatif::HumanDuration(std::time::Instant::now() - tic)
    );
    clusters
}

fn global_tags(sequences: &[String]) -> BTreeMap<String, Vec<TaggedMember>> {
    let mut result: BTreeMap<String, Vec<TaggedMember>> = BTreeMap::new();
    for sequence in sequences {
        let region = stats::motif_region(sequence).as_bytes();
        for position in 0..region.len() {
            let mut tag = region.to_vec();
            let residue = tag[position];
            tag[position] = b'%';
            result
                .entry(String::from_utf8_lossy(&tag).to_string())
                .or_default()
                .push(TaggedMember {
                    sequence: sequence.clone(),
                    residue,
                    position,
                });
        }
    }
    eprintln!(
        "Indexed {} global tags from {} sample sequences",
        result.len(),
        sequences.len()
    );
    for members in result.values_mut() {
        let mut seen = HashSet::new();
        members.retain(|member| seen.insert((member.sequence.clone(), member.position)));
    }
    result
}

fn reference_tag_counts(
    sequences: &[String],
    relevant_tags: &HashSet<String>,
) -> HashMap<String, HashMap<String, u8>> {
    let mut counts: HashMap<String, HashMap<String, u8>> = HashMap::new();
    for sequence in sequences {
        let region = stats::motif_region(sequence).as_bytes();
        for position in 0..region.len() {
            let mut tag = region.to_vec();
            let residue = tag[position];
            tag[position] = b'%';
            let tag = String::from_utf8_lossy(&tag).to_string();
            if relevant_tags.contains(&tag) {
                counts
                    .entry(tag)
                    .or_default()
                    .insert(sequence.clone(), residue);
            }
        }
    }

    eprintln!(
        "Indexed {} reference global tags from {} sample sequences",
        counts.len(),
        sequences.len()
    );
    counts
}

fn global_fisher_candidates(
    sample_sequences: &[String],
    reference_sequences: &[String],
    config: &GlobalFilterConfig,
    blosum: &Option<AminoAcidSubMatrix>,
) -> Vec<Cluster> {
    let sample_tags = global_tags(sample_sequences)
        .into_iter()
        .collect::<Vec<_>>();
    let log_fact_table = LogFactorialTable::new(reference_sequences.len() + sample_sequences.len());
    let relevant_tags: HashSet<String> = sample_tags.iter().map(|(tag, _)| tag.clone()).collect();
    let reference_tags = reference_tag_counts(reference_sequences, &relevant_tags);
    let clusters: Vec<Cluster> = sample_tags
        .into_par_iter()
        .flat_map_iter(|(tag, tagged_members)| {
            let mut clusters = Vec::new();
            let tag_reference = reference_tags.get(&tag);
            let tag_reference_count = tag_reference.map_or(0, HashMap::len);
            let reference_fraction =
                tag_reference_count as f64 / reference_sequences.len().max(1) as f64;

            if reference_fraction > config.max_reference_fraction {
                return Vec::new();
            }

            let sequences = tagged_members
                .iter()
                .map(|member| member.sequence.clone())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>();
            if sequences.len() < config.min_members {
                return Vec::new();
            }
            let mut union = UnionFind::new(sequences.len());
            let mut variable_residues: HashMap<String, BTreeSet<u8>> = HashMap::new();
            for member in tagged_members {
                variable_residues
                    .entry(member.sequence)
                    .or_default()
                    .insert(member.residue);
            }
            let residues = sequences
                .iter()
                .map(|sequence| {
                    variable_residues
                        .get(sequence.as_str())
                        .and_then(|values| values.iter().next())
                        .copied()
                        .unwrap_or_default()
                })
                .collect::<Vec<_>>();
            for left in 0..sequences.len() {
                for right in left + 1..sequences.len() {
                    let compatible = blosum.as_ref().is_none_or(|matrix| {
                        residues[left] == residues[right]
                            || blosum::blosum_nonnegative(residues[left], residues[right], matrix)
                    });
                    if compatible {
                        union.join(left, right);
                    }
                }
            }

            let mut components: BTreeMap<usize, BTreeSet<String>> = BTreeMap::new();
            for (index, sequence) in sequences.iter().enumerate() {
                let root = union.root(index);
                components.entry(root).or_default().insert(sequence.clone());
            }
            for component in components
                .into_values()
                .filter(|members| members.len() >= config.min_members)
            {
                // let fisher = stats::hypergeometric_tail(
                //     sample_sequences.len() + reference_sequences.len(),
                //     component.len() + reference_count,
                //     sample_sequences.len(),
                //     component.len(),
                // );
                // Only reference sequences the matrix would group with this component
                // count toward the background.
                let component_residues: BTreeSet<u8> = component
                    .iter()
                    .filter_map(|sequence| variable_residues.get(sequence))
                    .flatten()
                    .copied()
                    .collect();
                let reference_count = tag_reference.map_or(0, |reference| {
                    reference
                        .values()
                        .filter(|residue| {
                            blosum.as_ref().is_none_or(|matrix| {
                                component_residues.iter().any(|member| {
                                    member == *residue
                                        || blosum::blosum_nonnegative(*member, **residue, matrix)
                                })
                            })
                        })
                        .count()
                });
                let cont_table = ContingencyTable2x2::from_reported_counts(
                    component.len(),
                    sample_sequences.len(),
                    reference_count,
                    reference_sequences.len(),
                )
                .expect("Failed to create contingency table for global cluster");
                let (fisher, _) = log_fact_table.fishers_exact(&cont_table);
                if fisher > config.fisher_threshold {
                    continue;
                }

                clusters.push(Cluster::new(
                    tag.clone(),
                    ClusterKind::Global,
                    fisher,
                    component.into_iter().collect(),
                ));
            }
            clusters
        })
        .collect();
    //progress.finish_with_message("Global cluster evaluation complete");
    eprintln!("Evaluated {} global clusters", clusters.len());
    clusters
}

fn merged_gliph2_pattern(
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

fn merge_gliph2_candidates(candidates: Vec<Cluster>, _sample_sequences: &[String]) -> Vec<Cluster> {
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

fn add_singleton_clusters(
    mut clusters: Vec<Cluster>,
    sample_sequences: &[String],
    reference_sequences: &[String],
) -> Vec<Cluster> {
    let clustered_sequences = clusters
        .iter()
        .flat_map(|cluster| cluster.member_sequences.iter().cloned())
        .collect::<HashSet<_>>();
    let log_fact_table = LogFactorialTable::new(reference_sequences.len() + sample_sequences.len());
    let pb = progress_bar_with_message(sample_sequences.len(), "Adding singleton clusters")
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

    // for sequence in sample_sequences {
    //     if clustered_sequences.contains(sequence) {
    //         continue;
    //     }
    //     let observed = sample_sequences.iter().filter(|s| s.eq(&sequence)).count();
    //     let reference_count = reference_sequences
    //         .iter()
    //         .filter(|s| s.eq(&sequence))
    //         .count()
    //         + observed;
    //     let cont_table = ContingencyTable2x2::from_reported_counts(
    //         observed,
    //         sample_sequences.len(),
    //         reference_count,
    //         reference_sequences.len(),
    //     )
    //     .expect("Failed to create contingency table for global cluster");
    //     let (fisher, _) = log_fact_table.fishers_exact(&cont_table);
    //     // let fisher = stats::hypergeometric_tail(
    //     //     sample_sequences.len() + reference_sequences.len(),
    //     //     success,
    //     //     sample_sequences.len(),
    //     //     observed,
    //     // );
    //     pb.inc(1);
    //     clusters.push(Cluster::new(
    //         sequence.clone(),
    //         ClusterKind::Single,
    //         fisher,
    //         vec![sequence.clone()],
    //     ));
    // }

    clusters.append(&mut singletons);
    clusters
}

fn hamming(left: &str, right: &str) -> Option<usize> {
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

fn gliph1_candidates(
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

pub fn mark_hla(value: &str, allele: &str) -> String {
    let mut out = String::with_capacity(value.len() + 1);
    let mut token = String::new();
    let flush = |token: &mut String, out: &mut String| {
        let trimmed = token.trim();
        if !trimmed.is_empty() && trimmed == allele && !trimmed.ends_with('!') {
            let end = token.trim_end().len();
            out.push_str(&token[..end]);
            out.push('!');
            out.push_str(&token[end..]);
            token.clear();
        } else {
            out.push_str(token);
        }
        token.clear();
    };
    for ch in value.chars() {
        if matches!(ch, ',' | ';' | '|' | '/') {
            flush(&mut token, &mut out);
            out.push(ch);
        } else {
            token.push(ch);
        }
    }
    flush(&mut token, &mut out);
    out
}

pub fn analyze(
    input_records: &[ChainRecord],
    reference: &ReferenceData,
    algorithm: Algorithm,
    sim_depth: usize,
    hla_context: &HlaContext,
    exp_ctx: &ExpansionContext,
    blosum: Option<AminoAcidSubMatrix>,
) -> Result<Vec<OutputRow>> {
    if sim_depth == 0 {
        bail!("Simulation depth must be greater than zero");
    }
    let filter = GlobalFilterConfig::new(0.05, 2, 0.1);
    let records = filtered_sequences(input_records)?;
    let mut records_by_sequence: BTreeMap<String, Vec<ChainRecord>> = BTreeMap::new();
    for record in records.iter().cloned() {
        records_by_sequence
            .entry(record.beta.clone())
            .or_default()
            .push(record);
    }
    let sample_sequences: Vec<String> = records_by_sequence.keys().cloned().collect();
    let reference_sequences: Vec<String> = reference
        .sequences
        .iter()
        .filter(|sequence| {
            sequence.len() >= 8
                && sequence.starts_with('C')
                && sequence.ends_with('F')
                && is_amino_acid_sequence(sequence)
        })
        .cloned()
        .collect();
    if reference_sequences.len() <= sample_sequences.len() {
        bail!(
            "Reference repertoire has {} valid sequences, but the sample has {}; GLIPH requires a larger reference repertoire",
            reference_sequences.len(),
            sample_sequences.len()
        );
    }
    if reference.cdr3_length_frequencies.is_empty() || reference.v_gene_frequencies.is_empty() {
        bail!("Reference data is missing V-gene or CDR3-length frequency distributions");
    }

    let mut clusters = local_candidates(
        algorithm,
        &sample_sequences,
        &reference_sequences,
        sim_depth,
        &filter,
    );

    if algorithm == Algorithm::Gliph2 {
        let global =
            global_fisher_candidates(&sample_sequences, &reference_sequences, &filter, &blosum);
        clusters.extend(global);
        clusters = merge_gliph2_candidates(clusters, &sample_sequences);
    } else {
        let cutoff = if records.len() < 125 { 2 } else { 1 };
        clusters = gliph1_candidates(clusters, &sample_sequences, cutoff);
    }
    eprintln!(
        "Found {} {} groups after clustering",
        clusters.len(),
        algorithm
    );
    clusters = add_singleton_clusters(clusters, &sample_sequences, &reference_sequences);
    clusters.sort_unstable_by(|left, right| {
        left.fisher_score
            .total_cmp(&right.fisher_score)
            .then_with(|| left.pattern.cmp(&right.pattern))
            .then_with(|| left.kind.cmp(&right.kind))
    });

    let total_records = records.len().max(1) as f64;
    let precomputed_size_col = reference
        .cluster_size_sample_sizes
        .iter()
        .enumerate()
        .min_by_key(|&(_, sample_size)| {
            let ratio = *sample_size as f64 / total_records;
            ((1.0 - ratio).abs() * 1_000_000.0) as u64
        })
        .map(|(idx, _)| idx);

    let scoring_progress = progress_bar_with_message(clusters.len(), "Scoring clusters")?;
    let rows_by_cluster: Vec<Vec<OutputRow>> = clusters
        .par_iter()
        .enumerate()
        .map(|(cluster_index, cluster)| {
            let cluster_records: ClusterRecord = ClusterRecord::new(&records, cluster);

            if cluster_records.is_empty() {
                scoring_progress.inc(1);
                return Vec::new();
            }
            let sample_score = score::sample_scores(
                &cluster_records,
                reference,
                algorithm,
                sim_depth,
                precomputed_size_col,
                hla_context,
                exp_ctx,
            );

            let rows = cluster_records
                .clone()
                .records
                .into_iter()
                .map(|record| {
                    OutputRow::new_from_records(
                        cluster_index,
                        &cluster_records,
                        &sample_score,
                        record,
                    )
                })
                .collect();
            scoring_progress.inc(1);
            rows
        })
        .collect();
    scoring_progress.finish_and_clear();
    eprintln!(
        "Cluster scoring complete in {}",
        indicatif::HumanDuration(scoring_progress.elapsed())
    );
    let output = rows_by_cluster.into_iter().flatten().collect();
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gliph::score::expansion_score;
    use crate::gliph::types::{Algorithm, ChainRecord, Cluster, GlobalFilterConfig, ReferenceData};
    use pretty_assertions::assert_eq;
    use rstest::{fixture, rstest};

    #[test]
    fn marks_identified_hla_allele() {
        assert_eq!(
            mark_hla("DRB1*15:01, DRB1*03:01", "DRB1*15:01"),
            "DRB1*15:01!, DRB1*03:01"
        );
        assert_eq!(
            mark_hla("DPA1*01:03/DPA1*02:01", "DPA1*02:01"),
            "DPA1*01:03/DPA1*02:01!"
        );
        assert_eq!(mark_hla("DPA1*01:03", "DPA1*01:04"), "DPA1*01:03");
        assert_eq!(mark_hla("01:03", "01:03"), "01:03!");
    }

    #[test]
    fn gliph2_global_filter_removes_insignificant_candidates() {
        let sample_sequences = vec!["CASSIRSSYEQYF".to_owned(), "CASSIRSSYEQAF".to_owned()];
        let reference_sequences = vec![
            "CASSIRSSYEQVF".to_owned(),
            "CASSIRSSYEQGF".to_owned(),
            "CASSIRSSYEQDF".to_owned(),
            "CASSIRSSYEQEF".to_owned(),
            "CASSIRSSYEQNF".to_owned(),
            "CASSIRSSYEQHF".to_owned(),
            "CASSIRSSYEQKF".to_owned(),
            "CASSIRSSYEQRF".to_owned(),
        ];
        let config = GlobalFilterConfig::new(0.01, 2, 1.0);

        let clusters = global_fisher_candidates(
            &sample_sequences,
            &reference_sequences,
            &config,
            &Some(AminoAcidSubMatrix::Blosum62),
        );

        assert!(
            clusters.is_empty(),
            "common reference-supported global tags should not pass the Fisher threshold"
        );
    }

    #[test]
    fn gliph2_global_filter_removes_overly_common_reference_tags() {
        let sample_sequences = vec!["CASSIRSSYEQYF".to_owned(), "CASSIRSSYEQAF".to_owned()];
        let reference_sequences = vec![
            "CASSIRSSYEQVF".to_owned(),
            "CASSIRSSYEQGF".to_owned(),
            "CAAAAAAAAEQYF".to_owned(),
            "CAAAAAAAAEQAF".to_owned(),
        ];
        let config = GlobalFilterConfig::new(1.0, 2, 0.0);

        let clusters = global_fisher_candidates(
            &sample_sequences,
            &reference_sequences,
            &config,
            &Some(AminoAcidSubMatrix::Blosum62),
        );

        assert!(
            clusters.is_empty(),
            "tags present in the reference should be filtered when max_reference_fraction is zero"
        );
    }

    #[fixture]
    fn reference_sequences() -> Vec<String> {
        (0..300)
            .map(|index| {
                let amino_acids = b"ACDEFGHIKLMNPQTVWY";
                let first = amino_acids[index % amino_acids.len()] as char;
                let second = amino_acids[(index / amino_acids.len()) % amino_acids.len()] as char;
                format!("CAAAAA{first}{second}F")
            })
            .collect()
    }

    #[rstest]
    fn gliph2_global_filter_keeps_significant_rare_candidates(reference_sequences: Vec<String>) {
        let sample_sequences = vec![
            "CASSIRSSYEQYF".to_owned(),
            "CASSIRSSYEQAF".to_owned(),
            "CASSIRSSYEQVF".to_owned(),
        ];

        let config = GlobalFilterConfig::new(0.01, 2, 0.01);

        let clusters = global_fisher_candidates(
            &sample_sequences,
            &reference_sequences,
            &config,
            &Some(AminoAcidSubMatrix::Blosum62),
        );

        assert!(
            clusters
                .iter()
                .any(|cluster| cluster.member_sequences.len() == 3),
            "rare enriched sample-only global tags should pass filtering"
        );
        assert!(
            clusters
                .iter()
                .all(|cluster| cluster.fisher_score <= config.fisher_threshold)
        );
    }
    #[test]
    fn gliph2_merge_uses_minimal_local_pattern_label() {
        let sample_sequences = vec![
            "CASSQFLAAQTQYF".to_owned(),
            "CASSQFLGGDTQYF".to_owned(),
            "CASSQFLSSYEQYF".to_owned(),
        ];
        let candidates = vec![Cluster {
            pattern: "SQFL".to_owned(),
            kind: ClusterKind::Local,
            fisher_score: 0.001,
            member_sequences: sample_sequences.clone(),
        }];

        let merged = merge_gliph2_candidates(candidates, &sample_sequences);

        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].pattern, "SQFL");
    }

    #[test]
    fn gliph2_merge_prefers_longest_local_pattern_label() {
        let sample_sequences = vec![
            "CASSQFLAAQTQYF".to_owned(),
            "CASSQFLGGDTQYF".to_owned(),
            "CASSQFLSSYEQYF".to_owned(),
        ];
        let candidates = vec![
            Cluster {
                pattern: "QF".to_owned(),
                kind: ClusterKind::Local,
                fisher_score: 0.01,
                member_sequences: sample_sequences.clone(),
            },
            Cluster {
                pattern: "SQFL".to_owned(),
                kind: ClusterKind::Local,
                fisher_score: 0.001,
                member_sequences: sample_sequences.clone(),
            },
        ];

        let merged = merge_gliph2_candidates(candidates, &sample_sequences);

        assert_eq!(merged.len(), 2);
        assert_eq!(merged[0].pattern, "QF");
        assert_eq!(merged[1].pattern, "SQFL");
        assert_eq!(merged[0].fisher_score, 0.01);
        assert_eq!(merged[1].fisher_score, 0.001);
    }

    #[test]
    fn gliph2_merge_uses_trimmed_global_pattern_when_no_local_pattern_exists() {
        let sample_sequences = vec!["CASSQFLAAQTQYF".to_owned(), "CASSQFLGAQTQYF".to_owned()];
        let candidates = vec![Cluster {
            pattern: "SQFL%AQT".to_owned(),
            kind: ClusterKind::Global,
            fisher_score: 0.001,
            member_sequences: sample_sequences.clone(),
        }];

        let merged = merge_gliph2_candidates(candidates, &sample_sequences);

        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].pattern, "SQFL%AQT");
    }

    fn map_test_sample_to_chain_record(index: usize, beta: &str) -> ChainRecord {
        ChainRecord {
            subject: format!("S{}", index + 1),
            condition: "case".to_owned(),
            beta: (*beta).to_owned(),
            alpha: String::new(),
            v_gene: "TRBV1".to_owned(),
            j_gene: "TRBJ1".to_owned(),
            counts: None,
            hla: [
                "DPA1*01:03".to_owned(),
                "DPB1*04:01".to_owned(),
                "DQA1*05:01".to_owned(),
                "DQB1*02:01".to_owned(),
                "DRB1*15:01".to_owned(),
                String::new(),
                String::new(),
                String::new(),
            ],
        }
    }

    #[test]
    fn gliph2_builds_annotated_cluster_rows() {
        let sample = ["CASSIRSSYEQYF", "CASSIRSSYEQAF", "CASSIRSSYEQVF"];
        let _ = sample
            .iter()
            .enumerate()
            .map(|(index, beta)| map_test_sample_to_chain_record(index, beta))
            .collect::<Vec<_>>();
    }
    #[rstest]
    fn gliph2_global_filter_honors_min_members(reference_sequences: Vec<String>) {
        let sample_sequences = vec![
            "CASSIRSSYEQYF".to_owned(),
            "CASSIRSSYEQAF".to_owned(),
            "CASSIRSSYEQVF".to_owned(),
        ];

        let config = GlobalFilterConfig::new(1.0, 4, 1.0);

        let clusters = global_fisher_candidates(
            &sample_sequences,
            &reference_sequences,
            &config,
            &Some(AminoAcidSubMatrix::Blosum62),
        );

        assert!(
            clusters.is_empty(),
            "global candidates smaller than min_members should be filtered"
        );
    }

    #[test]
    fn gliph2_builds_annotated_cluster_rows2() {
        let sample = ["CASSIRSSYEQYF", "CASSIRSSYEQAF", "CASSIRSSYEQVF"];
        let records = sample
            .iter()
            .enumerate()
            .map(|(index, beta)| map_test_sample_to_chain_record(index, beta))
            .collect::<Vec<_>>();
        let hla_context = HlaContext::build(&records, 0.05);
        let exp_ctx = ExpansionContext::build(&records);
        let alphabet = b"ACDEFGHIKLMNPQTVWY";
        let reference_sequences = (0..30)
            .map(|index| {
                let first = alphabet[index % alphabet.len()] as char;
                let second = alphabet[(index / alphabet.len()) % alphabet.len()] as char;
                format!("CAAAA{first}{second}F")
            })
            .collect::<Vec<_>>();
        let reference = ReferenceData {
            sequences: reference_sequences,
            v_gene_frequencies: vec![("TRBV1".to_owned(), 1.0)],
            cdr3_length_frequencies: vec![(13, 1.0)],
            cluster_size_probabilities: Vec::new(),
            cluster_size_sample_sizes: Vec::new(),
        };

        let result = analyze(
            &records,
            &reference,
            Algorithm::Gliph2,
            20,
            &hla_context,
            &exp_ctx,
            Some(AminoAcidSubMatrix::Blosum62),
        )
        .unwrap();

        assert!(!result.is_empty());
        assert!(result.iter().all(|row| row.number_subject == 3));
        assert!(result.iter().all(|row| row.number_unique_cdr3 >= 2));
        assert!(result.iter().all(|row| !row.hla_dpa1.is_empty()));
        assert_eq!(
            result,
            analyze(
                &records,
                &reference,
                Algorithm::Gliph2,
                20,
                &hla_context,
                &exp_ctx,
                Some(AminoAcidSubMatrix::Blosum62)
            )
            .unwrap()
        );

        let gliph1_result = analyze(
            &records,
            &reference,
            Algorithm::Gliph1,
            20,
            &hla_context,
            &exp_ctx,
            Some(AminoAcidSubMatrix::Blosum62),
        )
        .unwrap();
        assert!(!gliph1_result.is_empty());
        assert!(gliph1_result.iter().all(|row| row.kind == "global"));
    }

    #[test]
    fn expansion_score_is_neutral_without_counts() {
        let records = vec![
            ChainRecord {
                beta: "CASSAAAAEQYF".to_owned(),
                ..ChainRecord::default()
            },
            ChainRecord {
                beta: "CASSBBBBEQYF".to_owned(),
                ..ChainRecord::default()
            },
        ];
        let members: HashSet<&str> = ["CASSAAAAEQYF"].into_iter().collect::<HashSet<_>>();
        let exp_ctx = ExpansionContext::build(&records);

        let score = expansion_score(&members, &exp_ctx, 100, "no-counts");

        assert_eq!(score, 1.0);
    }

    #[test]
    fn expansion_score_aggregates_counts_by_unique_cdr3() {
        let records = vec![
            ChainRecord {
                beta: "CASSAAAAEQYF".to_owned(),
                counts: Some(10.0),
                ..ChainRecord::default()
            },
            ChainRecord {
                beta: "CASSAAAAEQYF".to_owned(),
                counts: Some(15.0),
                ..ChainRecord::default()
            },
            ChainRecord {
                beta: "CASSBBBBEQYF".to_owned(),
                counts: Some(1.0),
                ..ChainRecord::default()
            },
            ChainRecord {
                beta: "CASSCCCCEQYF".to_owned(),
                counts: Some(1.0),
                ..ChainRecord::default()
            },
        ];
        let members: HashSet<&str> = ["CASSAAAAEQYF"].into_iter().collect::<HashSet<_>>();
        let exp_ctx = ExpansionContext::build(&records);

        let score = expansion_score(&members, &exp_ctx, 100, "aggregated-counts");

        assert!(
            score < 1.0,
            "expanded unique CDR3 should be more extreme than at least some background draws"
        );
    }

    #[test]
    fn expansion_score_is_deterministic() {
        let records = vec![
            ChainRecord {
                beta: "CASSAAAAEQYF".to_owned(),
                counts: Some(20.0),
                ..ChainRecord::default()
            },
            ChainRecord {
                beta: "CASSBBBBEQYF".to_owned(),
                counts: Some(1.0),
                ..ChainRecord::default()
            },
            ChainRecord {
                beta: "CASSCCCCEQYF".to_owned(),
                counts: Some(1.0),
                ..ChainRecord::default()
            },
            ChainRecord {
                beta: "CASSDDDDEQYF".to_owned(),
                counts: Some(1.0),
                ..ChainRecord::default()
            },
        ];
        let members: HashSet<&str> = ["CASSAAAAEQYF"].into_iter().collect::<HashSet<_>>();
        let exp_ctx = ExpansionContext::build(&records);

        assert_eq!(
            expansion_score(&members, &exp_ctx, 100, "deterministic"),
            expansion_score(&members, &exp_ctx, 100, "deterministic")
        );
    }
}
