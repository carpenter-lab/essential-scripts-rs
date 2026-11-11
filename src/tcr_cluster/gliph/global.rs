use super::*;
use crate::tcr_cluster::blosum::Chain;

pub(super) fn global_tags(sequences: &[String]) -> BTreeMap<String, Vec<TaggedMember>> {
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

pub(super) fn reference_tag_counts(
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

pub(super) fn global_fisher_candidates(
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
                            || blosum::blosum_nonnegative(
                                residues[left],
                                residues[right],
                                matrix,
                                &Chain::Beta,
                            )
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
                                        || blosum::blosum_nonnegative(
                                            *member,
                                            **residue,
                                            matrix,
                                            &Chain::Beta,
                                        )
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
