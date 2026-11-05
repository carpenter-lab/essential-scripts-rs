use crate::cache;
use crate::tcr_cluster::References;
use crate::tcr_cluster::model::ReferenceData;
use anyhow::{Context, Result, bail};
use rdstream::{Document, RObject, Rules, Value};
use std::path::Path;

const REFERENCE_URL: &str =
    "https://zenodo.org/records/18925758/files/reference_list.RData?download=1";
const REFERENCE_MD5: &str = "73cef8be5b9016d88343e9e4c38d85a8";
const REFERENCE_SIZE: usize = 19_001_532;
const CLUSTER_SIZES_URL: &str = "https://raw.githubusercontent.com/BorchLab/immGLIPH/eaf5485995b795744fc80580a11b0cf25233bcaa/data/ref_cluster_sizes.RData";

fn object<'a>(document: &'a Document, name: &str) -> Result<&'a RObject> {
    document
        .objects()
        .into_iter()
        .find_map(|(object_name, object)| (object_name == name).then_some(object))
        .with_context(|| format!("RData file does not contain object '{name}'"))
}

fn character_vector(object: &RObject, context: &str) -> Result<Vec<String>> {
    match &object.value {
        Value::Character(values) => values
            .iter()
            .map(|value| {
                value
                    .clone()
                    .with_context(|| format!("Missing value in {context}"))
            })
            .collect(),
        _ => bail!("Expected a character vector for {context}"),
    }
}

fn numeric_vector(object: &RObject, context: &str) -> Result<Vec<f64>> {
    match &object.value {
        Value::Real(values) => Ok(values.clone()),
        Value::Integer(values) => Ok(values.iter().map(|value| *value as f64).collect()),
        _ => bail!("Expected a numeric vector for {context}"),
    }
}

fn named_component<'a>(object: &'a RObject, name: &str, context: &str) -> Result<&'a RObject> {
    object
        .get(name)
        .with_context(|| format!("Missing '{name}' component in {context}"))
}

fn reference_data(document: &Document, reference: &References) -> Result<ReferenceData> {
    let reference_string = reference.to_string();
    let reference_name = reference_string.as_str();
    let all_references = object(document, "reference_list")?;
    let selected = named_component(all_references, reference_name, "reference_list")?;
    let seq_table = named_component(selected, "refseqs", reference_name)?;
    let sequences = character_vector(
        named_component(seq_table, "CDR3b", "refseqs")?,
        "refseqs$CDR3b",
    )?;
    let v_genes = character_vector(
        named_component(seq_table, "TRBV", "refseqs")?,
        "refseqs$TRBV",
    )?;
    if sequences.len() != v_genes.len() {
        bail!("Reference sequence and TRBV vectors have different lengths");
    }

    let v_table = named_component(selected, "vgene_frequencies", reference_name)?;
    let v_frequency_names = character_vector(
        named_component(v_table, "vgene", "vgene_frequencies")?,
        "vgene_frequencies$vgene",
    )?;
    let v_frequencies = numeric_vector(
        named_component(v_table, "freq", "vgene_frequencies")?,
        "vgene_frequencies$freq",
    )?;
    let cdr3_table = named_component(selected, "cdr3_length_frequencies", reference_name)?;
    let lengths = numeric_vector(
        named_component(cdr3_table, "len", "cdr3_length_frequencies")?,
        "cdr3_length_frequencies$len",
    )?;
    let length_frequencies = numeric_vector(
        named_component(cdr3_table, "freq", "cdr3_length_frequencies")?,
        "cdr3_length_frequencies$freq",
    )?;

    if v_frequency_names.len() != v_frequencies.len() || lengths.len() != length_frequencies.len() {
        bail!("Reference frequency tables contain mismatched columns");
    }

    Ok(ReferenceData {
        sequences,
        v_gene_frequencies: v_frequency_names.into_iter().zip(v_frequencies).collect(),
        cdr3_length_frequencies: lengths
            .into_iter()
            .zip(length_frequencies)
            .map(|(length, frequency)| (length as usize, frequency))
            .collect(),
        cluster_size_probabilities: Vec::new(),
        cluster_size_sample_sizes: Vec::new(),
    })
}

fn cluster_size_table(document: &Document) -> Result<(Vec<Vec<f64>>, Vec<usize>)> {
    let sizes = object(document, "ref_cluster_sizes")?;
    let original = named_component(sizes, "original", "ref_cluster_sizes")?;
    let Value::List(columns) = &original.value else {
        bail!("ref_cluster_sizes$original is not a data frame/list");
    };
    let names = original
        .names()
        .context("ref_cluster_sizes$original is missing column names")?;
    if names.len() != columns.len() {
        bail!("Cluster-size table column names do not match its columns");
    }
    let mut numeric_columns = Vec::new();
    for (name, column) in names.iter().zip(columns) {
        let Some(name) = name else {
            bail!("Cluster-size table contains an unnamed column");
        };
        let values = match &column.value {
            Value::Real(_) | Value::Integer(_) => numeric_vector(column, name)?,
            Value::Character(values) => values
                .iter()
                .map(|value| {
                    value
                        .as_deref()
                        .context("Missing cluster-size table cell")?
                        .parse::<f64>()
                        .with_context(|| format!("Invalid cluster-size probability in {name}"))
                })
                .collect::<Result<Vec<_>>>()?,
            _ => continue,
        };
        numeric_columns.push((name.as_str(), values));
    }
    if numeric_columns.len() < 2 {
        bail!("Cluster-size probability table has fewer than two numeric columns");
    }
    let first_name = numeric_columns[0].0.to_ascii_lowercase();
    if first_name.contains("size") || first_name == "x" || first_name == "v1" {
        numeric_columns.remove(0);
    }
    let row_count = numeric_columns
        .first()
        .context("Cluster-size probability table contains no score columns")?
        .1
        .len();
    if numeric_columns
        .iter()
        .any(|(_, values)| values.len() != row_count)
    {
        bail!("Cluster-size probability columns have different lengths");
    }
    let sizes = numeric_columns
        .iter()
        .enumerate()
        .map(|(index, (name, _))| {
            let parsed = name
                .chars()
                .filter(char::is_ascii_digit)
                .collect::<String>()
                .parse::<usize>()
                .ok();
            parsed.unwrap_or_else(|| {
                [125, 250, 500, 1000, 2000, 4000, 6000, 8000, 10000]
                    .get(index)
                    .copied()
                    .unwrap_or(index + 1)
            })
        })
        .collect::<Vec<_>>();
    let mut matrix = vec![vec![0.0; numeric_columns.len()]; row_count];
    for (column_index, (_, values)) in numeric_columns.into_iter().enumerate() {
        for (row_index, value) in values.into_iter().enumerate() {
            matrix[row_index][column_index] = value;
        }
    }
    Ok((matrix, sizes))
}

fn load_rdata(path: &Path) -> Result<rdstream::Parsed> {
    let mut rules = Rules::all();
    rdstream::parse_file(path, &mut rules)
        .map_err(anyhow::Error::new)
        .with_context(|| format!("Failed to parse RData file {}", path.display()))
}

pub fn load_reference(
    reference_name: &References,
    cache_dir: Option<&Path>,
    force_download: bool,
) -> Result<ReferenceData> {
    let reference_path = cache::cached_file(
        cache_dir,
        "reference_list.RData",
        REFERENCE_URL,
        force_download,
        true,
        REFERENCE_SIZE,
        REFERENCE_MD5,
    )?;
    let cluster_path = cache::cached_file(
        cache_dir,
        "ref_cluster_sizes.RData",
        CLUSTER_SIZES_URL,
        force_download,
        false,
        REFERENCE_SIZE,
        REFERENCE_MD5,
    )?;

    let mut reference = reference_data(&load_rdata(&reference_path)?.document, reference_name)?;
    (
        reference.cluster_size_probabilities,
        reference.cluster_size_sample_sizes,
    ) = cluster_size_table(&load_rdata(&cluster_path)?.document)?;
    Ok(reference)
}

#[cfg(test)]
mod tests {
    use super::{load_reference, numeric_vector};
    use crate::tcr_cluster::References;
    use rdstream::{RObject, Value};

    #[test]
    fn numeric_vector_accepts_integer_and_real_values() {
        let integers = RObject::new(Value::Integer(vec![1, 2]));
        let reals = RObject::new(Value::Real(vec![1.5, 2.5]));
        assert_eq!(numeric_vector(&integers, "ints").unwrap(), vec![1.0, 2.0]);
        assert_eq!(numeric_vector(&reals, "reals").unwrap(), vec![1.5, 2.5]);
    }

    #[test]
    #[ignore = "downloads published immGLIPH reference files"]
    fn loads_published_reference_files() {
        let cache = tempfile::tempdir().unwrap();
        let reference =
            load_reference(&References::HumanV2CD4CD8, Some(cache.path()), true).unwrap();
        assert!(!reference.sequences.is_empty());
        assert!(!reference.v_gene_frequencies.is_empty());
        assert!(!reference.cdr3_length_frequencies.is_empty());
        assert!(!reference.cluster_size_probabilities.is_empty());
        assert!(!reference.cluster_size_sample_sizes.is_empty());
    }
}
