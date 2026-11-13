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
    use super::*;
    use crate::tcr_cluster::References;
    use rdstream::{Document, Encoding, Format, RObject, Value};
    use std::path::Path;

    fn chr(values: &[&str]) -> RObject {
        RObject::new(Value::Character(
            values.iter().map(|v| Some((*v).to_owned())).collect(),
        ))
    }

    fn real(values: &[f64]) -> RObject {
        RObject::new(Value::Real(values.to_vec()))
    }

    fn named_list(items: Vec<(&str, RObject)>) -> RObject {
        let names = chr(&items.iter().map(|(n, _)| *n).collect::<Vec<_>>());
        let mut list = RObject::new(Value::List(items.into_iter().map(|(_, o)| o).collect()));
        list.attributes.push(("names".to_owned(), names));
        list
    }

    fn document(objects: Vec<(&str, RObject)>) -> Document {
        Document {
            format: Format::RData,
            encoding: Encoding::Xdr,
            version: 3,
            writer_version: (4, 4, 0),
            native_encoding: None,
            root: RObject::new(Value::Pairlist {
                sexptype: 2,
                items: objects
                    .into_iter()
                    .map(|(n, o)| (Some(n.to_owned()), o))
                    .collect(),
            }),
            environments: Vec::new(),
            bytes_read: 0,
        }
    }

    fn reference_doc(refseqs: RObject, vgenes: RObject, lengths: RObject) -> Document {
        let selected = named_list(vec![
            ("refseqs", refseqs),
            ("vgene_frequencies", vgenes),
            ("cdr3_length_frequencies", lengths),
        ]);
        document(vec![(
            "reference_list",
            named_list(vec![("human_v2.0_CD48", selected)]),
        )])
    }

    fn valid_refseqs() -> RObject {
        named_list(vec![
            ("CDR3b", chr(&["CASSF", "CASRF"])),
            ("TRBV", chr(&["TRBV2", "TRBV28"])),
        ])
    }

    fn valid_vgenes() -> RObject {
        named_list(vec![
            ("vgene", chr(&["TRBV2", "TRBV28"])),
            ("freq", real(&[0.25, 0.75])),
        ])
    }

    fn valid_lengths() -> RObject {
        named_list(vec![
            ("len", RObject::new(Value::Integer(vec![12, 13]))),
            ("freq", real(&[0.4, 0.6])),
        ])
    }

    fn sizes_doc(original: RObject) -> Document {
        document(vec![(
            "ref_cluster_sizes",
            named_list(vec![("original", original)]),
        )])
    }

    #[test]
    fn reference_data_reads_selected_reference() {
        let doc = reference_doc(valid_refseqs(), valid_vgenes(), valid_lengths());
        let data = reference_data(&doc, &References::HumanV2CD4CD8).unwrap();
        assert_eq!(data.sequences, vec!["CASSF", "CASRF"]);
        assert_eq!(
            data.v_gene_frequencies,
            vec![("TRBV2".to_owned(), 0.25), ("TRBV28".to_owned(), 0.75)]
        );
        assert_eq!(data.cdr3_length_frequencies, vec![(12, 0.4), (13, 0.6)]);
    }

    #[test]
    fn reference_data_reports_missing_pieces() {
        let doc = reference_doc(valid_refseqs(), valid_vgenes(), valid_lengths());
        let err = reference_data(&doc, &References::MouseV1CD4).unwrap_err();
        assert!(format!("{err:#}").contains("mouse_v1.0_CD4"));

        let err = reference_data(&document(vec![]), &References::HumanV2CD4CD8).unwrap_err();
        assert!(format!("{err:#}").contains("reference_list"));
    }

    #[test]
    fn reference_data_rejects_mismatched_columns() {
        let short_v = named_list(vec![
            ("CDR3b", chr(&["CASSF", "CASRF"])),
            ("TRBV", chr(&["TRBV2"])),
        ]);
        let doc = reference_doc(short_v, valid_vgenes(), valid_lengths());
        assert!(reference_data(&doc, &References::HumanV2CD4CD8).is_err());

        let short_freq = named_list(vec![
            ("vgene", chr(&["TRBV2", "TRBV28"])),
            ("freq", real(&[1.0])),
        ]);
        let doc = reference_doc(valid_refseqs(), short_freq, valid_lengths());
        assert!(reference_data(&doc, &References::HumanV2CD4CD8).is_err());
    }

    #[test]
    fn vector_helpers_reject_wrong_types_and_na() {
        assert!(numeric_vector(&chr(&["x"]), "ctx").is_err());
        assert!(character_vector(&real(&[1.0]), "ctx").is_err());
        let with_na = RObject::new(Value::Character(vec![Some("a".into()), None]));
        assert!(character_vector(&with_na, "ctx").is_err());
    }

    #[test]
    fn cluster_size_table_drops_size_column_and_parses_names() {
        let original = named_list(vec![
            ("cluster.size", real(&[1.0, 2.0])),
            ("125", real(&[0.9, 0.1])),
            ("250", chr(&["0.8", "0.2"])),
            ("label", RObject::new(Value::Logical(vec![1, 0]))),
        ]);
        let (matrix, sizes) = cluster_size_table(&sizes_doc(original)).unwrap();
        assert_eq!(sizes, vec![125, 250]);
        assert_eq!(matrix, vec![vec![0.9, 0.8], vec![0.1, 0.2]]);
    }

    #[test]
    fn cluster_size_table_falls_back_to_default_sizes() {
        let original = named_list(vec![("a", real(&[0.5])), ("b", real(&[0.5]))]);
        let (_, sizes) = cluster_size_table(&sizes_doc(original)).unwrap();
        assert_eq!(sizes, vec![125, 250]);
    }

    #[test]
    fn cluster_size_table_rejects_malformed_tables() {
        let bad_cases = vec![
            // not a list
            real(&[1.0]),
            // fewer than two numeric columns
            named_list(vec![("125", real(&[1.0]))]),
            // ragged columns
            named_list(vec![("125", real(&[1.0, 2.0])), ("250", real(&[1.0]))]),
            // unparsable character cell
            named_list(vec![("125", real(&[1.0])), ("250", chr(&["oops"]))]),
            // NA character cell
            named_list(vec![
                ("125", real(&[1.0])),
                ("250", RObject::new(Value::Character(vec![None]))),
            ]),
        ];
        for original in bad_cases {
            assert!(cluster_size_table(&sizes_doc(original)).is_err());
        }

        let mut unnamed = named_list(vec![("125", real(&[1.0])), ("250", real(&[1.0]))]);
        unnamed.attributes[0].1 = RObject::new(Value::Character(vec![Some("125".into()), None]));
        assert!(cluster_size_table(&sizes_doc(unnamed)).is_err());

        let mut nameless = named_list(vec![("125", real(&[1.0])), ("250", real(&[1.0]))]);
        nameless.attributes.clear();
        assert!(cluster_size_table(&sizes_doc(nameless)).is_err());

        let mut short_names = named_list(vec![("125", real(&[1.0])), ("250", real(&[1.0]))]);
        short_names.attributes[0].1 = chr(&["125"]);
        assert!(cluster_size_table(&sizes_doc(short_names)).is_err());
    }

    #[test]
    fn loads_reference_from_cached_fixture_files() {
        let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/gliph_reference");
        let cache = tempfile::tempdir().unwrap();
        for name in ["reference_list.RData", "ref_cluster_sizes.RData"] {
            std::fs::copy(fixtures.join(name), cache.path().join(name)).unwrap();
        }
        let reference =
            load_reference(&References::HumanV2CD4CD8, Some(cache.path()), false).unwrap();
        assert_eq!(reference.sequences.len(), 2000);
        assert!(!reference.v_gene_frequencies.is_empty());
        assert_eq!(reference.cluster_size_sample_sizes.len(), 9);
    }

    #[test]
    fn load_rdata_reports_unparsable_files() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("broken.RData");
        std::fs::write(&path, b"not an rdata file").unwrap();
        assert!(load_rdata(&path).is_err());
    }

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
