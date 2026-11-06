use crate::cache;
use crate::tcr_cluster::model::{ChainRecord, Organism, with_default_allele};
use anyhow::{Context, Result, bail};
use std::borrow::Cow;
use std::collections::HashMap;
use std::path::Path;

/// tcrdist3's gene database, pinned to the last commit that changed it.
const TCRDIST_DB_URL: &str = "https://raw.githubusercontent.com/kmayerb/tcrdist3/f56bd4aae07a437d180ffb795392b72031b0fa1d/tcrdist/db/alphabeta_gammadelta_db.tsv";
const TCRDIST_DB_FILENAME: &str = "alphabeta_gammadelta_db.tsv";
const TCRDIST_DB_MD5: &str = "b227df55fd0749fd53621cdaa8b2d600";
const TCRDIST_DB_SIZE: usize = 487_823;

/// tcrdist3 separates the CDRs within the `cdrs` column with this character.
const CDRS_SEP: char = ';';
/// IMGT alignment gap character used in the `cdrs` column.
pub(super) const GAP: char = '.';

/// Germline CDRs encoded by a V allele (tcrdist3 `TCR_Gene.cdrs[0..3]`).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct VGeneCdrs {
    pub cdr1: String,
    pub cdr2: String,
    /// CDR2.5, which tcrdist3 calls `pmhc`.
    pub pmhc: String,
}

impl VGeneCdrs {
    pub fn without_gaps(&self) -> Self {
        let strip = |s: &str| s.replace(GAP, "");
        Self {
            cdr1: strip(&self.cdr1),
            cdr2: strip(&self.cdr2),
            pmhc: strip(&self.pmhc),
        }
    }
}

/// V-region rows of tcrdist3's `alphabeta_gammadelta_db.tsv` for one organism,
/// keyed by allele id (e.g. `TRBV19*01`).
pub struct TcrDistDb {
    organism: Organism,
    genes: HashMap<String, VGeneCdrs>,
}

impl TcrDistDb {
    /// Load the database from the reference cache, downloading it first if it is
    /// missing or `force_download` is set.
    pub fn load(
        organism: Organism,
        cache_dir: Option<&Path>,
        force_download: bool,
    ) -> Result<Self> {
        let path = cache::cached_file(
            cache_dir,
            TCRDIST_DB_FILENAME,
            TCRDIST_DB_URL,
            force_download,
            true,
            TCRDIST_DB_SIZE,
            TCRDIST_DB_MD5,
        )?;
        let tsv = std::fs::read_to_string(&path)
            .with_context(|| format!("Failed to read tcrdist3 gene db {}", path.display()))?;
        Self::from_tsv(&tsv, organism)
    }

    pub fn from_tsv(tsv: &str, organism: Organism) -> Result<Self> {
        let mut lines = tsv.lines();
        let header: Vec<&str> = lines
            .next()
            .context("tcrdist3 gene db is empty")?
            .split('\t')
            .collect();
        let column = |name: &str| {
            header
                .iter()
                .position(|h| *h == name)
                .with_context(|| format!("tcrdist3 gene db is missing column '{name}'"))
        };
        let id_i = column("id")?;
        let organism_i = column("organism")?;
        let region_i = column("region")?;
        let cdrs_i = column("cdrs")?;

        let organism_name = organism.to_string();
        let mut genes = HashMap::new();
        for line in lines.filter(|l| !l.is_empty()) {
            let fields: Vec<&str> = line.split('\t').collect();
            let field = |i: usize| fields.get(i).copied().unwrap_or("");
            if field(organism_i) != organism_name || field(region_i) != "V" {
                continue;
            }
            // [CDR1, CDR2, CDR2.5, partial CDR3]; pseudogenes have no CDRs and are skipped,
            // matching tcrdist3, where they end up with an empty `cdrs` list.
            let cdrs: Vec<&str> = field(cdrs_i).split(CDRS_SEP).collect();
            let [cdr1, cdr2, pmhc, ..] = cdrs[..] else {
                continue;
            };
            if cdr1.is_empty() {
                continue;
            }
            genes.insert(
                field(id_i).to_owned(),
                VGeneCdrs {
                    cdr1: cdr1.to_owned(),
                    cdr2: cdr2.to_owned(),
                    pmhc: pmhc.to_owned(),
                },
            );
        }
        if genes.is_empty() {
            bail!("tcrdist3 gene db has no V genes for organism '{organism}'");
        }
        Ok(Self { organism, genes })
    }

    pub fn organism(&self) -> Organism {
        self.organism
    }

    /// Look up a V allele. tcrdist3 requires the allele; a name without one
    /// (`TRBV19`) is looked up as `*01` (`TRBV19*01`).
    fn get(&self, v_gene: &str) -> Option<&VGeneCdrs> {
        self.genes.get(with_default_allele(v_gene).as_ref())
    }
}

/// Germline CDRs for both chains of a record. Used internally for distances only;
/// they are never stored on the record or written to the output.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ChainCdrs {
    pub beta: Option<VGeneCdrs>,
    /// `None` when the input had no `TRAV` column (beta-only run) or the gene is unknown.
    pub alpha: Option<VGeneCdrs>,
}

/// Infer the CDR1/CDR2/CDR2.5 of a record's beta and (if given) alpha V genes.
pub fn infer_record_cdrs(db: &TcrDistDb, record: &ChainRecord, imgt_aligned: bool) -> ChainCdrs {
    let infer = |gene: Option<Cow<'_, str>>| {
        gene.and_then(|gene| infer_cdr_seqs_from_v_gene(db, &gene, imgt_aligned))
    };
    ChainCdrs {
        beta: infer(record.vb_allele()),
        alpha: infer(record.va_allele()),
    }
}

/// Port of tcrdist3's `TCRrep.infer_cdrs_from_v_gene` for a single V gene.
/// With `imgt_aligned = false` the IMGT gap characters are removed.
fn infer_cdr_seqs_from_v_gene(
    db: &TcrDistDb,
    v_gene: &str,
    imgt_aligned: bool,
) -> Option<VGeneCdrs> {
    let cdrs = db.get(v_gene)?;
    Some(if imgt_aligned {
        cdrs.clone()
    } else {
        cdrs.without_gaps()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "id\torganism\tchain\tregion\tnucseq\tframe\taligned_protseq\tcdr_columns\tcdrs
TRAV1*01\tmouse\tA\tV\tggg\t+1\tGQG\t28-39;57-66;82-88;106-111\tTSG......FNG;VVL....DGL;SRSN.GY;CAVR..
TRBV19*01\thuman\tB\tV\tggg\t+1\tDGG\t27-38;56-65;81-86;104-111\tMN.......HEY;SVG....AGI;TDQGEV;CASSI...
TRBV1*01\thuman\tB\tV\tggg\t+1\tDGG\t\t
TRBJ1-1*01\thuman\tB\tJ\tggg\t+1\tNTE\t\t
";

    #[test]
    fn parses_v_genes_for_requested_organism() {
        let db = TcrDistDb::from_tsv(SAMPLE, Organism::Human).unwrap();
        assert_eq!(db.genes.len(), 1, "mouse, J and CDR-less rows are skipped");
        let cdrs = infer_cdr_seqs_from_v_gene(&db, "TRBV19*01", true).unwrap();
        assert_eq!(cdrs.cdr1, "MN.......HEY");
        assert_eq!(cdrs.cdr2, "SVG....AGI");
        assert_eq!(cdrs.pmhc, "TDQGEV");
    }

    #[test]
    fn strips_gaps_and_falls_back_to_first_allele() {
        let db = TcrDistDb::from_tsv(SAMPLE, Organism::Human).unwrap();
        let cdrs = infer_cdr_seqs_from_v_gene(&db, "TRBV19", false).unwrap();
        assert_eq!(cdrs.cdr1, "MNHEY");
        assert!(infer_cdr_seqs_from_v_gene(&db, "TRBV19*02", true).is_none());
    }

    #[test]
    fn errors_when_organism_has_no_v_genes() {
        let only_mouse = SAMPLE.lines().take(2).collect::<Vec<_>>().join("\n");
        assert!(TcrDistDb::from_tsv(&only_mouse, Organism::Human).is_err());
    }

    #[test]
    #[ignore = "downloads the tcrdist3 gene database"]
    fn loads_published_tcrdist_db() {
        let cache = tempfile::tempdir().unwrap();
        let db = TcrDistDb::load(Organism::Human, Some(cache.path()), true).unwrap();
        let cdrs = db.get("TRBV19*01").unwrap();
        assert!(!cdrs.cdr1.is_empty() && !cdrs.pmhc.is_empty());
    }

    #[test]
    fn infers_alpha_and_beta_cdrs_from_record_genes_without_alleles() {
        let tsv = format!(
            "{SAMPLE}TRAV12-2*01\thuman\tA\tV\tggg\t+1\tKEV\t27-38;56-65;81-86;104-111\tDRG......SQS;IYS....NGD;NKGSN;CAVN...\n"
        );
        let db = TcrDistDb::from_tsv(&tsv, Organism::Human).unwrap();
        let record = ChainRecord {
            v_gene: Some("TRBV19".to_owned()),
            va_gene: Some("TRAV12-2".to_owned()),
            ..ChainRecord::default()
        };
        let cdrs = infer_record_cdrs(&db, &record, true);
        assert_eq!(cdrs.beta.unwrap().pmhc, "TDQGEV");
        assert_eq!(cdrs.alpha.unwrap().cdr1, "DRG......SQS");
        // The record itself keeps the user's gene names.
        assert_eq!(record.v_gene.as_deref(), Some("TRBV19"));
    }

    #[test]
    fn beta_only_record_has_no_alpha_cdrs() {
        let db = TcrDistDb::from_tsv(SAMPLE, Organism::Human).unwrap();
        let record = ChainRecord {
            v_gene: Some("TRBV19*01".to_owned()),
            ..ChainRecord::default()
        };
        let cdrs = infer_record_cdrs(&db, &record, false);
        assert_eq!(cdrs.beta.unwrap().cdr1, "MNHEY");
        assert!(cdrs.alpha.is_none());
    }
}
