use super::chain::ChainRecord;
use super::output::OutputRow;
use super::score::ClusterScore;
use std::collections::HashSet;
use std::fmt::{Display, Formatter};
use std::hash::Hash;
use std::sync::Arc;

#[derive(Clone, Debug, PartialEq, PartialOrd, Ord, Eq, Copy, Hash)]
pub enum ClusterKind {
    Local,
    Global,
    Single,
}

impl Display for ClusterKind {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            ClusterKind::Local => write!(f, "local"),
            ClusterKind::Global => write!(f, "global"),
            ClusterKind::Single => write!(f, "single"),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Cluster {
    pub pattern: String,
    pub kind: ClusterKind,
    pub fisher_score: f64,
    pub member_sequences: Vec<String>,
}

impl Cluster {
    pub fn new(
        pattern: String,
        kind: ClusterKind,
        fisher_score: f64,
        member_sequences: Vec<String>,
    ) -> Self {
        Self {
            pattern,
            kind,
            fisher_score,
            member_sequences,
        }
    }
}

impl PartialEq for Cluster {
    fn eq(&self, other: &Self) -> bool {
        self.pattern == other.pattern
            && self.kind == other.kind
            && self.fisher_score == other.fisher_score
            && self.member_sequences == other.member_sequences
    }
}

impl Hash for Cluster {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.pattern.hash(state);
        self.kind.hash(state);
        self.fisher_score.to_bits().hash(state);
        self.member_sequences.hash(state);
    }
}

#[derive(Clone, Debug)]
pub struct ClusterMembers<'a> {
    pub records: Vec<&'a ChainRecord>,
    pub subject_count: usize,
    pub cluster: &'a Cluster,
    pub member_sequences: Arc<HashSet<&'a str>>,
    pub cluster_subjects: Arc<HashSet<&'a str>>,
}

impl<'a> ClusterMembers<'a> {
    /// Members are exactly `records[i]` for each of `indices`.
    pub fn from_indices(
        records: &'a [ChainRecord],
        cluster: &'a Cluster,
        indices: &[usize],
    ) -> Self {
        let cluster_records: Vec<&ChainRecord> = indices.iter().map(|&i| &records[i]).collect();
        let member_set: HashSet<&str> = cluster_records.iter().map(|r| r.beta.as_str()).collect();
        Self::from_records(cluster_records, cluster, member_set)
    }

    fn from_records(
        cluster_records: Vec<&'a ChainRecord>,
        cluster: &'a Cluster,
        member_set: HashSet<&'a str>,
    ) -> Self {
        let cluster_subjects: HashSet<&'a str> = cluster_records
            .iter()
            .map(|record| record.subject.as_str())
            .collect();
        let subject_count = cluster_subjects.len();
        Self {
            records: cluster_records,
            subject_count,
            cluster,
            member_sequences: Arc::from(member_set),
            cluster_subjects: Arc::from(cluster_subjects),
        }
    }
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
}

/// Anything that wraps a shared [`Cluster`] together with its member records
/// (GLIPH clusters or tcrdist paired clusters).
pub trait HasCluster {
    fn cluster(&self) -> &Cluster;

    /// Indices of the member records, in ascending order.
    fn member_indices(&self) -> &[usize];

    /// Motif shared by the members' CDR3α sequences, if any (the `pattern_alpha` column).
    fn alpha_pattern(&self) -> Option<&str> {
        None
    }

    fn members<'a>(&'a self, records: &'a [ChainRecord]) -> ClusterMembers<'a> {
        ClusterMembers::from_indices(records, self.cluster(), self.member_indices())
    }
}

/// A cluster together with its score and the output rows generated for its members.
#[derive(Debug)]
pub struct ScoredCluster<C> {
    pub cluster: C,
    pub score: ClusterScore,
    pub rows: Vec<OutputRow>,
}
