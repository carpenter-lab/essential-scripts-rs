#[derive(Debug, Clone)]
pub struct GlobalFilterConfig {
    pub fisher_threshold: f64,
    pub min_members: usize,
    pub max_reference_fraction: f64,
}

impl GlobalFilterConfig {
    pub fn new(fisher_threshold: f64, min_members: usize, max_reference_fraction: f64) -> Self {
        Self {
            fisher_threshold,
            min_members,
            max_reference_fraction,
        }
    }
}

impl Default for GlobalFilterConfig {
    fn default() -> Self {
        Self::new(0.01, 2, 0.01)
    }
}

/// Organism whose germline V-gene CDRs are used by the tcrdist "all CDRs" modes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum Organism {
    #[default]
    Human,
    Mouse,
}

impl std::fmt::Display for Organism {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Human => "human",
            Self::Mouse => "mouse",
        })
    }
}

/// Which chains and CDRs contribute to the tcrdist distance between two TCRs.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum DistanceMode {
    /// CDR3β only.
    Beta,
    /// CDR3α + CDR3β.
    AlphaBeta,
    /// CDR3β + germline CDR1β, CDR2β, and CDR2.5β (from TRBV).
    BetaFull,
    /// CDR3α + CDR3β + germline CDR1, CDR2, and CDR2.5 of both chains (from TRAV and TRBV).
    AlphaBetaFull,
}

impl DistanceMode {
    pub fn uses_alpha(self) -> bool {
        matches!(self, Self::AlphaBeta | Self::AlphaBetaFull)
    }

    /// The largest number of CDR3 edits clustered in this mode. Chosen on a 160k-record
    /// repertoire: at these budgets single linkage already joins a large share of the
    /// records into one component in most modes (shared-neighbor linkage is needed to split
    /// it), and an edit or two more puts nearly every record into it.
    pub fn max_cdr3_edits(self) -> u8 {
        match self {
            Self::Beta => 2,
            Self::BetaFull => 3,
            Self::AlphaBeta => 5,
            Self::AlphaBetaFull => 8,
        }
    }

    /// Whether the germline CDR1/CDR2/CDR2.5 inferred from the V genes are included.
    pub fn uses_germline_cdrs(self) -> bool {
        matches!(self, Self::BetaFull | Self::AlphaBetaFull)
    }

    /// The mode used when none is given: alpha + beta CDR3 if any record has a CDR3α,
    /// otherwise beta only. Germline CDRs are never enabled implicitly.
    pub fn infer(records: &[super::ChainRecord]) -> Self {
        if records.iter().any(super::ChainRecord::has_alpha_cdr3) {
            Self::AlphaBeta
        } else {
            Self::Beta
        }
    }
}

impl std::fmt::Display for DistanceMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Beta => "CDR3β",
            Self::AlphaBeta => "CDR3α + CDR3β",
            Self::BetaFull => "CDR3β with all CDRs",
            Self::AlphaBetaFull => "CDR3α + CDR3β with all CDRs",
        })
    }
}

/// How the tcrdist distances are split into clusters at each threshold.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum ClusterMethod {
    /// Records within the threshold are linked only when their neighborhoods overlap
    /// enough (Jaccard index of their closed neighborhoods at least
    /// `--min-neighbor-overlap`); clusters are the connected components of the kept links.
    /// A single link between two dense groups does not join them, so chaining is limited
    /// by the same rule everywhere instead of a cluster size cap.
    #[default]
    SharedNeighbor,
    /// Connected components (single linkage) found with a union-find straight from the
    /// distance matrix. One pass covers every threshold and no graph is built.
    SingleLinkage,
    /// Connected components of an igraph graph built for each threshold. Same clusters as
    /// `single-linkage` but slower and memory-hungry at large thresholds
    ConnectedComponents,
    /// Louvain clustering of an igraph graph clustered at multiple resolutions
    #[value(skip)]
    Louvain,
}

impl ClusterMethod {
    /// Whether the method needs an igraph graph per threshold.
    pub fn needs_graph(self) -> bool {
        matches!(self, Self::ConnectedComponents | Self::Louvain)
    }
}

impl std::fmt::Display for ClusterMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::SharedNeighbor => "shared-neighbor",
            Self::SingleLinkage => "single linkage",
            Self::ConnectedComponents => "graph connected components",
            Self::Louvain => "Louvain",
        })
    }
}

/// Parameters of the tcrdist joint distance:
/// `cdr3_weight · Σ CDR3 + germline_weight · Σ (CDR1 + CDR2 + CDR2.5)` over the chains in `mode`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DistanceParams {
    pub mode: DistanceMode,
    /// tcrdist3 uses 3 for CDR3 and 1 for germline CDRs, so that the CDR3 dominates the distance.
    pub cdr3_weight: u8,
    pub germline_weight: u8,
    /// Penalty for a gap (CDR3 alignment gap, or residue vs IMGT gap in germline CDRs).
    pub gap_penalty: u8,
    /// Mismatch distance cap for the substitution matrix (tcrdist uses 4).
    pub max_mismatch: u8,
}

impl DistanceParams {
    pub fn new(mode: DistanceMode) -> Self {
        Self {
            mode,
            cdr3_weight: 3,
            germline_weight: 1,
            gap_penalty: 4,
            max_mismatch: 4,
        }
    }
    /// The tcrdist cutoff for `edits` full CDR3 edits (a mismatch costing `max_mismatch`, or
    /// one gap residue), summed over the chains in the mode and weighted like the CDR3s.
    /// Germline CDR differences count against the same budget.
    ///
    /// # Panics
    /// If the cutoff does not fit below 255 (the saturated "too far" distance).
    pub fn edit_cutoff(&self, edits: u8) -> u8 {
        edits
            .checked_mul(self.max_mismatch)
            .and_then(|x| x.checked_mul(self.cdr3_weight))
            .filter(|&x| x < u8::MAX)
            .unwrap_or_else(|| {
                panic!(
                    "{edits} CDR3 edits do not fit the u8 distance scale; \
                     use fewer edits or a smaller cdr3_weight"
                )
            })
    }
}
