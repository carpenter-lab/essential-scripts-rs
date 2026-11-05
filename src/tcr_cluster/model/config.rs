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

/// Parameters of the tcrdist joint distance:
/// `cdr3_weight · Σ CDR3 + germline_weight · Σ (CDR1 + CDR2 + CDR2.5)` over the chains in `mode`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DistanceParams {
    pub mode: DistanceMode,
    /// tcrdist3 uses 3; 1 keeps the distances on the scale of the clustering thresholds.
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
            cdr3_weight: 1,
            germline_weight: 1,
            gap_penalty: 4,
            max_mismatch: 4,
        }
    }
}
