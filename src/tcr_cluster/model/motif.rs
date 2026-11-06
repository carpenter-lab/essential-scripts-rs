use std::collections::BTreeSet;

#[derive(Clone, Debug)]
pub struct MotifStats {
    pub members: BTreeSet<String>,
    pub fisher: f64,
    pub fold_change: f64,
}

#[derive(Clone, Debug)]
pub struct TaggedMember {
    pub sequence: String,
    pub residue: u8,
    pub position: usize,
}
