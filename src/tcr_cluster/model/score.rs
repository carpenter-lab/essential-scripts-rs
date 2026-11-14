#[derive(Clone, Debug)]
pub struct ClusterScore {
    pub final_score: f64,
    pub hla_score: f64,
    pub vb_score: f64,
    pub expansion_score: f64,
    pub length_score: f64,
    pub cluster_size_score: f64,
    pub hla_hit: Vec<(usize, String)>,
}
