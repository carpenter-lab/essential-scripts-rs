pub mod aggregate;
pub mod copy_cellranger_outs;
pub mod enrich;
pub mod geo_submission;
#[cfg(feature = "base_cmd")]
pub mod io;
pub mod plate_reader;
pub mod split;
pub mod tcr_align;

pub mod dryad;
pub mod progress;
#[cfg(test)]
pub mod test_helpers;

pub fn underlying_clap_error_kind(error: &anyhow::Error) -> Option<clap::error::ErrorKind> {
    for cause in error.chain() {
        if let Some(e) = cause.downcast_ref::<clap::error::Error>() {
            return Some(e.kind());
        }
    }
    None
}
