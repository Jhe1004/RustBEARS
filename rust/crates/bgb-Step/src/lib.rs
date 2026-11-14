//! Single-stratum tree likelihood. All R APIs stay outside this crate.
mod tree;
pub use tree::{CompiledModelContext, CompiledStratifiedContext};
pub const IMPLEMENTATION_STAGE: &str = "M2_single_stratum";
