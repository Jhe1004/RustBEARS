//! Branch histories conditional on observed endpoints.
mod bridge;
mod legacy;

pub use bridge::{sample_batch, BridgeJob, Generator, Xoshiro256PlusPlus};
pub use legacy::{reject, LegacyResult};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Event {
    pub from: usize,
    pub to: usize,
    pub time: f64,
}

mod tree;
pub use tree::{sample_tree, LayerInput, PreparedTreeSampler, SegmentInput, TreeInput};
