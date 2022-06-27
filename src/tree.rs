use crate::ogtree::*;
// use crate::ogtree::TreeCollection;
use serde::Serialize;
use tabled::Tabled;

#[derive(Debug, Tabled, Serialize)]
pub struct TreeStats {
    pub ntaxa: usize,
    pub num_inodes: usize,
    pub num_iedges: usize,
    pub resolution: f64,
    pub rooted: bool,
}

impl TreeStats {
    pub fn new(tree: &Tree) -> Self {
        let inodes = tree.num_nodes() - tree.ntaxa;
        let iedges = inodes - 1;
        return TreeStats {
            ntaxa: tree.ntaxa,
            num_inodes: inodes,
            num_iedges: iedges,
            resolution: iedges as f64 / (tree.ntaxa - 3) as f64,
            rooted: tree.fake_root,
        };
    }
}

pub fn tree_stats(collection: &TreeCollection) -> Vec<TreeStats> {
    collection.trees.iter().map(|t| TreeStats::new(t)).collect()
}
