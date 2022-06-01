use clap::ArgEnum;
use rand::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, HashSet};
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::{borrow::Borrow, cmp::max, cmp::min, collections::HashMap};
use tabled::{Table, Tabled};

#[derive(Debug)]
pub struct TaxonSet {
    pub to_id: HashMap<String, usize>,
    pub names: Vec<String>,
    last: usize,
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

impl TaxonSet {
    pub fn request(&mut self, taxon_name: String) -> usize {
        self.to_id
            .entry(taxon_name.clone())
            .or_insert_with(|| {
                self.names.push(taxon_name);
                self.last += 1;
                self.last - 1
            })
            .clone()
    }

    pub fn retreive(&self, taxon_name: String) -> usize {
        self.to_id.get(&taxon_name).unwrap().clone()
    }

    pub fn new() -> Self {
        TaxonSet {
            to_id: HashMap::new(),
            names: Vec::new(),
            last: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.last
    }
}

#[derive(Debug)]
pub struct Tree {
    pub taxa: Vec<i32>,
    pub parents: Vec<i32>,
    // pub support: Vec<f64>, // branch support
    // pub lengths: Vec<f64>, // branch lengths
    pub firstchild: Vec<i32>,
    pub nextsib: Vec<i32>,
    pub childcount: Vec<u32>,
    pub fake_root: bool,
}

impl Tree {
    pub fn children(&self, node: usize) -> ChildrenIterator {
        ChildrenIterator::new(self, node)
    }

    pub fn postorder(&self) -> PostorderIterator {
        PostorderIterator::new(self)
    }

    pub fn is_leaf(&self, node: usize) -> bool {
        if self.childcount[node] == 0 {
            return true;
        }
        false
    }

    pub fn is_root(&self, node: usize) -> bool {
        return node == 0;
    }

    pub fn xor_clades(&self, transl: &Vec<u64>, universe: u64) -> HashSet<u64> {
        let mut bips = HashSet::new();
        let mut bit_reprs = vec![0u64; self.taxa.len()];
        let mut calculated_fake_root = false;
        for node in self.postorder() {
            if self.is_leaf(node) {
                bit_reprs[node] = transl[self.taxa[node as usize] as usize];
            } else {
                if self.is_root(node) {
                    continue;
                }
                if self.is_root(self.parents[node as usize] as usize) & self.fake_root {
                    if calculated_fake_root
                        || self
                            .children(self.parents[node as usize] as usize)
                            .any(|it| self.is_leaf(it))
                    {
                        continue;
                    } else {
                        calculated_fake_root = true;
                    }
                }
                let clade = self
                    .children(node)
                    .map(|c| bit_reprs[c])
                    .reduce(|acc, i| acc ^ i)
                    .unwrap();
                bit_reprs[node] = clade;
                bips.insert(min(clade, universe ^ clade));
            }
        }
        return bips;
    }
}

pub fn compare_tree_pair(taxon_set: &TaxonSet, ref_tree: &Tree, est_tree: &Tree) -> RFOutput {
    let mut rng = rand::thread_rng();
    let n = taxon_set.len();
    let mut transl = vec![0u64; n];
    for i in 0..n {
        transl[i] = rng.gen();
    }
    let universe = transl.iter().fold(0u64, |acc, i| acc ^ i);
    let ref_bips = ref_tree.xor_clades(&transl, universe);
    let est_bips = est_tree.xor_clades(&transl, universe);
    let shared_bips = ref_bips.intersection(&est_bips).count();
    let fn_bips = ref_bips.len() - shared_bips;
    let fp_bips = est_bips.len() - shared_bips;
    return RFOutput {
        ntaxa: n,
        fp_edges: fp_bips,
        fn_edges: fn_bips,
        ref_edges: ref_bips.len(),
        est_edges: est_bips.len(),
    };
}

pub fn compare_two_trees(collection: &TreeCollection) -> RFOutput {
    // we are assuming that the first one is reference, the second one is estimated
    let ref_tree = &collection.trees[0];
    let est_tree = &collection.trees[1];
    return compare_tree_pair(&collection.taxon_set, ref_tree, est_tree);
}

pub fn compare_amplified(taxon_set: &TaxonSet, ref_tree: &Tree, est_tree: &Tree) -> RFOutput {
    let n = taxon_set.len();
    let max_tries = max(1, (n as f64).log10().floor() as u32 * 2);
    let mut prev_results = HashSet::new();
    for _ in 0..max_tries {
        let result = compare_tree_pair(taxon_set, ref_tree, est_tree);
        if prev_results.contains(&result) {
            return result;
        } else {
            prev_results.insert(result);
        }
    }
    return prev_results
        .iter()
        .max_by_key(|it| it.fn_edges + it.fp_edges)
        .unwrap()
        .clone();
}

pub fn compare_two_trees_amplified(collection: &TreeCollection) -> RFOutput {
    return compare_amplified(
        &collection.taxon_set,
        &collection.trees[0],
        &collection.trees[1],
    );
}

pub fn compare_one2many(collection: &TreeCollection) -> Vec<RFOutput> {
    return (1..(collection.ngenes()))
        .map(|i| {
            compare_amplified(
                &collection.taxon_set,
                &collection.trees[0],
                &collection.trees[i],
            )
        })
        .collect();
}

pub fn compare_many2many(collection: &TreeCollection) -> Vec<RFOutput> {
    let k = collection.ngenes() / 2;
    return (0..k)
        .map(|i| {
            compare_amplified(
                &collection.taxon_set,
                &collection.trees[i],
                &collection.trees[k + i],
            )
        })
        .collect();
}

pub fn compare_allpairs(collection: &TreeCollection) -> Vec<RFOutput> {
    let mut res = Vec::new();
    for i in 0..(collection.ngenes() - 1) {
        for j in (i + 1)..collection.ngenes() {
            res.push(compare_amplified(
                &collection.taxon_set,
                &collection.trees[i],
                &collection.trees[j],
            ));
        }
    }
    return res;
}

#[derive(Debug)]
pub struct TreeCollection {
    pub taxon_set: TaxonSet,
    pub trees: Vec<Tree>,
}

impl TreeCollection {
    pub fn new() -> Self {
        TreeCollection {
            taxon_set: TaxonSet::new(),
            trees: Vec::new(),
        }
    }

    pub fn from_multiple<P>(filenames: &[P]) -> Result<Self, String>
    where
        P: AsRef<Path> + std::fmt::Display,
    {
        if filenames.is_empty() {
            return Err("No files provided".to_string());
        } else {
            let mut tree_col = TreeCollection::from_newick(&filenames[0])?;
            for i in 1..filenames.len() {
                let added = tree_col.add_trees(&filenames[i])?;
                if added == 0 {
                    return Err(format!("{} contains no trees.", filenames[i]));
                } else if added > 1 {
                    return Err(format!(
                        "{} contains {} trees, but only one is expected.",
                        filenames[i], added
                    ));
                }
            }
            Ok(tree_col)
        }
    }

    pub fn from_newick<P>(filename: P) -> Result<Self, &'static str>
    where
        P: AsRef<Path>,
    {
        let mut trees: Vec<Tree> = vec![];
        let mut taxon_set = TaxonSet::new();
        if let Ok(lines) = read_lines(filename) {
            for line in lines {
                if let Ok(newick) = line {
                    let parsed = parse_newick(&mut taxon_set, newick.as_str());
                    trees.push(parsed);
                } else {
                    return Err("Error reading file");
                }
            }
            return Ok(TreeCollection { taxon_set, trees });
        } else {
            return Err("Could not read file");
        }
    }

    pub fn ngenes(&self) -> usize {
        self.trees.len()
    }

    pub fn ntaxa(&self) -> usize {
        self.taxon_set.len()
    }

    pub fn add_trees<P>(&mut self, filename: P) -> Result<usize, &'static str>
    where
        P: AsRef<Path>,
    {
        if let Ok(lines) = read_lines(filename) {
            let mut lines_read = 0;
            for line in lines {
                if let Ok(newick) = line {
                    let parsed = parse_newick(&mut self.taxon_set, newick.as_str());
                    self.trees.push(parsed);
                } else {
                    return Err("Error reading file");
                }
                lines_read += 1;
            }
            return Ok(lines_read);
        } else {
            return Err("Could not read file");
        }
    }
}

pub struct PostorderIterator {
    // s1 : Vec<usize>,
    s2: Vec<usize>,
}

pub struct ChildrenIterator<'a> {
    tree: &'a Tree,
    current: i32,
}

pub struct PostorderEdgeIterator<'a> {
    tree: &'a Tree,
    traversed_fakeroot: bool,
    postorder_iterator: PostorderIterator,
}

impl<'a> ChildrenIterator<'a> {
    pub fn new(tree: &'a Tree, node: usize) -> Self {
        ChildrenIterator {
            tree: tree,
            current: tree.firstchild[node],
        }
    }
}

impl<'a> PostorderEdgeIterator<'a> {
    pub fn new(tree: &'a Tree) -> Self {
        PostorderEdgeIterator {
            tree,
            traversed_fakeroot: false,
            postorder_iterator: PostorderIterator::new(tree),
        }
    }
}

impl<'a> Iterator for ChildrenIterator<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current == -1 {
            None
        } else {
            let res = self.current as usize;
            self.current = self.tree.nextsib[self.current as usize];
            Some(res)
        }
    }
}

// impl<'a> Iterator for PostorderEdgeIterator<'a> {
//     type Item = usize;

//     fn next(&mut self) -> Option<Self::Item> {
//         let candidate = self.s2.pop();
//         while let Some(node) = candidate {

//         }
//     }
// }

impl PostorderIterator {
    pub fn new(tree: &Tree) -> Self {
        let mut s1 = Vec::new();
        let mut s2 = Vec::new();
        s1.push(0usize);
        while let Some(n) = s1.pop() {
            s2.push(n);
            tree.children(n).for_each(|c| s1.push(c));
        }
        PostorderIterator {
            // s1,
            s2,
        }
    }
}

impl Iterator for PostorderIterator {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        self.s2.pop()
    }
}

#[derive(Hash, Eq, PartialEq, Debug, Copy, Clone, Serialize, Tabled)]
pub struct RFOutput {
    pub ntaxa: usize,
    pub fp_edges: usize,
    pub fn_edges: usize,
    pub ref_edges: usize,
    pub est_edges: usize,
}

#[derive(Debug, PartialEq, Serialize)]
pub struct RFPrettyOutput {
    pub raw: RFOutput,
    pub n_rf: f64,
    pub rf_rate: f64,
    pub fp_rate: f64,
    pub fn_rate: f64,
}

impl RFPrettyOutput {
    pub fn new(raw: RFOutput) -> RFPrettyOutput {
        let differing_edges = raw.fn_edges + raw.fp_edges;
        RFPrettyOutput {
            raw,
            n_rf: (differing_edges as f64) / (2 * raw.ntaxa - 6) as f64,
            rf_rate: (differing_edges as f64) / (raw.est_edges + raw.ref_edges) as f64,
            fp_rate: raw.fp_edges as f64 / raw.est_edges as f64,
            fn_rate: raw.fn_edges as f64 / raw.ref_edges as f64,
        }
    }
}

pub fn parse_newick(taxon_set: &mut TaxonSet, newick: &str) -> Tree {
    let mut taxa: Vec<i32> = vec![-42];
    let mut parents: Vec<i32> = vec![0];
    // let mut support: Vec<f64> = vec![-1.0];
    // let mut lengths: Vec<f64> = vec![-1.0];
    let mut childcount: Vec<u32> = vec![0];
    let mut firstchild: Vec<i32> = vec![-1];
    let mut nextsib: Vec<i32> = vec![-1];
    // we just reuse TreeSwift's logic
    let mut n: usize = 0; // the current node
    let mut chars = newick.chars().fuse().peekable();
    while let Some(c) = chars.next() {
        if c == ';' {
            break;
        } else if c == '(' {
            taxa.push(-1);
            childcount[n as usize] += 1;
            parents.push(n as i32);
            // support.push(0.0);
            // lengths.push(0.0);
            childcount.push(0);
            firstchild.push(-1);
            nextsib.push(-1);
            firstchild[n] = (taxa.len() - 1) as i32;
            n = taxa.len() - 1;
        } else if c == ')' {
            n = parents[n] as usize;
        } else if c == ',' {
            nextsib[n] = (taxa.len()) as i32;
            n = parents[n] as usize;
            taxa.push(-1);
            childcount[n as usize] += 1;
            parents.push(n as i32);
            // support.push(0.0);
            // lengths.push(0.0);
            childcount.push(0);
            firstchild.push(-1);
            nextsib.push(-1);
            n = taxa.len() - 1;
        } else if c == ':' {
            let mut ls = "".to_string();
            loop {
                match chars.peek() {
                    Some(',') | Some(')') | Some(';') | Some('[') | None => {
                        break;
                    }
                    Some(_) => {
                        ls.push(chars.next().unwrap());
                    }
                }
            }
            // if ls.len() > 0 {
            //     lengths[n as usize] = ls.parse::<f64>().unwrap();
            // }
        } else {
            let mut ts = c.to_string();
            loop {
                match chars.peek() {
                    Some(':') | Some(',') | Some(')') | Some(';') | Some('[') | None => {
                        break;
                    }
                    Some(_) => {
                        ts.push(chars.next().unwrap());
                    }
                }
            }
            if childcount[n] == 0 {
                taxa[n] = taxon_set.request(ts) as i32;
            }
        }
    }

    let mut fake_root = false;
    if childcount[0] == 2 {
        // let c = firstchild[0] as usize;
        // let c2 = nextsib[c] as usize;
        fake_root = true;
    }
    let tree = Tree {
        taxa,
        parents,
        // support,
        // lengths,
        firstchild,
        nextsib,
        childcount,
        fake_root,
    };
    return tree;
}
