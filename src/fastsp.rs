use itertools::Itertools;

use ahash::AHashMap;
use autocompress::{
    create, iothread::IoThread, suggest_format_from_path, CompressionLevel, Encoder,
};
use rand::{seq::*, Rng};
use rayon::prelude::*;
use seq_io::fasta::Reader;
use seq_io::policy::BufPolicy;
use seq_io::{prelude::*, PositionStore};
use serde::Serialize;
use std::collections::{BTreeMap, HashMap};
use std::fmt::Display;
use std::io::{BufWriter, Read};
use std::path::PathBuf;
use tabled::Tabled;

const MISSING_POS: u32 = u32::MAX;

#[derive(PartialEq, Debug, Copy, Clone, Serialize, Tabled)]
pub struct SpResult {
    pub shared_homologies: u64,
    pub ref_homologies: u64,
    pub est_homologies: u64,
    pub ref_columns: u64,
    pub est_columns: u64,
    pub spfp: f64,
    pub spfn: f64,
    pub expansion: f64,
}

impl SpResult {
    pub fn new(shared: u64, ref_homs: u64, est_homs: u64, ref_columns: u64, est_columns: u64) -> Self {
        let spfn = ((ref_homs - shared) as f64) / (ref_homs as f64);
        let spfp = ((est_homs - shared) as f64) / (est_homs as f64);
        let expansion = est_columns as f64 / ref_columns as f64;
        Self {
            shared_homologies: shared,
            ref_homologies: ref_homs,
            est_homologies: est_homs,
            ref_columns,
            est_columns,
            spfp,
            spfn,
            expansion,
        }
    }
}

// See Siavash Mirarab, Tandy Warnow, FASTSP: linear time calculation of alignment accuracy, Bioinformatics
pub fn calc_fpfn(reference: &PathBuf, estimated: &PathBuf, ignore_case : bool) -> SpResult {
    let thread_pool = IoThread::new(2);
    let mut est_reader = Reader::new(thread_pool.open(estimated).unwrap());
    let mut it = est_reader.records().peekable();
    let r = it.peek();
    let width = r.unwrap().as_ref().unwrap().full_seq().iter().count(); // I am pretty sure this is a horrible idea...
    drop(it);
    est_reader = Reader::new(thread_pool.open(estimated).unwrap());

    let mut s: Vec<Vec<u32>> = vec![]; // s[i, j] = column in the estimated alignment of the i-th row, j-th non-gap char
    let mut i: usize = 0;
    let mut rows = 0u32;
    let mut gaps = vec![0u32; width];
    let mut est_has_upper = vec![false; width];
    let mut est_num_lower = vec![0u32; width];
    let mut names: AHashMap<String, usize> = AHashMap::new();
    while let Some(result) = est_reader.next() {
        let lined_seqs = result.unwrap();
        // let mut j: usize = 0; // j is the index of the current non-gap char in the current row
        let mut x: usize = 0; // x is the current site
        let name = String::from_utf8(lined_seqs.head().iter().copied().collect_vec()).unwrap();
        names.insert(name, i);
        s.push(vec![]);
        for l in lined_seqs.seq_lines() {
            for c in l {
                if *c != b'-' {
                    let ind = if !ignore_case && c.is_ascii_lowercase() {
                        gaps[x] += 1;
                        est_num_lower[x] += 1;
                        MISSING_POS
                    } else {
                        est_has_upper[x] = true;
                        x as u32
                    };
                    s[i].push(ind);
                } else {
                    // is gap
                    gaps[x] += 1;
                }
                x += 1;
            }
        }
        i += 1;
        rows += 1;
    }

    let est_true_columns =
        est_num_lower.iter().sum::<u32>() + est_has_upper.iter().map(|it| *it as u32).sum::<u32>();

    let est_homologies = gaps
        .iter()
        .map(|g| (rows - *g) as u64 * (rows - *g - 1) as u64 / 2)
        .sum::<u64>();
    drop(gaps);
    drop(est_reader); // we have no use for the estimated anymore
    let mut ref_reader = Reader::new(thread_pool.open(reference).unwrap());
    let ref_width = ref_reader
        .records()
        .peekable()
        .peek()
        .unwrap()
        .as_ref()
        .unwrap()
        .full_seq()
        .iter()
        .count(); // this is horrible
    ref_reader = Reader::new(thread_pool.open(reference).unwrap());
    let mut columns = vec![AHashMap::new(); ref_width]; // columns[i] contains the N(i, j) char
    let mut ref_gaps = vec![0u32; ref_width];
    let mut ref_rows = 0u32;
    let mut ref_has_upper = vec![false; width];
    let mut ref_num_lower = vec![0u32; width];
    while let Some(result) = ref_reader.next() {
        let lined_seqs = result.unwrap();
        let name = String::from_utf8(lined_seqs.head().iter().copied().collect_vec()).unwrap();
        i = names[&name];
        let mut j: usize = 0; // j is the index of the current non-gap char in the current row
        let mut x: usize = 0; // x is the current site
        s.push(vec![]);
        for l in lined_seqs.seq_lines() {
            for c in l {
                if *c != b'-' {
                    if !ignore_case && c.is_ascii_lowercase() {
                        // lower case letter. We just ignore it
                        ref_num_lower[x] += 1;
                        ref_gaps[x] += 1;
                    } else {
                        // actual letter. Put it in the column
                        ref_has_upper[x] = true;
                        let entry = columns[x].entry(s[i][j]).or_insert(0);
                        *entry += 1;
                    };
                    j += 1;
                } else {
                    // is gap
                    ref_gaps[x] += 1;
                }
                x += 1;
            }
        }
        // i += 1;
        ref_rows += 1;
    }
    let ref_true_columns = ref_num_lower.iter().sum::<u32>()
        + ref_has_upper.iter().map(|it| *it as u32).sum::<u32>();
    let ref_homologies = ref_gaps
        .iter()
        .map(|g| (ref_rows - *g) as u64 * (ref_rows - *g - 1) as u64 / 2)
        .sum::<u64>();
    let shared_homologies = columns
        .par_iter()
        .map(|counter| {
            let mut cnt = 0;
            for v in counter.values() {
                cnt += *v as u64 * (*v as u64 - 1) / 2;
            }
            cnt
        })
        .sum::<u64>();
    return SpResult::new(shared_homologies, ref_homologies, est_homologies, ref_true_columns.into(), est_true_columns.into());
}
