use itertools::Itertools;

use autocompress::{
    create, iothread::IoThread, suggest_format_from_path, CompressionLevel, Encoder,
};
use rand::{seq::*, Rng};
use rayon::prelude::*;
use seq_io::fasta::Reader;
use seq_io::policy::BufPolicy;
use seq_io::{prelude::*, PositionStore};
use std::collections::HashMap;
use std::fmt::Display;
use std::io::{BufWriter, Read};
use serde::Serialize;
use tabled::Tabled;
use std::path::PathBuf;

const MISSING_POS: u32 = u32::MAX;

#[derive(PartialEq, Debug, Copy, Clone, Serialize, Tabled)]
pub struct SpResult {
    pub shared_homologies : u64,
    pub ref_homologies : u64,
    pub est_homologies : u64,
    pub spfp : f64,
    pub spfn : f64,
}

impl SpResult {
    pub fn new(shared : u64, ref_homs : u64, est_homs : u64) -> Self {
        let spfn = ((ref_homs - shared) as f64) / (ref_homs as f64);
        let spfp = ((est_homs - shared) as f64) / (est_homs as f64);
        Self {
            shared_homologies: shared,
            ref_homologies: ref_homs,
            est_homologies: est_homs,
            spfp,
            spfn,
        }
    }
}

// See Siavash Mirarab, Tandy Warnow, FASTSP: linear time calculation of alignment accuracy, Bioinformatics
pub fn calc_fpfn(reference: &PathBuf, estimated: &PathBuf) -> SpResult {
    let thread_pool = IoThread::new(2);
    let mut est_reader = Reader::new(thread_pool.open(estimated).unwrap());
    let mut it = est_reader.records().peekable();
    let r = it.peek();
    let width = r
        .unwrap()
        .as_ref()
        .unwrap()
        .full_seq()
        .iter()
        .filter(|it| !it.is_ascii_lowercase())
        .count(); // I am pretty sure this is a horrible idea...
    drop(it);
    est_reader = Reader::new(thread_pool.open(estimated).unwrap());

    let mut s: Vec<Vec<u32>> = vec![]; // s[i, j] = column in the estimated alignment of the i-th row, j-th non-gap char
    let mut i: usize = 0;
    let mut rows = 0u32;
    let mut gaps = vec![0u32; width];
    let mut names : HashMap<String, usize> = HashMap::new();
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
                    let ind = if c.is_ascii_lowercase() {
                        gaps[x] += 1;
                        MISSING_POS
                    } else {
                        x as u32
                    };
                    s[i].push(ind);
                    // x += 1;
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
    // println!("{:?}", gaps);
    // println!("{:?}", rows);
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
        .filter(|it| !it.is_ascii_lowercase())
        .count(); // this is horrible
    ref_reader = Reader::new(thread_pool.open(reference).unwrap());
    let mut columns: Vec<Vec<(u32, u32)>> = vec![vec![]; ref_width]; // columns[i] contains the N(i, j) char
    let mut ref_gaps = vec![0u32; ref_width];
    let mut ref_rows = 0u32;
    // i = 0;
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
                    if c.is_ascii_lowercase() {
                        // lower case letter. We just ignore it
                        ref_gaps[x] += 1;
                    } else {
                        // actual letter. Put it in the column
                        columns[x].push((i as u32, j as u32));
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
    let ref_homologies = ref_gaps
        .iter()
        .map(|g| (ref_rows - *g) as u64 * (ref_rows - *g - 1) as u64 / 2)
        .sum::<u64>();
    let shared_homologies = columns.par_iter().map(|it| {
        let mut counter : HashMap<u32, u32> = HashMap::new();
        for (u, v) in it.iter() {
            let c = s[*u as usize][*v as usize];
            if c == MISSING_POS {
                continue;
            }
            let entry = counter.entry(c).or_insert(0);
            *entry += 1;
        }
        let mut cnt = 0;
        for v in counter.values() {
            cnt += *v as u64 * (*v as u64 - 1) / 2;
        }
        cnt
    }).sum::<u64>();
    // println!("{:?}", s);
    // println!("{:?}", &columns);
    // let d = columns.iter().map(|it| it.iter().map(|(u, v)| s[*u as usize][*v as usize]).collect_vec()).collect_vec();
    // println!("{:?}", d);
    return SpResult::new(shared_homologies, ref_homologies, est_homologies);
}