use ahash::AHashSet;
use anyhow::bail;
use autocompress::{iothread::IoThread, CompressionLevel};
use clap::ArgEnum;
use itertools::Itertools;
use rand::prelude::ThreadRng;
use rand::{seq::*, Rng};
use rayon::prelude::*;
use seq_io::fasta::{Reader, RefRecord};
use seq_io::{prelude::*, PositionStore};
use serde::Serialize;
use similar::TextDiff;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Display;
use tabled::builder::Builder;

// needed to import necessary traits
use std::{
    io::Write,
    path::{Path, PathBuf},
};
use tabled::{Style, Table, Tabled};

const AA_WILDCARD: u8 = b'X';
const NC_WILDCARD: u8 = b'N';

pub fn write_record<W>(writer: &mut W, head: &[u8], seq: &[u8]) -> anyhow::Result<()>
where
    W: Write,
{
    writer.write_all(b">")?;
    writer.write_all(head)?;
    writer.write_all(b"\n")?;
    seq.chunks(60)
        .try_for_each::<_, anyhow::Result<()>>(|chunk| {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
            Ok(())
        })?;
    Ok(())
}

#[derive(Debug, Eq, PartialEq, Clone, Copy, Serialize)]
pub enum Alphabet {
    AminoAcid,
    Nucleotide,
}

impl Display for Alphabet {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Alphabet::AminoAcid => write!(f, "AminoAcid"),
            Alphabet::Nucleotide => write!(f, "Nucleotide"),
        }
    }
}

pub fn p_distance(lhs: &[u8], rhs: &[u8], ambiguous: u8) -> Option<(usize, usize, f64)> {
    let mut tt = 0;
    let mut diff = 0;
    lhs.iter().zip(rhs.iter()).for_each(|(l, r)| {
        match (l, r) {
            (b'-', _) => {}
            (_, b'-') => {}
            (a, _) if *a == ambiguous => {}
            (_, b) if *b == ambiguous => {}
            (a, b) => {
                tt += 1;
                diff += if *a == *b { 0 } else { 1 };
            }
        };
    });
    if tt == 0 {
        None
    } else {
        Some((diff, tt, diff as f64 / tt as f64))
    }
}

pub fn all_pairs_p_distance(matrix: &[Vec<u8>], alph: Alphabet) -> Option<(f64, f64)> {
    // let mut p_dis = vec![];
    let wildcard = if alph == Alphabet::AminoAcid {
        AA_WILDCARD
    } else {
        NC_WILDCARD
    };
    let p_dis: Vec<(usize, usize, f64)> = (0..matrix.len())
        .combinations(2)
        .par_bridge()
        .flat_map(|v| p_distance(&matrix[v[0]], &matrix[v[1]], wildcard))
        .collect();

    if p_dis.is_empty() {
        None
    } else {
        let top = p_dis.iter().map(|i| i.0).sum::<usize>();
        let bottom = p_dis.iter().map(|i| i.1).sum::<usize>();
        let avg = top as f64 / bottom as f64;
        let max = p_dis.iter().fold(0f64, |acc, x| acc.max(x.2));
        Some((avg, max))
    }
}

#[derive(Debug, Tabled)]
pub struct AlnSimpleStats {
    pub alph: Alphabet,
    pub total_cells: usize,
    pub gap_cells: usize,
    pub width: usize,
    pub rows: usize,
    pub avg_sequence_length: f64,
    pub median_sequence_length: f64,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum, Debug)]
pub enum Approx {
    Auto,
    Yes,
    No,
}

pub struct AlignmentSampler {
    pub rng: ThreadRng,
    pub records: Vec<Vec<u8>>,
    pub max_capacity: Option<usize>,
    pub i: usize,
}

impl AlignmentSampler {
    pub fn new(max_capacitiy: Option<usize>) -> Self {
        Self {
            rng: rand::thread_rng(),
            records: vec![],
            max_capacity: max_capacitiy,
            i: 0,
        }
    }

    pub fn is_full(&self) -> bool {
        self.max_capacity.map_or(false, |b| self.i >= b)
    }

    pub fn see(&mut self, record: &RefRecord) {
        match self.max_capacity {
            None => {
                self.records.push(record.to_owned_record().seq);
            }
            Some(s) => {
                if self.i < s {
                    self.records.push(record.to_owned_record().seq);
                } else {
                    let j = self.rng.gen_range(0..(self.i + 1));
                    if j < s {
                        self.records[j] = record.to_owned_record().seq;
                    }
                }
            }
        }
        self.i += 1;
    }

    pub fn compute_pdis(
        &mut self,
        alph: Alphabet,
        subsamples: bool,
    ) -> Result<PdisResult, &'static str> {
        let mut avg_pdis_samples: Vec<f64> = vec![];
        let mut max_pdis_samples: Vec<f64> = vec![];
        if subsamples {
            self.records.shuffle(&mut self.rng);
            let samples = self
                .records
                .chunks(1000)
                .flat_map(|it| all_pairs_p_distance(it, alph));
            samples.for_each(|(l, r)| {
                avg_pdis_samples.push(l);
                max_pdis_samples.push(r);
            });
        } else {
            all_pairs_p_distance(&self.records, alph)
                .iter()
                .for_each(|(a, b)| {
                    avg_pdis_samples.push(*a);
                    max_pdis_samples.push(*b);
                });
        }

        if avg_pdis_samples.is_empty() {
            Err("No valid samples obtained")
        } else {
            // compute the median of the avg pdis
            // the max of the max pdis
            avg_pdis_samples.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let avg_pdis = avg_pdis_samples[avg_pdis_samples.len() / 2];
            let max_pdis = max_pdis_samples
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max);
            Ok(PdisResult {
                avg_pdis,
                max_pdis,
                approx: subsamples,
            })
        }
    }
}

#[derive(Debug, Serialize)]
pub struct CombinedAlnStats {
    alph: Alphabet,
    columns: usize,
    rows: usize,
    gap_ratio: f64,
    avg_seq_length: f64,
    median_seq_length: f64,
    avg_p_dis: Option<f64>,
    max_p_dis: Option<f64>,
    p_dis_approx: Option<bool>,
}

impl CombinedAlnStats {
    pub fn new(stats: &AlnSimpleStats, pdis: &Option<PdisResult>) -> Self {
        Self {
            alph: stats.alph,
            columns: stats.width,
            rows: stats.rows,
            gap_ratio: stats.gap_cells as f64 / stats.total_cells as f64,
            avg_seq_length: stats.avg_sequence_length,
            median_seq_length: stats.median_sequence_length,
            avg_p_dis: pdis.as_ref().map(|p| p.avg_pdis),
            max_p_dis: pdis.as_ref().map(|p| p.max_pdis),
            p_dis_approx: pdis.as_ref().map(|p| p.approx),
        }
    }
}

fn median_of(coll: &mut [u64]) -> f64 {
    if coll.len() % 2 != 0 {
        coll.select_nth_unstable(coll.len() / 2);
        coll[coll.len() / 2] as u32 as f64
    } else {
        let l = coll.len() / 2 - 1;
        coll.select_nth_unstable(l);
        let l_v = coll[l] as u32;
        coll.select_nth_unstable(l + 1);
        let r_v = coll[l + 1] as u32;
        (l_v + r_v) as f64 / 2.0
    }
}

pub fn aln_linear_stats<P>(
    filename: P,
    p_dis: bool,
    approx: Approx,
) -> (AlnSimpleStats, Option<PdisResult>)
where
    P: AsRef<Path>,
{
    let mut tt_cells = 0u64;
    let mut gap_cells = 0u64;
    let mut width = 0;
    let mut rows = 0;
    let mut seq_lengths: Vec<u64> = vec![];
    let thread_pool = IoThread::new(3);
    let mut reader = Reader::new(thread_pool.open(filename).unwrap());
    let mut guesser = AlphabetGuesser::new();
    let mut sampler = AlignmentSampler::new(match (p_dis, approx) {
        (false, _) => None,
        (true, Approx::No) => None,
        _ => Some(9000),
    });
    while let Some(result) = reader.next() {
        let mut record_width = 0;
        let mut local_gaps = 0u64;
        let rec = result.unwrap();
        for l in rec.seq_lines() {
            tt_cells += l.len() as u64;
            record_width += l.len();
            for c in l {
                guesser.see(c);
                if *c == b'-' {
                    local_gaps += 1;
                }
            }
        }
        if p_dis {
            sampler.see(&rec);
        }
        gap_cells += local_gaps;
        if width == 0 {
            width = record_width;
        } else {
            assert_eq!(width, record_width);
        }
        rows += 1;
        seq_lengths.push(record_width as u64 - local_gaps as u64);
    }
    let stats = AlnSimpleStats {
        alph: guesser.alph(),
        total_cells: tt_cells as usize,
        gap_cells: gap_cells as usize,
        width,
        rows,
        avg_sequence_length: seq_lengths.iter().sum::<u64>() as f64 / seq_lengths.len() as f64,
        median_sequence_length: median_of(&mut seq_lengths),
    };
    let pdis_result = match (p_dis, approx) {
        (false, _) => None,
        (true, Approx::No) => Some(sampler.compute_pdis(guesser.alph(), false).unwrap()),
        (true, Approx::Auto) if !sampler.is_full() => {
            Some(sampler.compute_pdis(guesser.alph(), false).unwrap())
        }
        _ => Some(sampler.compute_pdis(guesser.alph(), true).unwrap()),
    };
    (stats, pdis_result)
}

pub struct AlphabetGuesser {
    pub seen_aa_char: bool,
}

impl AlphabetGuesser {
    pub fn new() -> Self {
        AlphabetGuesser {
            seen_aa_char: false,
        }
    }

    pub fn see(&mut self, c: &u8) {
        match c.to_ascii_uppercase() {
            b'A' | b'C' | b'T' | b'G' | b'-' | b'N' | b'U' => {}
            _ => {
                self.seen_aa_char = true;
            }
        }
    }

    pub fn alph(&self) -> Alphabet {
        if self.seen_aa_char {
            Alphabet::AminoAcid
        } else {
            Alphabet::Nucleotide
        }
    }
}

#[derive(Debug, Tabled, Serialize)]
pub struct MaskResult {
    pub total_columns: usize,
    pub masked_columns: usize,
    pub rest_columns: usize,
    pub total_rows: usize,
}

#[derive(Debug, Tabled, Serialize)]
pub struct WhereResult {
    pub total_rows: usize,
    pub matched_rows: usize,
    #[tabled(display_with = "display_ratio")]
    pub force_included: (usize, usize),
}

fn display_ratio(v: &(usize, usize)) -> String {
    format!("{}/{}", v.0, v.1)
}

#[derive(Debug, Tabled, Clone)]
pub struct PdisResult {
    pub avg_pdis: f64,
    pub max_pdis: f64,
    pub approx: bool,
}

pub fn parse_sequence_set(files: &[PathBuf]) -> AHashSet<String> {
    let mut res = AHashSet::new();
    let thread_pool = IoThread::new(2);
    for f in files {
        let mut reader = Reader::new(thread_pool.open(f).unwrap());
        while let Some(result) = reader.next() {
            let name =
                String::from_utf8(result.unwrap().head().iter().copied().collect_vec()).unwrap();
            res.insert(name);
        }
    }
    res
}

pub fn aln_where(
    filename: &PathBuf,
    length_lb: Option<usize>,
    length_ub: Option<usize>,
    inclusion: &[PathBuf],
    outfile: &PathBuf,
) -> WhereResult {
    let force_inclusion_set = parse_sequence_set(inclusion);
    let mut force_included = 0usize;
    let thread_pool = IoThread::new(2);
    let mut reader = Reader::new(thread_pool.open(filename).unwrap());
    let mut rows = 0usize;
    let mut matched = 0usize;
    let mut writer = thread_pool
        .create(outfile, CompressionLevel::Default)
        .unwrap();
    while let Some(result) = reader.next() {
        let mut nongaps = 0usize;
        let mut buf: Vec<u8> = vec![];
        let r = result.unwrap();
        let mut force_include = false;
        if !force_inclusion_set.is_empty() {
            let name = String::from_utf8(r.head().iter().copied().collect_vec()).unwrap();
            if force_inclusion_set.contains(&name) {
                force_include = true;
            }
        }
        rows += 1;
        for l in r.seq_lines() {
            for c in l {
                buf.push(*c);
                if *c != b'-' {
                    nongaps += 1;
                }
            }
        }
        let a = length_lb.map(|x| nongaps >= x);
        let b = length_ub.map(|x| nongaps <= x);
        match (a, b, force_include) {
            (_, _, true) => {
                force_included += 1;
            }
            (Some(false), _, _) => continue,
            (_, Some(false), _) => continue,
            (_, _, _) => {}
        }
        write_record(&mut writer, r.head(), &buf).unwrap();
        matched += 1;
    }
    WhereResult {
        total_rows: rows,
        matched_rows: matched,
        force_included: (force_included, force_inclusion_set.len()),
    }
}

pub fn aln_mask(
    filename: &PathBuf,
    num_sites: usize,
    percent_gappy: f64,
    ignore_case: bool,
    outfile: &PathBuf,
) -> MaskResult {
    let thread_pool = IoThread::new(3);
    let mut reader = Reader::new(thread_pool.open(filename).unwrap());
    let mut it = reader.records().peekable();
    let r = it.peek();
    let width = r.unwrap().as_ref().unwrap().full_seq().len();
    reader = Reader::new(thread_pool.open(filename).unwrap());
    let mut height = 0usize;
    let mut gaps = vec![0u64; width];
    let mut n_char = vec![0u64; width];
    let mut x_char = vec![0u64; width];
    let mut guesser = AlphabetGuesser::new();
    while let Some(result) = reader.next() {
        let mut pos = 0usize;
        for l in result.unwrap().seq_lines() {
            for c in l {
                if !ignore_case && c.is_ascii_lowercase() {
                    gaps[pos] += 1;
                    pos += 1;
                    continue;
                }

                guesser.see(c);
                match c {
                    b'-' => {
                        gaps[pos] += 1;
                    }
                    b'N' => {
                        n_char[pos] += 1;
                    }
                    b'X' => {
                        x_char[pos] += 1;
                    }
                    _ => {}
                }
                pos += 1;
            }
        }
        height += 1;
    }

    let remove: Vec<bool> = (0..width)
        .map(|i| {
            let gap_count = match guesser.alph() {
                Alphabet::AminoAcid => gaps[i] + x_char[i],
                Alphabet::Nucleotide => gaps[i] + n_char[i],
            };
            (gap_count >= height as u64 - num_sites as u64)
                || (gap_count as f64 >= height as f64 * percent_gappy)
        })
        .collect();

    let mut r2 = Reader::new(thread_pool.open(filename).unwrap());
    let mut writer = thread_pool
        .create(outfile, CompressionLevel::Default)
        .unwrap();
    while let Some(result) = r2.next() {
        let mut pos = 0usize;
        let mut buf: Vec<u8> = vec![];
        let r = result.unwrap();
        for l in r.seq_lines() {
            for c in l {
                if !remove[pos] {
                    if !ignore_case && c.is_ascii_lowercase() {
                        buf.push(b'-');
                    } else {
                        buf.push(*c);
                    }
                }
                pos += 1;
            }
        }
        write_record(&mut writer, r.head(), &buf).unwrap();
    }
    let removed = remove.iter().filter(|x| **x).count();
    MaskResult {
        total_columns: width,
        masked_columns: removed,
        rest_columns: width - removed,
        total_rows: height,
    }
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum, Debug)]
pub enum AlnDiffMode {
    Udiff,
    Set,
    Fragile,
}

#[derive(Clone, PartialEq, Serialize)]
pub struct AlnDiffSummary {
    different: Vec<String>,
    lhs_only: Vec<String>,
    rhs_only: Vec<String>,
}

pub fn aln_diff_summary_table(summary: &AlnDiffSummary) -> Table {
    let mut b = Builder::default().set_columns(["Name", "Discrepancy"]);
    for d in &summary.different {
        b = b.add_record([d, "different"]);
    }
    for d in &summary.lhs_only {
        b = b.add_record([d, "lhs only"]);
    }
    for d in &summary.rhs_only {
        b = b.add_record([d, "rhs only"]);
    }
    b.build().with(Style::modern())
}

pub fn aln_diff(
    lhs: &PathBuf,
    rhs: &PathBuf,
    degap: bool,
    mode: AlnDiffMode,
) -> anyhow::Result<Option<AlnDiffSummary>> {
    let mut lhs_seqs: BTreeMap<String, String> = BTreeMap::new();
    let mut rhs_seqs: BTreeMap<String, String> = BTreeMap::new();
    let thread_pool = IoThread::new(2);
    let mut r1 = Reader::new(thread_pool.open(lhs)?);
    while let Some(r) = r1.next() {
        let record = r?;
        let name = String::from_utf8(record.head().to_vec())?;
        let seq = if degap {
            String::from_utf8(
                record
                    .full_seq()
                    .iter()
                    .filter(|it| **it != b'-')
                    .copied()
                    .collect_vec(),
            )?
        } else {
            String::from_utf8(record.full_seq().to_vec())?
        };
        lhs_seqs.insert(name, seq);
    }
    let mut r2 = Reader::new(thread_pool.open(rhs)?);
    while let Some(r) = r2.next() {
        let record = r?;
        let name = String::from_utf8(record.head().to_vec())?;
        let seq = if degap {
            String::from_utf8(
                record
                    .full_seq()
                    .iter()
                    .filter(|it| **it != b'-')
                    .copied()
                    .collect_vec(),
            )?
        } else {
            String::from_utf8(record.full_seq().to_vec())?
        };
        rhs_seqs.insert(name, seq);
    }
    match mode {
        AlnDiffMode::Udiff => {
            let mut lhs_buf = String::new();
            let mut rhs_buf = String::new();
            for (k, v) in lhs_seqs {
                lhs_buf.push_str(&format!(">{}\n{}\n", k, v));
            }
            for (k, v) in rhs_seqs {
                rhs_buf.push_str(&format!(">{}\n{}\n", k, v));
            }
            let diff = TextDiff::from_lines(&lhs_buf, &rhs_buf);
            print!(
                "{}",
                diff.unified_diff().context_radius(1).header("lhs", "rhs")
            );
        }
        AlnDiffMode::Set => {
            let lhs_keys = BTreeSet::from_iter(lhs_seqs.keys());
            let rhs_keys = BTreeSet::from_iter(rhs_seqs.keys());
            let shared = lhs_keys.intersection(&rhs_keys);
            let diff_lhs = lhs_keys.difference(&rhs_keys);
            let diff_rhs = rhs_keys.difference(&lhs_keys);
            let mut different: Vec<String> = vec![];
            let mut lhs_only: Vec<String> = vec![];
            let mut rhs_only: Vec<String> = vec![];
            for s in shared {
                if lhs_seqs.get(*s) != rhs_seqs.get(*s) {
                    different.push(s.to_string());
                }
            }
            for s in diff_lhs {
                lhs_only.push(s.to_string());
            }
            for s in diff_rhs {
                rhs_only.push(s.to_string());
            }
            return Ok(Some(AlnDiffSummary {
                different,
                lhs_only,
                rhs_only,
            }));
        }
        AlnDiffMode::Fragile => {
            if lhs_seqs != rhs_seqs {
                bail!("Two alignments are not equal");
            }
        }
    }
    Ok(None)
}
