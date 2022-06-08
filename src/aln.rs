use autocompress::iothread::ThreadReader;
use itertools::Itertools;
use rand::distributions::uniform::SampleBorrow;
use rand::{seq::*, Rng};
use rayon::prelude::*;
// use seq_io::core::BufReader;
use seq_io::fasta::{LineStore, OwnedRecord, Reader, RecordsIter, RefRecord};
// use lz4::{Decoder, EncoderBuilder, Encoder};

use seq_io::policy::{BufPolicy, StdPolicy};
// use std::fmt::Display;
use autocompress::{create, iothread::IoThread, open, CompressionLevel, Decoder, Encoder};
use seq_io::{prelude::*, Position, PositionStore};
use std::fmt::Display;
use std::io::{BufReader, BufWriter, Read};
// needed to import necessary traits
use std::{
    fs::File,
    io::Write,
    path::{Path, PathBuf},
};
use tabled::Tabled;
// use thread_io::read::reader;

const AA_WILDCARD: u8 = b'X';
const NC_WILDCARD: u8 = b'N';

#[derive(Debug, Eq, PartialEq, Clone, Copy)]
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

pub trait AlnReader<S>
where
    S: PositionStore,
{
    fn next(
        &mut self,
    ) -> std::option::Option<Result<seq_io::fasta::RefRecord<'_, S>, seq_io::fasta::Error>>;
    fn peek_size(&mut self) -> usize;
    // fn records(&mut self) -> RecordsIter<'_, R, P, S>;
}

impl<R, P, S> AlnReader<S> for Reader<R, P, S>
where
    R: Read,
    P: BufPolicy,
    S: PositionStore,
{
    fn next(
        &mut self,
    ) -> std::option::Option<Result<seq_io::fasta::RefRecord<'_, S>, seq_io::fasta::Error>> {
        return self.next();
    }

    fn peek_size(&mut self) -> usize {
        return self
            .records()
            .peekable()
            .peek()
            .unwrap()
            .as_ref()
            .unwrap()
            .full_seq()
            .len();
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
}

pub fn aln_linear_stats<P>(filename: P) -> AlnSimpleStats
where
    P: AsRef<Path>,
{
    let mut tt_cells = 0u64;
    let mut gap_cells = 0u64;
    let mut width = 0;
    let mut rows = 0;
    let mut seq_lengths: Vec<u64> = vec![];
    let mut reader = Reader::from_path(filename).unwrap();
    let mut guesser = AlphabetGuesser::new();
    while let Some(result) = reader.next() {
        let mut record_width = 0;
        let mut local_gaps = 0u64;
        for l in result.unwrap().seq_lines() {
            tt_cells += l.len() as u64;
            record_width += l.len();
            for c in l {
                guesser.see(c);
                if *c == b'-' {
                    local_gaps += 1;
                }
            }
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
    AlnSimpleStats {
        alph: guesser.alph(),
        total_cells: tt_cells as usize,
        gap_cells: gap_cells as usize,
        width,
        rows,
        avg_sequence_length: seq_lengths.iter().sum::<u64>() as f64 / seq_lengths.len() as f64,
    }
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
        match c {
            b'A' | b'C' | b'T' | b'G' | b'-' | b'N' => {}
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

#[derive(Debug, Tabled)]
pub struct MaskResult {
    pub total_columns: usize,
    pub masked_columns: usize,
    pub total_rows: usize,
}

#[derive(Debug, Tabled)]
pub struct WhereResult {
    pub total_rows: usize,
    pub matched_rows: usize,
}

#[derive(Debug, Tabled, Clone)]
pub struct PdisResult {
    pub avg_pdis: f64,
    pub max_pdis: f64,
}

pub fn approx_pdis(filename: &PathBuf, alph: Alphabet) -> Result<PdisResult, &'static str> {
    let mut rng = rand::thread_rng();
    let mut reader = Reader::from_path(filename).unwrap();
    let mut i = 0;
    let s = 9000;
    let mut records: Vec<Vec<u8>> = vec![];
    //reservoir sampling
    while let Some(result) = reader.next() {
        let r = result.unwrap();
        if i < s {
            records.push(r.to_owned_record().seq);
        } else {
            let j = rng.gen_range(0..(i + 1));
            if j < s {
                records[j] = r.to_owned_record().seq;
            }
        }
        i += 1;
    }

    let mut avg_pdis_samples: Vec<f64> = vec![];
    let mut max_pdis_samples: Vec<f64> = vec![];
    if i <= s {
        all_pairs_p_distance(&records, alph).map(|(a, b)| {
            avg_pdis_samples.push(a);
            max_pdis_samples.push(b);
        });
    } else {
        records.shuffle(&mut rng);
        let samples = records
            .chunks(1000)
            .flat_map(|it| all_pairs_p_distance(it, alph));
        samples.for_each(|(l, r)| {
            avg_pdis_samples.push(l);
            max_pdis_samples.push(r);
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
        Ok(PdisResult { avg_pdis, max_pdis })
    }
}

// pub fn retrieve_aln_reader<'a>(filename : &PathBuf, thread_pool : &'a IoThread) ->  {
//     let x = thread_pool.open(filename).unwrap();
//     let y = Reader::new(x);
//     return y;
// }

pub fn retrieve_writer(filename: &PathBuf) -> BufWriter<Encoder<File>> {
    let x = create(filename, CompressionLevel::Default).unwrap();
    return BufWriter::new(x);
}

pub fn aln_where(
    filename: &PathBuf,
    length_lb: Option<usize>,
    length_ub: Option<usize>,
    outfile: &PathBuf,
) -> WhereResult {
    let thread_pool = IoThread::new(2);
    let mut reader = Reader::new(thread_pool.open(filename).unwrap());
    // let mut reader = retrieve_aln_reader(filename);
    let mut rows = 0usize;
    let mut matched = 0usize;
    let mut writer = thread_pool
        .create(outfile, CompressionLevel::Default)
        .unwrap();
    while let Some(result) = reader.next() {
        let mut nongaps = 0usize;
        let mut buf: Vec<u8> = vec![];
        let r = result.unwrap();
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
        match (a, b) {
            (Some(false), _) => continue,
            (_, Some(false)) => continue,
            (_, _) => {
                writer.write_all(b">").unwrap();
                writer.write_all(r.head()).unwrap();
                writer.write_all(b"\n").unwrap();
                buf.chunks(60).for_each(|chunk| {
                    writer.write_all(chunk).unwrap();
                    writer.write_all(b"\n").unwrap();
                });
                matched += 1;
            }
        }
    }
    WhereResult {
        total_rows: rows,
        matched_rows: matched,
    }
}

pub fn aln_mask(
    filename: &PathBuf,
    num_sites: usize,
    percent_gappy: f64,
    outfile: &PathBuf,
) -> MaskResult {
    let thread_pool = IoThread::new(2);
    let mut reader = Reader::new(thread_pool.open(filename).unwrap());
    let mut it = reader.records().peekable();
    let r = it.peek();
    let width = r.unwrap().as_ref().unwrap().full_seq().len();
    let mut height = 0usize;
    let mut gaps = vec![0u64; width];
    let mut n_char = vec![0u64; width];
    let mut x_char = vec![0u64; width];
    let mut guesser = AlphabetGuesser::new();
    while let Some(result) = reader.next() {
        let mut pos = 0usize;
        for l in result.unwrap().seq_lines() {
            for c in l {
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
                    buf.push(*c);
                }
                pos += 1;
            }
        }
        writer.write_all(b">").unwrap();
        writer.write_all(r.head()).unwrap();
        writer.write_all(b"\n").unwrap();
        buf.chunks(60).for_each(|chunk| {
            writer.write_all(chunk).unwrap();
            writer.write_all(b"\n").unwrap();
        });
    }
    MaskResult {
        total_columns: width,
        masked_columns: remove.iter().filter(|x| **x).count(),
        total_rows: height,
    }
}
