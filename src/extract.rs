use autocompress::iothread::IoThread;
use clap::ArgEnum;
use itertools::Itertools;
use seq_io::{fasta::Reader, BaseRecord};
use serde::Serialize;
use std::{fmt::Display, io::Write, path::Path, str::FromStr};

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum, Debug, Serialize)]
pub enum InfoType {
    Names,
}

impl FromStr for InfoType {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "names" => Ok(InfoType::Names),
            _ => Err(format!("Unknown info type: {}", s)),
        }
    }
}

impl Display for InfoType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            InfoType::Names => write!(f, "names"),
        }
    }
}

#[derive(Debug, Serialize)]
pub struct AlnInfo {
    names: Option<Vec<String>>,
}

pub fn aln_extract<P>(filename: P, _stats: &[InfoType]) -> AlnInfo
where
    P: AsRef<Path>,
{
    let thread_pool = IoThread::new(2);
    let mut reader = Reader::new(thread_pool.open(filename).unwrap());
    let mut names: Option<Vec<String>> = Some(vec![]);
    while let Some(result) = reader.next() {
        let lined_seqs = result.unwrap();
        let name = String::from_utf8(lined_seqs.head().iter().copied().collect_vec()).unwrap();
        if let Some(names) = &mut names {
            names.push(name);
        }
    }
    let info = AlnInfo { names };
    return info;
}
