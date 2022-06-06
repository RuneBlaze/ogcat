mod aln;
mod ogcat;

use clap::{ArgEnum, Parser, Subcommand};
use ogcat::TreeCollection;
use serde::Serialize;
use serde_json::json;
use std::fmt::Display;
use std::path::PathBuf;
use tabled::{builder::Builder, Style};
use tabled::{Table, Tabled};

use crate::ogcat::RFPrettyOutput;

#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    #[clap(subcommand)]
    cmd: SubCommand,
    #[clap(long,arg_enum, default_value_t = Format::Human)]
    format: Format,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum, Debug)]
enum Format {
    Human,
    Json,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum, Debug)]
enum Approx {
    Auto,
    Yes,
    No,
}

#[derive(Subcommand, Debug)]
enum SubCommand {
    Rf {
        #[clap(short, long = "ref")]
        reference: PathBuf,
        #[clap(short, long = "est")]
        estimated: PathBuf,
        #[clap(short, long)]
        fast: bool,
    },

    RfMulti {
        #[clap(short, long = "ref")]
        reference: PathBuf,
        #[clap(short, long = "est")]
        estimated: PathBuf,
    },

    RfAllpairs {
        #[clap(short, long = "trees", multiple_values = true)]
        trees: Vec<String>,
    },

    AlnStats {
        #[clap()]
        input: PathBuf,
        /// compute p-distance (avg, max)
        #[clap(short, long)]
        pdis: bool,
    },

    Mask {
        #[clap()]
        input: PathBuf,
        #[clap(short, long)]
        output: PathBuf,
        #[clap(short, long, default_value_t = 1f64)]
        percent: f64,
    },
}

fn two_pass_mean_stddev(x: &[f64]) -> (f64, f64) {
    let s1 = x.iter().sum::<f64>();
    let mean = s1 / x.len() as f64;
    let s2 = x.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>();
    let sample_var = (s2 as f64 / (x.len() as f64 - 1f64)) as f64;
    let sample_stddev = sample_var.sqrt();
    (mean, sample_stddev)
}

#[derive(Debug, Serialize)]
pub enum RFMultiType {
    OneToMany,
    ManyPairs,
}

impl Display for RFMultiType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RFMultiType::OneToMany => write!(f, "OneToMany"),
            RFMultiType::ManyPairs => write!(f, "ManyPairs"),
        }
    }
}

#[derive(Debug, Serialize)]
pub struct RFPairResult {
    pub pair: (usize, usize),
    pub result: ogcat::RFPrettyOutput,
}

#[derive(Debug, Serialize)]
pub struct RFPairResults {
    pub files: Option<Vec<String>>,
    pub results: Vec<RFPairResult>,
}

#[derive(Debug, Serialize)]
pub struct RFMultiResult {
    pub nrecords: usize,
    pub multi_type: RFMultiType,
    pub nrf_summary: (f64, f64),
    pub rf_rate_summary: (f64, f64),
    pub fp_rate_summary: (f64, f64),
    pub fn_rate_summary: (f64, f64),
    pub results: Vec<RFPrettyOutput>,
}

impl RFMultiResult {
    pub fn new(multi_type: RFMultiType, results: Vec<RFPrettyOutput>) -> Self {
        let mut nrfs: Vec<f64> = vec![];
        let mut rfs: Vec<f64> = vec![];
        let mut fps: Vec<f64> = vec![];
        let mut fns: Vec<f64> = vec![];
        for result in &results {
            nrfs.push(result.n_rf);
            rfs.push(result.rf_rate);
            fps.push(result.fp_rate);
            fns.push(result.fn_rate);
        }
        let nrf_summary = two_pass_mean_stddev(&nrfs);
        let rf_rate_summary = two_pass_mean_stddev(&rfs);
        let fp_rate_summary = two_pass_mean_stddev(&fps);
        let fn_rate_summary = two_pass_mean_stddev(&fns);
        RFMultiResult {
            nrecords: results.len(),
            multi_type,
            nrf_summary,
            rf_rate_summary,
            fp_rate_summary,
            fn_rate_summary,
            results,
        }
    }
}

fn execute_rf_allpairs(trees: &[String], format: Format) {
    let coll = if trees.len() == 1 {
        TreeCollection::from_newick(&trees[0]).unwrap()
    } else {
        TreeCollection::from_multiple(trees).unwrap()
    };
    assert!(coll.ngenes() > 1);
    let res = ogcat::compare_allpairs(&coll);
    let mut paired_results: Vec<RFPairResult> = vec![];
    let mut cnt = 0;
    for i in 0..(coll.ngenes() - 1) {
        for j in (i + 1)..coll.ngenes() {
            paired_results.push(RFPairResult {
                pair: (i, j),
                result: RFPrettyOutput::new(res[cnt]),
            });
            cnt += 1;
        }
    }
    let output = RFPairResults {
        files: if trees.len() == 1 {
            None
        } else {
            Some(trees.iter().cloned().collect())
        },
        results: paired_results,
    };
    assert_eq!(
        Format::Json,
        format,
        "none JSON formats not yet implemented"
    );
    println!("{}", serde_json::to_string_pretty(&output).unwrap());
}

fn execute_rf(reference: &PathBuf, estimated: &PathBuf, fast: bool, format: Format) {
    let mut trees = TreeCollection::from_newick(reference).unwrap();
    let ntaxa = trees.ntaxa();
    trees.add_trees(estimated).unwrap();
    assert_eq!(
        ntaxa,
        trees.ntaxa(),
        "Number of taxa in reference and estimated trees do not match"
    );
    let res = if !fast {
        ogcat::compare_two_trees_amplified(&trees)
    } else {
        ogcat::compare_two_trees(&trees)
    };
    let pretty = ogcat::RFPrettyOutput::new(res);
    match format {
        Format::Human => {
            let table = Builder::default()
                .set_columns([
                    "ntaxa", "fp_edges", "fn_edges", "n_rf", "rf_rate", "fp_rate", "fn_rate",
                ])
                .add_record([
                    pretty.raw.ntaxa.to_string(),
                    pretty.raw.fp_edges.to_string(),
                    pretty.raw.fn_edges.to_string(),
                    format!("{:.2}%", pretty.n_rf * 100.0),
                    format!("{:.2}%", pretty.rf_rate * 100.0),
                    format!("{:.2}%", pretty.fp_rate * 100.0),
                    format!("{:.2}%", pretty.fn_rate * 100.0),
                ])
                .build()
                .with(Style::modern());
            println!("{}", table);
        }
        Format::Json => {
            println!("{}", serde_json::to_string_pretty(&pretty).unwrap());
        }
    }
}

fn format_uncertainty_percentage(d: &(f64, f64)) -> String {
    let (mean, stddev) = d;
    format!("{:.2}% ± {:.2}%", mean * 100f64, stddev * 100f64)
}

fn execute_rf_multi(reference: &PathBuf, estimated: &PathBuf, format: Format) {
    let mut trees = TreeCollection::from_newick(reference).unwrap();
    let ntaxa = trees.ntaxa();
    let ngenes = trees.ngenes();
    trees.add_trees(estimated).unwrap();
    assert_eq!(
        ntaxa,
        trees.ntaxa(),
        "Number of taxa in reference and estimated trees do not match"
    );

    if ngenes > 1 {
        assert_eq!(
            ngenes * 2,
            trees.ngenes(),
            "Number of genes in reference and estimated trees do not match"
        );
    }

    let rf_outputs = if ngenes > 1 {
        ogcat::compare_many2many(&trees)
    } else {
        ogcat::compare_one2many(&trees)
    };
    let pretty = rf_outputs
        .iter()
        .map(|it| RFPrettyOutput::new(*it))
        .collect::<Vec<_>>();
    let multiresult = RFMultiResult::new(
        if ngenes > 1 {
            RFMultiType::ManyPairs
        } else {
            RFMultiType::OneToMany
        },
        pretty,
    );
    match format {
        Format::Human => {
            let table = Builder::default()
                .set_columns([
                    "nrecords",
                    "multi_type",
                    "n_rf",
                    "rf_rate",
                    "fp_rate",
                    "fn_rate",
                ])
                .add_record([
                    multiresult.nrecords.to_string(),
                    multiresult.multi_type.to_string(),
                    format_uncertainty_percentage(&multiresult.nrf_summary),
                    format_uncertainty_percentage(&multiresult.rf_rate_summary),
                    format_uncertainty_percentage(&multiresult.fp_rate_summary),
                    format_uncertainty_percentage(&multiresult.fn_rate_summary),
                ])
                .build()
                .with(Style::modern());
            println!("{}", table.to_string());
        }
        Format::Json => {
            println!("{}", serde_json::to_string_pretty(&multiresult).unwrap());
        }
    }
}

fn main() {
    let args = Args::parse();
    match args.cmd {
        SubCommand::Rf {
            reference,
            estimated,
            fast,
        } => {
            execute_rf(&reference, &estimated, fast, args.format);
        }
        SubCommand::Mask {
            input,
            output,
            percent,
        } => {
            let res = if percent >= 1f64 {
                aln::aln_mask(&input, 1, 1f64, &output)
            } else {
                aln::aln_mask(&input, 0, percent, &output)
            };
            println!("{}", Table::new([res]).with(Style::modern()).to_string());
        }
        SubCommand::AlnStats { input, pdis } => {
            let res = aln::aln_linear_stats(&input);
            let p_result = if pdis {
                Some(aln::approx_pdis(&input, (res.alph).clone()).unwrap())
            } else {
                None
            };
            let table = Builder::default()
                .set_columns([
                    "alph",
                    "columns",
                    "rows",
                    "gap_ratio",
                    "avg_seq_length",
                    "avg_p_dis",
                    "max_p_dis",
                ])
                .add_record([
                    res.alph.to_string(),
                    res.width.to_string(),
                    res.rows.to_string(),
                    format!(
                        "{:.3}%",
                        (res.gap_cells as f64) / res.total_cells as f64 * 100.0
                    ),
                    format!("{:.4}", res.avg_sequence_length),
                    format!("{:.4}", p_result.as_ref().unwrap().avg_pdis),
                    format!("{:.4}", p_result.unwrap().max_pdis),
                ])
                .build();
            println!("{}", table.with(Style::modern()).to_string());
        }
        SubCommand::RfAllpairs { trees } => {
            execute_rf_allpairs(&trees, args.format);
        }
        SubCommand::RfMulti {
            reference,
            estimated,
        } => {
            execute_rf_multi(&reference, &estimated, args.format);
        }
        _ => {
            panic!("Unsupported subcommand");
        }
    }
}
