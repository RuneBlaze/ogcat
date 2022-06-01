mod aln;
mod ogcat;

use clap::{ArgEnum, Parser, Subcommand};
use ogcat::TreeCollection;
use serde::Serialize;
use serde_json::json;
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
                    format!("{:.3}%", res.avg_sequence_length),
                ]);
            // if pdis {

            // } else {

            // }
            // println!("{}", Table::new([&res]).with(Style::modern()).to_string());
            // let pdis = ;
            // println!("{}", Table::new([&pdis]).with(Style::modern()).to_string());
        }
        SubCommand::RfAllpairs { trees } => {
            execute_rf_allpairs(&trees, args.format);
        },
        _ => {
            panic!("Unsupported subcommand");
        }
    }
}
