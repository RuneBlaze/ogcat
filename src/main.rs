pub mod aln;
pub mod extract;
pub mod fastsp;
pub mod ogtree;
pub mod tree;

use crate::fastsp::pretty_spresults;
use crate::ogtree::*;
use aln::{Approx, CombinedAlnStats};
use anyhow::Ok;
use clap::{ArgEnum, Parser, Subcommand};
use itertools::Itertools;
use ogcat::aln::{aln_diff, aln_diff_summary_table, AlnDiffMode};
use once_cell::sync::Lazy;
use serde::Serialize;
use shadow_rs::shadow;
use std::fmt::Display;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use tabled::Table;
use tabled::{builder::Builder, Style};
shadow!(build);

#[derive(Parser, Debug)]
#[clap(author, version, about, long_version = build::CLAP_LONG_VERSION)]
struct Args {
    #[clap(subcommand)]
    cmd: SubCommand,
    #[clap(long,arg_enum, default_value_t = Format::Human)]
    format: Format,
    #[clap(short, long)]
    output: Option<PathBuf>,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum, Debug)]
enum Format {
    Human,
    Json,
}

#[derive(Subcommand, Debug)]
enum SubCommand {
    /// Compute the RF and related rates of two trees
    Rf {
        #[clap(short, long = "ref")]
        reference: PathBuf,
        #[clap(short, long = "est")]
        estimated: PathBuf,
        #[clap(short, long)]
        fast: bool,
    },

    /// Compute pairwise (one-to-many, many-pairs) RF rates
    RfMulti {
        #[clap(short, long = "ref")]
        reference: PathBuf,
        #[clap(short, long = "est")]
        estimated: PathBuf,
    },

    /// Compute all pairs RF rates from a collection of trees
    RfAllpairs {
        #[clap(short, long = "trees", multiple_values = true)]
        trees: Vec<String>,
    },

    /// Statistics about a Newick unrooted tree (# taxa, degree of resolution, etc.)
    TreeStats {
        #[clap()]
        input: PathBuf,
    },

    /// Centroid-edge decomposition of a tree
    TreeDecomp {
        #[clap()]
        input: PathBuf,
        #[clap(short = 'c', long)]
        max_count: Option<usize>,
        #[clap(short = 's', long)]
        max_size: Option<usize>,
    },

    /// Compute the sum-of-pairs error of two alignments a la FastSP.
    Sp {
        /// Path to reference alignment in FASTA format
        #[clap(short, long = "ref")]
        reference: PathBuf,
        /// Path to estimated alignment in FASTA format
        #[clap(short, long = "est")]
        estimated: PathBuf,
        /// Treat lower-case letters NOT as insertion columns for both reference and estimated
        #[clap(long)]
        ignore_case: bool,
        /// "ignore-case" only for the estimated alignment
        #[clap(long)]
        ignore_case_est: bool,
        /// "ignore-case" only for the reference alignment
        #[clap(long)]
        ignore_case_ref: bool,
        /// Allow comparison when the reference is a subset of the estimated
        #[clap(long)]
        restricted: bool,
        /// Specify an additional character denoting missing data (gap)
        #[clap(short, long = "missing")]
        missing_char: Option<char>,
    },

    /// Statistics about a FASTA alignment (gap ratio, p-distance, etc.)
    AlnStats {
        #[clap()]
        input: PathBuf,
        /// skip the computation of p-distance (avg, max)
        #[clap(long)]
        skip_pdis: bool,
        /// decide if the alignment should be subsampled for computing the p-distance
        #[clap(short, long, arg_enum, default_value_t = Approx::Auto)]
        approx: Approx,
    },

    /// Filter an alignment record-wise
    AlnWhere {
        #[clap()]
        input: PathBuf,
        #[clap(short, long)]
        output: PathBuf,
        #[clap(long = "len_lb")]
        length_lb: Option<usize>,
        #[clap(long = "len_ub")]
        length_ub: Option<usize>,
        #[clap(long = "include", multiple_values = true)]
        inclusion: Vec<PathBuf>,
    },

    /// Extract names or other information from an alignment
    AlnExtract {
        #[clap()]
        input: PathBuf,
        #[clap(short, long, multiple_values = true, default_value = "names")]
        types: Vec<extract::InfoType>,
    },

    /// Compares alignment identity (by default sans gaps)
    AlnDiff {
        #[clap(multiple_values = true)]
        input: Vec<PathBuf>,
        #[clap(short, long)]
        keep_gaps: bool,
        #[clap(long, arg_enum, default_value_t = AlnDiffMode::Udiff)]
        mode: AlnDiffMode,
    },

    /// Mask gappy sites in an alignment
    Mask {
        #[clap()]
        input: PathBuf,
        #[clap(short, long)]
        output: PathBuf,
        #[clap(short, long, default_value_t = 1f64)]
        percent: f64,
        /// Treat lower-case letters NOT as insertion columns
        #[clap(long)]
        ignore_case: bool,
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
    pub result: ogtree::RFPrettyOutput,
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
    let res = compare_allpairs(&coll);
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
            Some(trees.to_vec())
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
        ogtree::compare_two_trees_amplified(&trees)
    } else {
        ogtree::compare_two_trees(&trees)
    };
    let pretty = ogtree::RFPrettyOutput::new(res);
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
        ogtree::compare_many2many(&trees)
    } else {
        ogtree::compare_one2many(&trees)
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
            println!("{}", table);
        }
        Format::Json => {
            println!("{}", serde_json::to_string_pretty(&multiresult).unwrap());
        }
    }
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();
    let mut out = Lazy::new(|| match args.output {
        Some(x) => Box::new(BufWriter::new(File::create(x).unwrap())),
        None => Box::new(io::stdout()) as Box<dyn Write>,
    });
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
            ignore_case,
        } => {
            let res = if percent >= 1f64 {
                aln::aln_mask(&input, 1, 1f64, ignore_case, &output)
            } else {
                aln::aln_mask(&input, 0, percent, ignore_case, &output)
            };
            writeln!(&mut out, "{}", Table::new([res]).with(Style::modern())).unwrap();
        }
        SubCommand::AlnStats {
            input,
            skip_pdis,
            approx,
        } => {
            let (stats, p_result) = aln::aln_linear_stats(&input, !skip_pdis, approx);
            match args.format {
                Format::Human => {
                    let table = Builder::default()
                        .set_columns([
                            "alph",
                            "columns",
                            "rows",
                            "gap_ratio",
                            "avg_seq_length",
                            "avg_p_dis",
                            "max_p_dis",
                            "p_dis_approx",
                        ])
                        .add_record([
                            stats.alph.to_string(),
                            stats.width.to_string(),
                            stats.rows.to_string(),
                            format!(
                                "{:.3}%",
                                (stats.gap_cells as f64) / stats.total_cells as f64 * 100.0
                            ),
                            format!("{:.4}", stats.avg_sequence_length),
                            p_result.as_ref().map_or("Skipped".to_string(), |res| {
                                format!("{:.4}", res.avg_pdis)
                            }),
                            p_result.as_ref().map_or("Skipped".to_string(), |res| {
                                format!("{:.4}", res.max_pdis)
                            }),
                            p_result
                                .as_ref()
                                .map_or("Skipped".to_string(), |res| format!("{}", res.approx)),
                        ])
                        .build();
                    writeln!(&mut out, "{}", table.with(Style::modern())).unwrap();
                }
                Format::Json => {
                    let formatted_stats = CombinedAlnStats::new(&stats, &p_result);
                    writeln!(
                        &mut out,
                        "{}",
                        serde_json::to_string_pretty(&formatted_stats).unwrap()
                    )
                    .unwrap();
                }
            }
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
        SubCommand::AlnWhere {
            input,
            output,
            length_lb,
            length_ub,
            inclusion,
        } => {
            let res = aln::aln_where(&input, length_lb, length_ub, &inclusion, &output);
            match args.format {
                Format::Human => {
                    writeln!(&mut out, "{}", Table::new([res]).with(Style::modern())).unwrap();
                }
                Format::Json => {
                    writeln!(&mut out, "{}", serde_json::to_string_pretty(&res).unwrap()).unwrap();
                }
            }
        }
        SubCommand::TreeStats { input } => {
            let collection = TreeCollection::from_newick(input).unwrap();
            let stats = tree::tree_stats(&collection);
            match args.format {
                Format::Human => {
                    writeln!(&mut out, "{}", Table::new(&stats).with(Style::modern())).unwrap();
                }
                Format::Json => {
                    writeln!(
                        &mut out,
                        "{}",
                        serde_json::to_string_pretty(&stats).unwrap()
                    )
                    .unwrap();
                }
            }
        }
        SubCommand::TreeDecomp {
            input,
            max_count,
            max_size,
        } => {
            assert_eq!(
                Format::Json,
                args.format,
                "none JSON formats not yet implemented"
            );
            let collection = TreeCollection::from_newick(input).unwrap();
            let decomp = ogtree::centroid_edge_decomp(&collection.trees[0], &max_count, &max_size);
            let labels = ogtree::cuts_to_subsets(&collection.trees[0], &decomp);
            let mut humanized_clusters: Vec<Vec<&str>> = vec![vec![]; decomp.len() as usize];
            let ts = &collection.taxon_set;
            for i in 0..ts.len() {
                humanized_clusters[labels[i]].push(&ts.names[i]);
            }
            let json = serde_json::to_string_pretty(&humanized_clusters).unwrap();
            writeln!(&mut out, "{}", json).unwrap();
        }
        SubCommand::Sp {
            reference,
            estimated,
            ignore_case,
            ignore_case_est,
            ignore_case_ref,
            restricted,
            missing_char,
        } => {
            let mut b = [0; 1];
            let char = missing_char.map(|it| {
                it.to_ascii_uppercase().encode_utf8(&mut b);
                b[0]
            });
            let ig_est = ignore_case || ignore_case_est;
            let ig_ref = ignore_case || ignore_case_ref;
            let fastsp_result = if restricted {
                let oppo =
                    fastsp::calc_fpfn(&estimated, &reference, ig_ref, ig_est, restricted, char);
                oppo.flip()
            } else {
                fastsp::calc_fpfn(&reference, &estimated, ig_est, ig_ref, restricted, char)
            };
            match args.format {
                Format::Human => {
                    writeln!(&mut out, "{}", pretty_spresults(&fastsp_result)).unwrap();
                }
                Format::Json => {
                    writeln!(
                        &mut out,
                        "{}",
                        serde_json::to_string_pretty(&fastsp_result).unwrap()
                    )
                    .unwrap();
                }
            }
        }
        SubCommand::AlnExtract { input, types } => {
            assert_eq!(
                Format::Json,
                args.format,
                "none JSON formats not yet implemented"
            );
            let res = extract::aln_extract(&input, &types);
            match args.format {
                Format::Human => {
                    unimplemented!();
                }
                Format::Json => {
                    writeln!(&mut out, "{}", serde_json::to_string_pretty(&res).unwrap()).unwrap();
                }
            }
        }
        SubCommand::AlnDiff {
            input,
            keep_gaps,
            mode,
        } => {
            assert_eq!(input.len(), 2, "only two files are currently supported");
            for (l, r) in input.iter().tuple_combinations() {
                let r = aln_diff(l, r, !keep_gaps, mode)?;
                if let Some(s) = r {
                    match args.format {
                        Format::Human => {
                            writeln!(&mut out, "{}", aln_diff_summary_table(&s))?;
                        }
                        Format::Json => {
                            writeln!(&mut out, "{}", serde_json::to_string_pretty(&s).unwrap())?;
                        }
                    }
                }
            }
        }
        _ => {
            panic!("Unsupported subcommand");
        }
    }
    Ok(())
}
