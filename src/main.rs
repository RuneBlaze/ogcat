mod ogcat;

use clap::{Parser, Subcommand, ArgEnum};
use ogcat::TreeCollection;
use tabled::{builder::Builder, Style};
use serde::Serialize;
use serde_json::json;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    #[clap(subcommand)]
    cmd: SubCommand,
    #[clap(long,arg_enum, default_value_t = Format::Human)]
    format : Format,
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
}

fn execute_rf(reference: &PathBuf, estimated: &PathBuf, fast: bool, format : Format) {
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
            .set_columns(["ntaxa", "fp_edges", "fn_edges", "n_rf", "rf_rate", "fp_rate", "fn_rate"])
            .add_record([pretty.raw.ntaxa.to_string(), 
            pretty.raw.fp_edges.to_string(),
            pretty.raw.fn_edges.to_string(),
            format!("{:.2}%", pretty.n_rf * 100.0),
            format!("{:.2}%", pretty.rf_rate * 100.0),
            format!("{:.2}%", pretty.fp_rate * 100.0),
            format!("{:.2}%", pretty.fn_rate * 100.0),])
            .build()
            .with(Style::modern());
            println!("{}", table);
        },
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
    }
}
