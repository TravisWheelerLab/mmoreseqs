mod command_ext;
mod external_steps;
mod pipeline;

use crate::external_steps::{check_hmmer_installed, check_mmseqs_installed};
use crate::pipeline::{align, prep, search, seed};
use anyhow::Result;
use clap::{Parser, Subcommand};
use std::fs::create_dir_all;
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "mmoreseqs")]
#[command(
    about = "Using MMseqs2 to find rough alignment seeds, perform bounded profile HMM sequence alignment"
)]
pub struct Cli {
    #[command(subcommand)]
    command: SubCommands,
}

#[derive(Debug, Parser)]
struct CommonArgs {
    /// Path for alignment output
    /// The number of threads to use
    #[arg(long, default_value_t = 1usize)]
    threads: usize,
}

/// Doc comment
#[derive(Debug, Subcommand)]
enum SubCommands {
    #[command(about = "Prepare a query (MSA) file and target (fasta) file for the seed step")]
    Prep {
        /// Query MSA file
        query: String,
        /// Target fasta file
        target: String,
        /// Where to place output files
        #[arg(short, long, default_value = "./prep/")]
        output_dir: String,
        #[command(flatten)]
        common: CommonArgs,
    },
    #[command(about = "Use MMseqs2 to create a set of alignment seeds for the align step")]
    Seed {
        /// Query MMseqs2 profile database
        query_db: String,
        /// Query P7 profile HMM
        query_hmm: String,
        /// Target MMseqs 2 sequence database
        target: String,
        /// Where to place the seeds output
        #[arg(short, long, default_value = "seeds.tsv")]
        output_file: String,
        /// Where to place intermediate files
        #[arg(short, long, default_value = "./tmp/")]
        work_dir: String,
        #[command(flatten)]
        common: CommonArgs,
    },
    #[command(
        about = "Search with the query (HMM) against the target (fasta), using alignment seeds"
    )]
    Align {
        /// Query P7 HMM file
        query: String,
        /// Target fasta file
        target: String,
        /// Seed file (result of mmoreseqs seed)
        seeds: String,
        /// Only report hits with an E-value above this value
        #[arg(short = 'E', default_value_t = 10.0)]
        evalue_cutoff: f32,
        #[command(flatten)]
        common: CommonArgs,
    },
    #[command(about = "Search a query (MSA) file and target (fasta) file")]
    Search {
        /// Query MSA file
        query: String,
        /// Target fasta file
        target: String,
        /// Only report hits with an E-value above this value
        #[arg(short = 'E', default_value_t = 10.0)]
        evalue_cutoff: f32,
        /// Where to place the results
        #[arg(short, long, default_value = "results.tsv")]
        output_file: String,
        /// Where to place intermediate files
        #[arg(long, default_value = "./tmp/")]
        work_dir: String,
        #[command(flatten)]
        common: CommonArgs,
    },
}

impl Cli {
    fn args(self) -> Args {
        let mut args = Args::default();
        match self.command {
            SubCommands::Prep {
                query,
                target,
                output_dir,
                common,
            } => {
                args.set_common(&common);
                args.command = Command::Prep;
                args.paths.query_msa = PathBuf::from(query);
                args.paths.target_fasta = PathBuf::from(target);

                let output_dir = PathBuf::from(output_dir);

                create_dir_all(&output_dir).expect("failed to create output directory");
                args.paths.query_msa_db = output_dir.join("msaDB");
                args.paths.query_db = output_dir.join("queryDB");
                args.paths.target_db = output_dir.join("targetDB");
                args.paths.query_hmm = output_dir.join("query.hmm");
            }
            SubCommands::Seed {
                query_db,
                query_hmm,
                target,
                output_file,
                work_dir,
                common,
            } => {
                args.set_common(&common);
                args.command = Command::Seed;

                args.paths.query_db = PathBuf::from(&query_db);
                args.paths.query_db_index = PathBuf::from(format!("{}.index", query_db));
                args.paths.query_db_h = PathBuf::from(format!("{}_h", query_db));
                args.paths.query_db_h_index = PathBuf::from(format!("{}_h.index", query_db));
                args.paths.query_hmm = PathBuf::from(query_hmm);
                args.paths.target_db = PathBuf::from(target);

                let work_dir = PathBuf::from(work_dir);
                create_dir_all(&work_dir).expect("failed to create working directory");

                args.paths.prefilter_db = work_dir.join("prefilterDB");
                args.paths.align_db = work_dir.join("alignDB");

                args.paths.seeds = PathBuf::from(output_file);
            }
            SubCommands::Align {
                query,
                target,
                seeds,
                evalue_cutoff,
                common,
            } => {
                args.set_common(&common);
                args.command = Command::Align;
                args.paths.query_hmm = PathBuf::from(query);
                args.paths.target_fasta = PathBuf::from(target);
                args.paths.seeds = PathBuf::from(seeds);
                args.evalue_cutoff = evalue_cutoff;
            }
            SubCommands::Search {
                query,
                target,
                evalue_cutoff,
                output_file,
                work_dir,
                common,
            } => {
                args.set_common(&common);

                args.command = Command::Search;
                args.paths.query_msa = PathBuf::from(query);
                args.paths.target_fasta = PathBuf::from(target);

                let work_dir = PathBuf::from(work_dir);

                create_dir_all(&work_dir).expect("failed to create working directory");

                args.paths.query_msa_db = work_dir.join("msaDB");
                args.paths.query_db = work_dir.join("queryDB");
                args.paths.query_db_index = work_dir.join("queryDB.index");
                args.paths.query_db_h = work_dir.join("queryDB_h");
                args.paths.query_db_h_index = work_dir.join("queryDB_h.index");
                args.paths.target_db = work_dir.join("targetDB");
                args.paths.prefilter_db = work_dir.join("prefilterDB");
                args.paths.align_db = work_dir.join("alignDB");
                args.paths.seeds = work_dir.join("seeds.tsv");
                args.paths.query_hmm = work_dir.join("query.hmm");

                args.evalue_cutoff = evalue_cutoff;
                args.paths.results = PathBuf::from(output_file);
            }
        }
        args
    }
}

#[derive(Default)]
pub struct FilePaths {
    pub query_hmm: PathBuf,
    pub query_msa: PathBuf,
    pub target_fasta: PathBuf,
    pub query_msa_db: PathBuf,
    pub query_db: PathBuf,
    pub query_db_index: PathBuf,
    pub query_db_h: PathBuf,
    pub query_db_h_index: PathBuf,
    pub target_db: PathBuf,
    pub prefilter_db: PathBuf,
    pub align_db: PathBuf,
    pub seeds: PathBuf,
    pub results: PathBuf,
}

#[derive(Default)]
pub enum Command {
    Prep,
    Seed,
    Align,
    Search,
    #[default]
    CommandNotSet,
}

#[derive(Default)]
pub struct Args {
    pub command: Command,
    pub paths: FilePaths,
    pub threads: usize,
    pub evalue_cutoff: f32,
}

impl Args {
    fn set_common(&mut self, args: &CommonArgs) {
        self.threads = args.threads;
    }
}

fn main() -> Result<()> {
    let args = Cli::parse().args();

    check_hmmer_installed()?;
    check_mmseqs_installed()?;

    match args.command {
        Command::Prep => {
            prep(&args)?;
        }
        Command::Seed => {
            seed(&args)?;
        }
        Command::Align => {
            align(&args)?;
        }
        Command::Search => {
            search(&args)?;
        }
        Command::CommandNotSet => {
            unreachable!()
        }
    }

    Ok(())
}
