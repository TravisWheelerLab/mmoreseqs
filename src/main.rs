mod command_ext;
mod external_steps;
mod pipeline;

use crate::pipeline::{align, prep, search, seed};
use anyhow::Result;
use clap::{ArgAction, Parser, Subcommand};
use std::env;
use std::fs::create_dir_all;
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "mmoreseqs")]
#[command(about = "Some stuff should go here", long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    command: SubCommands,
    /// Path for alignment output
    /// The number of threads to use
    #[arg(long, default_value_t = 1usize)]
    threads: usize,
    /// Allow output files to be overwritten
    #[arg(long, action = ArgAction::SetTrue)]
    allow_overwrite: Option<bool>,
}

const PREP_ABOUT: &str = "\
Prepare a query (MSA) file and target (fasta) file\
";

const PREP_LONG_ABOUT: &str = "\
------------------
| mmoreseqs prep |
------------------

This command is useful for caching the first step in the mmoreseqs pipeline.

For the query:  produce both a P7 profile HMM and an MMseqs2 profile database.
For the target: produce an MMseqs2 sequence database.\
";

const SEED_ABOUT: &str = "\
\
";

const SEED_LONG_ABOUT: &str = "\
------------------
| mmoreseqs seed |
------------------
";

const ALIGN_ABOUT: &str = "\
\
";

const ALIGN_LONG_ABOUT: &str = "\
-------------------
| mmoreseqs align |
-------------------
";

const SEARCH_ABOUT: &str = "\
\
";

const SEARCH_LONG_ABOUT: &str = "\
--------------------
| mmoreseqs search |
--------------------
";

#[derive(Debug, Parser)]
struct CommonArgs {
    /// Path for alignment output
    /// The number of threads to use
    #[arg(long, default_value_t = 1usize)]
    threads: usize,
    /// Allow output files to be overwritten
    #[arg(long, action = ArgAction::SetTrue)]
    allow_overwrite: Option<bool>,
}

/// Doc comment
#[derive(Debug, Subcommand)]
enum SubCommands {
    #[command(about = PREP_ABOUT, long_about = PREP_LONG_ABOUT)]
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
    #[command(about = SEED_ABOUT, long_about = SEED_LONG_ABOUT)]
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
    #[command(about = ALIGN_ABOUT, long_about = ALIGN_LONG_ABOUT)]
    Align {
        /// Query P7 HMM file
        query: String,
        /// Target fasta file
        target: String,
        /// Seed file (result of mmoreseqs seed)
        seeds: String,
        #[command(flatten)]
        common: CommonArgs,
    },
    #[command(about = SEARCH_ABOUT, long_about = SEARCH_LONG_ABOUT)]
    Search {
        /// Query MSA file
        query: String,
        /// Target fasta file
        target: String,
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
                common,
            } => {
                args.set_common(&common);
                args.command = Command::Align;
                args.paths.query_hmm = PathBuf::from(query);
                args.paths.target_fasta = PathBuf::from(target);
                args.paths.seeds = PathBuf::from(seeds);
            }
            SubCommands::Search {
                query,
                target,
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

                args.paths.results = PathBuf::from("results.tsv");
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
    pub allow_overwrite: bool,
}

impl Args {
    fn set_common(&mut self, args: &CommonArgs) {
        self.threads = args.threads;
        self.allow_overwrite = args.allow_overwrite.unwrap_or(false);
    }
}

fn main() -> Result<()> {
    let args = Cli::parse().args();

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
            panic!("somehow failed to set the command")
        }
    }

    Ok(())
}
