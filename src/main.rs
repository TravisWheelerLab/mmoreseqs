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
        #[command(flatten)]
        common: CommonArgs,
    },
    #[command(about = SEED_ABOUT, long_about = SEED_LONG_ABOUT)]
    Seed {
        /// Query MMseqs2 profile database
        query: String,
        /// Target MMseqs 2 sequence database
        target: String,
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
                common,
            } => {
                args.command = Command::Prep;
                args.set_common(&common);
                args.paths.query_msa = PathBuf::from(query);
                args.paths.target_fasta = PathBuf::from(target);
            }
            SubCommands::Seed {
                query,
                target,
                common,
            } => {
                args.command = Command::Seed;
                args.set_common(&common);
                args.paths.query_db = PathBuf::from(&query);
                args.paths.query_db_index = PathBuf::from(format!("{}.index", query));
                args.paths.query_db_h = PathBuf::from(format!("{}_h", query));
                args.paths.query_db_h_index = PathBuf::from(format!("{}_h.index", query));
                args.paths.target_db = PathBuf::from(target);
            }
            SubCommands::Align {
                query,
                target,
                seeds,
                common,
            } => {
                args.command = Command::Align;
                args.set_common(&common);
                args.paths.query_hmm = PathBuf::from(query);
                args.paths.target_fasta = PathBuf::from(target);
                args.paths.seeds = PathBuf::from(seeds);
            }
            SubCommands::Search {
                query,
                target,
                common,
            } => {
                args.command = Command::Search;
                args.set_common(&common);
                args.paths.query_msa = PathBuf::from(query);
                args.paths.target_fasta = PathBuf::from(target);
            }
        }
        args
    }
}

pub struct FilePaths {
    pub root: PathBuf,
    // prep
    pub query_msa: PathBuf,
    pub target_fasta: PathBuf,
    // seed
    pub query_msa_db: PathBuf,
    pub query_db: PathBuf,
    pub query_db_index: PathBuf,
    pub query_db_h: PathBuf,
    pub query_db_h_index: PathBuf,
    pub target_db: PathBuf,
    pub prefilter_db: PathBuf,
    pub align_db: PathBuf,
    // align
    pub query_hmm: PathBuf,
    pub seeds: PathBuf,
}

impl Default for FilePaths {
    fn default() -> Self {
        let root_path = env::current_dir().unwrap().join("tmp");
        Self {
            root: env::current_dir().unwrap().join("tmp"),
            query_msa: Default::default(),
            target_fasta: Default::default(),
            query_msa_db: root_path.join("msaDB"),
            query_db: root_path.join("queryDB"),
            query_db_index: root_path.join("queryDB.index"),
            query_db_h: root_path.join("queryDB_h"),
            query_db_h_index: root_path.join("queryDB_h.index"),
            target_db: root_path.join("targetDB"),
            prefilter_db: root_path.join("prefilterDB"),
            align_db: root_path.join("alignDB"),
            seeds: root_path.join("seeds.tsv"),
            query_hmm: root_path.join("query.hmm"),
        }
    }
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

    create_dir_all(&args.paths.root)?;

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
