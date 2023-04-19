mod extension_traits;
mod pipeline;

use crate::extension_traits::CommandExt;
use crate::pipeline::{align, prep, search, seed};

use std::fs::{create_dir_all, File};
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use thiserror::Error;

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
        /// Where to place the prepared files
        #[arg(short, long, default_value = "./prep/")]
        prep_dir: String,
        #[command(flatten)]
        common: CommonArgs,
    },
    #[command(about = "Use MMseqs2 to create a set of alignment seeds for the align step")]
    Seed {
        /// The directory containing output from running mmoreseqs prep
        prep_dir: String,
        /// Where to place the seeds output
        #[arg(short, long, default_value = "seeds.json")]
        output_file: String,
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
        /// Where to place the results
        #[arg(short, long, default_value = "results.tsv")]
        output_file: String,
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
        #[arg(long, default_value = "./prep/")]
        prep_dir: String,
        #[command(flatten)]
        common: CommonArgs,
    },
}

impl Cli {
    fn args(self) -> Result<Args> {
        let mut args = Args::default();
        match self.command {
            SubCommands::Prep {
                query,
                target,
                prep_dir,
                common,
            } => {
                args.set_common(&common);
                args.command = Command::Prep;

                args.paths.query = PathBuf::from(query);
                args.paths.target = PathBuf::from(target);

                args.paths.prep_dir = PathBuf::from(prep_dir);
                create_dir_all(&args.paths.prep_dir).context(format!(
                    "failed to create prep output directory: {}",
                    args.paths.prep_dir.to_string_lossy()
                ))?;
            }
            SubCommands::Seed {
                prep_dir,
                output_file,
                common,
            } => {
                args.set_common(&common);
                args.command = Command::Seed;

                args.paths.prep_dir = PathBuf::from(prep_dir);
                args.paths.seeds = PathBuf::from(output_file);
            }
            SubCommands::Align {
                query,
                target,
                seeds,
                evalue_cutoff,
                output_file,
                common,
            } => {
                args.set_common(&common);
                args.command = Command::Align;

                args.paths.query = PathBuf::from(query);
                args.paths.target = PathBuf::from(target);
                args.paths.seeds = PathBuf::from(seeds);
                args.evalue_cutoff = evalue_cutoff;
                args.paths.results = PathBuf::from(output_file);
            }
            SubCommands::Search {
                query,
                target,
                evalue_cutoff,
                output_file,
                prep_dir,
                common,
            } => {
                args.set_common(&common);

                args.command = Command::Search;
                args.paths.query = PathBuf::from(query);
                args.paths.target = PathBuf::from(target);

                args.paths.prep_dir = PathBuf::from(prep_dir);
                create_dir_all(&args.paths.prep_dir).context(format!(
                    "failed to create prep output directory: {}",
                    args.paths.prep_dir.to_string_lossy()
                ))?;

                args.evalue_cutoff = evalue_cutoff;
                args.paths.seeds = args.paths.prep_dir.join("seeds.json");
                args.paths.results = PathBuf::from(output_file);
            }
        }
        match args.command {
            Command::Prep | Command::Align | Command::Search => args.guess_query_format()?,
            _ => {}
        }
        Ok(args)
    }
}

#[derive(Default)]
pub struct Paths {
    /// The query file provided at the command line
    pub query: PathBuf,
    /// The target file provided at the command line
    pub target: PathBuf,
    /// The directory under which prepared/intermediate files are placed
    pub prep_dir: PathBuf,
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
    pub paths: Paths,
    pub query_format: FileFormat,
    pub threads: usize,
    pub evalue_cutoff: f32,
}

impl Args {
    fn guess_query_format(&mut self) -> Result<()> {
        let file = File::open(&self.paths.query).context(format!(
            "failed to open query file: {}",
            &self.paths.query.to_string_lossy()
        ))?;

        let mut reader = BufReader::new(file);
        let mut first_line = String::new();
        reader.read_line(&mut first_line)?;

        if &first_line[0..1] == ">" {
            self.query_format = FileFormat::Fasta;
        } else if &first_line[0..11] == "# STOCKHOLM" {
            self.query_format = FileFormat::Stockholm;
        } else if &first_line[0..5] == "HMMER" {
            self.query_format = FileFormat::Hmm;
        } else {
            return Err(UnrecognizedFileFormatError).context(format!(
                "couldn't guess the format of query file: {}",
                &self.paths.query.to_string_lossy()
            ));
        };
        Ok(())
    }

    fn set_common(&mut self, args: &CommonArgs) {
        self.threads = args.threads;
    }
    /// Produce a path to the query P7 HMM file
    fn query_hmm(&self) -> PathBuf {
        self.paths.prep_dir.join("query.hmm")
    }
    /// Produce a path to the MMseqs2 query database.
    ///
    /// If a fasta target was provided, this will be a sequence database.
    /// If a stockholm target was provided, this will be a profile database.
    fn mmseqs_query_db(&self) -> PathBuf {
        self.paths.prep_dir.join("queryDB")
    }
    /// Produce a path to the MMseqs2 query database index
    fn mmseqs_query_db_index(&self) -> PathBuf {
        self.paths.prep_dir.join("queryDB.index")
    }
    /// Produce a path to the MMseqs2 query database h file
    fn mmseqs_query_db_h(&self) -> PathBuf {
        self.paths.prep_dir.join("queryDB_h")
    }
    /// Produce a path to the MMseqs2 query database h file index
    fn mmseqs_query_db_h_index(&self) -> PathBuf {
        self.paths.prep_dir.join("queryDB_h.index")
    }
    /// Produce a path to the MMseqs2 target database.
    ///
    /// This will always be a sequence database.
    fn mmseqs_target_db(&self) -> PathBuf {
        self.paths.prep_dir.join("targetDB")
    }
    /// Produce a path to the MMseqs2 prefilter database.
    ///
    /// This is the result of running `mmseqs prefilter` on the query and target databases.
    fn mmseqs_prefilter_db(&self) -> PathBuf {
        self.paths.prep_dir.join("prefilterDB")
    }
    /// Produce a path to the MMseqs2 alignment database.
    ///
    /// This is the result of running `mmseqs align` on the query, target, and prefilter databases.
    fn mmseqs_align_db(&self) -> PathBuf {
        self.paths.prep_dir.join("alignDB")
    }
    /// Produce a path to the MMseqs2 alignment output.
    ///
    /// This is the result of running `mmseqs convertalis` on the query, target, and align databases.
    fn mmseqs_align_tsv(&self) -> PathBuf {
        self.paths.prep_dir.join("align.tsv")
    }
}

#[derive(Default)]
pub enum FileFormat {
    Fasta,
    Stockholm,
    Hmm,
    #[default]
    Unset,
}

#[derive(Error, Debug)]
#[error("can't guess file format")]
pub struct UnrecognizedFileFormatError;

fn check_hmmer_installed() -> Result<()> {
    std::process::Command::new("hmmbuild")
        .arg("-h")
        .run()
        .context("hmmbuild does not appear to be in the system path")
}

fn check_mmseqs_installed() -> Result<()> {
    std::process::Command::new("mmseqs")
        .arg("-h")
        .run()
        .context("mmseqs2 does not appear to be in the system path")
}

fn main() -> Result<()> {
    let args = Cli::parse().args()?;

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
            align(&args, None, None)?;
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
