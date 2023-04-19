use crate::args::{Args, Command};

use std::fs::create_dir_all;
use std::path::PathBuf;

use anyhow::Context;
use clap::{Parser, Subcommand};

#[derive(Debug, Subcommand)]
enum SubCommands {
    #[command(about = "Run the entire mmoreseqs pipeline: prep, seed, & align")]
    Search {
        /// Query file
        #[arg(value_name = "QUERY.[fasta:hmm:sto]")]
        query: String,
        /// Target file
        #[arg(value_name = "TARGET.fasta")]
        target: String,
        /// Only report hits with an E-value above this value
        #[arg(short = 'E', default_value_t = 10.0)]
        evalue_cutoff: f32,
        /// Where to place the results
        #[arg(short, default_value = "results.tsv")]
        output_file: String,
        /// Where to place intermediate files
        #[arg(short, default_value = "./prep/")]
        prep_dir: String,
        /// The number of threads to use
        #[arg(short, long, default_value_t = 8usize, value_name = "n")]
        threads: usize,
    },
    #[command(about = "Prepare a query (MSA) file and target (fasta) file for the seed step")]
    Prep {
        /// Query file
        #[arg(value_name = "QUERY.[fasta:hmm:sto]")]
        query: String,
        /// Target file
        #[arg(value_name = "TARGET.fasta")]
        target: String,
        /// Where to place the prepared files
        #[arg(short, long, default_value = "./prep/")]
        prep_dir: String,
        #[arg(short, long, default_value_t = 8usize, value_name = "n")]
        threads: usize,
    },
    #[command(about = "Use MMseqs2 to create a set of alignment seeds for the align step")]
    Seed {
        /// The location of files prepared with mmoreseqs prep
        #[arg()]
        prep_dir: String,
        /// Where to place the seeds output file
        #[arg(short, long, default_value = "seeds.json")]
        seeds: String,
        /// The number of threads to use
        #[arg(short, long, default_value_t = 8usize, value_name = "n")]
        threads: usize,
    },
    #[command(about = "Search with the query against the target, using alignment seeds")]
    Align {
        /// Query file
        #[arg(value_name = "QUERY.[fasta:hmm:sto]")]
        query: String,
        /// Target file
        #[arg(value_name = "TARGET.fasta")]
        target: String,
        /// Alignment seeds from running mmoreseqs seed (or elsewhere)
        #[arg(value_name = "SEEDS.json")]
        seeds: String,
        /// Only report hits with an E-value above this value
        #[arg(short = 'E', default_value_t = 10.0)]
        evalue_cutoff: f32,
        /// Where to place the results
        #[arg(short, long, default_value = "results.tsv")]
        output_file: String,
        /// The number of threads to use
        #[arg(short, long, default_value_t = 8usize, value_name = "n")]
        threads: usize,
    },
}

#[derive(Debug, Parser)]
#[command(name = "mmoreseqs")]
#[command(
    about = "Using MMseqs2 to find rough alignment seeds, perform bounded profile HMM sequence alignment"
)]
pub struct Cli {
    #[command(subcommand)]
    command: SubCommands,
}

impl Cli {
    pub fn args(self) -> anyhow::Result<Args> {
        let mut args = Args::default();
        match self.command {
            SubCommands::Prep {
                query,
                target,
                prep_dir,
                threads,
            } => {
                args.threads = threads;
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
                seeds,
                threads,
            } => {
                args.threads = threads;
                args.command = Command::Seed;

                args.paths.prep_dir = PathBuf::from(prep_dir);
                args.paths.seeds = PathBuf::from(seeds);
            }
            SubCommands::Align {
                query,
                target,
                seeds,
                evalue_cutoff,
                output_file,
                threads,
            } => {
                args.threads = threads;
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
                threads,
            } => {
                args.threads = threads;

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
            Command::Seed => args.get_query_format_from_mmseqs_file()?,
            Command::NotSet => {
                panic!("command not set")
            }
        }
        Ok(args)
    }
}
