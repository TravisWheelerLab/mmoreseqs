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
        #[arg(value_name = "QUERY.[fasta:sto]")]
        query: String,
        /// Target file
        #[arg(value_name = "TARGET.fasta")]
        target: String,
        /// Only report hits with an E-value above this value
        #[arg(short = 'E', default_value_t = 10.0)]
        evalue_cutoff: f64,
        /// Where to place the results
        #[arg(short, default_value = "results.tsv")]
        output_file: String,
        /// Where to place intermediate files
        #[arg(short, default_value = "./prep/")]
        prep_dir: String,
        /// The number of threads to use
        #[arg(short, long, default_value_t = 8usize, value_name = "n")]
        threads: usize,
        /// MMseqs2 prefilter: k-mer length (0: automatically set to optimum)
        #[arg(long, default_value_t = 0usize)]
        mmseqs_k: usize,
        /// MMseqs2 prefilter: k-mer threshold for generating similar k-mer lists
        #[arg(long, default_value_t = 80usize)]
        mmseqs_k_score: usize,
        /// MMseqs2 prefilter: Accept only matches with ungapped alignment score above threshold
        #[arg(long, default_value_t = 15usize)]
        mmseqs_min_ungapped_score: usize,
        /// MMseqs2 prefilter: Maximum results per query sequence allowed to pass the prefilter
        #[arg(long, default_value_t = 1000usize)]
        mmseqs_max_seqs: usize,
        /// MMseqs2 align: Include matches below this E-value as seeds
        #[arg(long, default_value_t = 1000f64)]
        mmseqs_e: f64,
    },
    #[command(about = "Prepare a query (MSA) file and target (fasta) file for the seed step")]
    Prep {
        /// Query file
        #[arg(value_name = "QUERY.[fasta:sto]")]
        query: String,
        /// Target file
        #[arg(value_name = "TARGET.fasta")]
        target: String,
        /// Where to place the prepared files
        #[arg(short, long, default_value = "./prep/")]
        prep_dir: String,
        /// The number of threads to use
        #[arg(short, long, default_value_t = 8usize, value_name = "n")]
        threads: usize,
        /// Don't build a profile HMM with the input MSA
        #[arg(long, action)]
        skip_hmmbuild: bool,
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
        /// MMseqs2 prefilter: k-mer length (0: automatically set to optimum)
        #[arg(long, default_value_t = 0usize)]
        mmseqs_k: usize,
        /// MMseqs2 prefilter: k-mer threshold for generating similar k-mer lists
        #[arg(long, default_value_t = 80usize)]
        mmseqs_k_score: usize,
        /// MMseqs2 prefilter: Accept only matches with ungapped alignment score above threshold
        #[arg(long, default_value_t = 15usize)]
        mmseqs_min_ungapped_score: usize,
        /// MMseqs2 prefilter: Maximum results per query sequence allowed to pass the prefilter
        #[arg(long, default_value_t = 1000usize)]
        mmseqs_max_seqs: usize,
        /// MMseqs2 align: Include matches below this E-value as seeds
        #[arg(long, default_value_t = 1000f64)]
        mmseqs_e: f64,
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
        evalue_cutoff: f64,
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
                skip_hmmbuild,
            } => {
                args.threads = threads;
                args.command = Command::Prep;
                args.build_hmm = !skip_hmmbuild;

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
                mmseqs_k,
                mmseqs_k_score,
                mmseqs_min_ungapped_score,
                mmseqs_max_seqs,
                mmseqs_e,
            } => {
                args.threads = threads;
                args.command = Command::Seed;
                args.mmseqs_args.k = mmseqs_k;
                args.mmseqs_args.k_score = mmseqs_k_score;
                args.mmseqs_args.min_ungapped_score = mmseqs_min_ungapped_score;
                args.mmseqs_args.max_seqs = mmseqs_max_seqs;
                args.mmseqs_args.e = mmseqs_e;

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
                mmseqs_k,
                mmseqs_k_score,
                mmseqs_min_ungapped_score,
                mmseqs_max_seqs,
                mmseqs_e,
            } => {
                args.threads = threads;
                args.mmseqs_args.k = mmseqs_k;
                args.mmseqs_args.k_score = mmseqs_k_score;
                args.mmseqs_args.min_ungapped_score = mmseqs_min_ungapped_score;
                args.mmseqs_args.max_seqs = mmseqs_max_seqs;
                args.mmseqs_args.e = mmseqs_e;

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
