use crate::pipeline::{AlignArgs, PrepArgs, SearchArgs, SeedArgs};
use clap::{Parser, Subcommand};

#[derive(Subcommand)]
pub enum SubCommands {
    #[command(about = "Run the entire mmoreseqs pipeline: prep, seed, & align")]
    Search(SearchArgs),
    #[command(about = "Prepare a query (MSA) file and target (fasta) file for the seed step")]
    Prep(PrepArgs),
    #[command(about = "Use MMseqs2 to create a set of alignment seeds for the align step")]
    Seed(SeedArgs),
    #[command(about = "Search with the query against the target, using alignment seeds")]
    Align(AlignArgs),
}

#[derive(Parser)]
#[command(name = "mmoreseqs")]
#[command(
    about = "Using MMseqs2 to find rough alignment seeds, perform bounded profile HMM sequence alignment"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: SubCommands,
}
