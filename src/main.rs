mod args;
mod cli;
mod extension_traits;
mod pipeline;

use args::MmoreCommand;
use cli::Cli;
use extension_traits::CommandExt;
use pipeline::{align, prep, search, seed};

use anyhow::{Context, Result};
use clap::Parser;

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
        MmoreCommand::Prep => {
            prep(&args)?;
        }
        MmoreCommand::Seed => {
            seed(&args)?;
        }
        MmoreCommand::Align => {
            align(&args, None, None)?;
        }
        MmoreCommand::Search => {
            search(&args)?;
        }
        MmoreCommand::NotSet => {
            unreachable!()
        }
    }

    Ok(())
}
