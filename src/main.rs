mod args;
mod cli;
mod extension_traits;
mod pipeline;

use args::Command;
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
        Command::NotSet => {
            unreachable!()
        }
    }

    Ok(())
}
