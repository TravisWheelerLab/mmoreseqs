use std::process::Command;

use anyhow::{Context, Result};
use thiserror::Error;

#[derive(Error, Debug)]
#[error("command exited without success")]
struct CommandExitStatusError;

/// An extension trait that is intended to add a run method to the std::process::Command struct.
pub trait CommandExt {
    fn run(&mut self) -> Result<()>;
}

impl CommandExt for Command {
    fn run(&mut self) -> Result<()> {
        let output = self.output().context("failed to start command")?;

        match output.status.success() {
            true => Ok(()),
            false => {
                let stdout = std::str::from_utf8(&output.stdout)
                    .context("failed to convert sdtout to UTF8")?;
                let stderr = std::str::from_utf8(&output.stderr)
                    .context("failed to convert sdterr to UTF8")?;
                println!("stdout: {stdout}");
                println!("stderr: {stderr}");
                Err(CommandExitStatusError.into())
            }
        }
    }
}
