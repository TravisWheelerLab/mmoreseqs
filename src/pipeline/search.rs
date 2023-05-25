use crate::args::Args;
use crate::extension_traits::PathBufExt;
use crate::pipeline::{align, prep, seed};

pub fn search(args: &Args) -> anyhow::Result<()> {
    {
        // quickly make sure we can write the results
        args.paths.results.open(true)?;
    }
    prep(args)?;
    let (profiles, seed_map) = seed::seed(args)?;
    align(args, Some(profiles), Some(seed_map))?;
    Ok(())
}
