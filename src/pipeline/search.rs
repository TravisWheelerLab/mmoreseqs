use crate::extension_traits::PathBufExt;
use crate::pipeline::{align, prep, seed};
use crate::Args;

pub fn search(args: &Args) -> anyhow::Result<()> {
    {
        // quickly make sure we can write the results
        args.paths.results.open(true)?;
    }
    prep(args)?;
    let (profiles, seed_map) = seed::seed(args)?;
    // align(args, Some(profiles))?;
    Ok(())
}
