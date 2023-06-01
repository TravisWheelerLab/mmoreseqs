use crate::extension_traits::PathBufExt;
use crate::pipeline::{align, prep, seed, AlignArgs, MmseqsArgs, PrepArgs, SeedArgs};
use clap::Args;
use std::path::PathBuf;

#[derive(Args)]
pub struct SearchArgs {
    /// Query file
    #[arg(value_name = "QUERY.[fasta:sto]")]
    pub query_path: PathBuf,
    /// Target file
    #[arg(value_name = "TARGET.fasta")]
    pub target_path: PathBuf,
    /// The path to a pre-built P7HMM file
    #[arg(short = 'q', long = "query-hmm", value_name = "QUERY.hmm")]
    pub prebuilt_query_hmm_path: Option<PathBuf>,
    /// Only report hits with an E-value above this value
    #[arg(short = 'E', default_value_t = 10.0)]
    pub evalue_threshold: f64,
    /// Where to place tabular output
    #[arg(short = 'T', long = "tab_output", default_value = "results.tsv")]
    pub tsv_results_path: PathBuf,
    /// Where to place alignment output
    #[arg(short = 'O', long = "output")]
    pub ali_results_path: Option<PathBuf>,
    /// Where to place intermediate files
    #[arg(short, default_value = "./prep/")]
    pub prep_dir_path: PathBuf,
    /// The number of threads to use
    #[arg(
        short = 't',
        long = "threads",
        default_value_t = 8usize,
        value_name = "n"
    )]
    pub num_threads: usize,
    #[command(flatten)]
    pub mmseqs_args: MmseqsArgs,
}

pub fn search(args: &SearchArgs) -> anyhow::Result<()> {
    {
        // quickly make sure we can write the results
        args.tsv_results_path.open(true)?;
    }

    let seeds_path = args.prep_dir_path.join("./seeds.json");
    let prep_args = PrepArgs {
        query_path: args.query_path.clone(),
        target_path: args.target_path.clone(),
        prep_dir_path: args.prep_dir_path.clone(),
        num_threads: args.num_threads,
        skip_hmmbuild: args.prebuilt_query_hmm_path.is_some(),
    };

    let seed_args = SeedArgs {
        prep_dir_path: args.prep_dir_path.clone(),
        seeds_path: seeds_path.clone(),
        prebuilt_query_hmm_path: args.prebuilt_query_hmm_path.clone(),
        num_threads: args.num_threads,
        mmseqs_args: args.mmseqs_args.clone(),
    };

    let align_args = AlignArgs {
        query_path: args.query_path.clone(),
        target_path: args.target_path.clone(),
        seeds_path,
        evalue_threshold: args.evalue_threshold,
        tsv_results_path: args.tsv_results_path.clone(),
        ali_results_path: args.ali_results_path.clone(),
        num_threads: args.num_threads,
    };

    prep(&prep_args)?;
    let (profiles, seed_map) = seed(&seed_args)?;

    align(&align_args, Some(profiles), Some(seed_map))?;
    Ok(())
}
