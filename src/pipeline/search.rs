use clap::Args;

#[derive(Args)]
pub struct SearchArgs {
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
}

pub fn search(args: &SearchArgs) -> anyhow::Result<()> {
    todo!();
    // {
    //     // quickly make sure we can write the results
    //     args.paths.results.open(true)?;
    // }
    // prep(args)?;
    // let (profiles, seed_map) = seed::seed(args)?;
    // align(args, Some(profiles), Some(seed_map))?;
    Ok(())
}
