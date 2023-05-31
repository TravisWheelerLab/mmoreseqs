use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Write};
use std::path::PathBuf;

use crate::args::{guess_query_format_from_query_file, FileFormat};
use crate::extension_traits::PathBufExt;
use crate::pipeline::prep::{build_hmm_from_fasta, build_hmm_from_stockholm};
use crate::pipeline::seed::SeedMap;

use nale::align::bounded::structs::{
    CloudBoundGroup, CloudMatrixLinear, CloudSearchParams, DpMatrixSparse, RowBounds,
};
use nale::align::bounded::{
    backward_bounded, cloud_search_backward, cloud_search_forward, forward_bounded, null1_score,
    null2_score, optimal_accuracy_bounded, posterior_bounded, traceback_bounded,
};
use nale::structs::alignment::ScoreParams;
use nale::structs::hmm::parse_hmms_from_p7hmm_file;
use nale::structs::{Alignment, Profile, Sequence, Trace};

use crate::pipeline::{
    align_threaded_a, align_threaded_b, align_threaded_c, align_threaded_d, align_threaded_e,
    align_threaded_f,
};

use anyhow::Context;
use clap::Args;
use thiserror::Error;

#[derive(Error, Debug)]
#[error("no profile with name: {profile_name}")]
pub struct ProfileNotFoundError {
    profile_name: String,
}

#[derive(Error, Debug)]
#[error("no target with name: {target_name}")]
pub struct TargetNotFoundError {
    target_name: String,
}

#[derive(Debug, Args)]
pub struct AlignArgs {
    /// Query file
    #[arg(value_name = "QUERY.[fasta:hmm:sto]")]
    pub query_path: PathBuf,
    /// Target file
    #[arg(value_name = "TARGET.fasta")]
    pub target_path: PathBuf,
    /// Alignment seeds from running mmoreseqs seed (or elsewhere)
    #[arg(value_name = "SEEDS.json")]
    pub seeds_path: PathBuf,
    /// Only report hits with an E-value above this value
    #[arg(short = 'E', default_value_t = 10.0)]
    pub evalue_cutoff: f64,
    /// Where to place tabular output
    #[arg(short, long, default_value = "results.tsv")]
    pub tsv_results_path: PathBuf,
    /// Where to place alignment output
    #[arg(short, long, default_value = "results.out")]
    pub ali_results_path: PathBuf,
    /// The number of threads to use
    #[arg(short, long, default_value_t = 8usize, value_name = "n")]
    pub num_threads: usize,
    // TODO: this is a temporary dev argument
    /// Thread strategy
    #[arg(short, long, default_value_t = 1usize, value_name = "n")]
    pub thread_strategy: usize,
}

pub fn align(
    args: &AlignArgs,
    profiles: Option<Vec<Profile>>,
    seed_map: Option<SeedMap>,
) -> anyhow::Result<()> {
    let profiles = match profiles {
        // if we happened to run the seed step before
        // this, the profiles will be passed in
        Some(profiles) => profiles,
        None => {
            let query_format = guess_query_format_from_query_file(&args.query_path)?;
            let hmm_path = match query_format {
                FileFormat::Fasta => {
                    let hmm_path = args.query_path.with_extension("hmm");
                    build_hmm_from_fasta(&args.query_path, &hmm_path, args.num_threads)?;
                    hmm_path
                }
                FileFormat::Stockholm => {
                    let hmm_path = args.query_path.with_extension("hmm");
                    build_hmm_from_stockholm(&args.query_path, &hmm_path, args.num_threads)?;
                    hmm_path
                }
                FileFormat::Hmm => args.query_path.clone(),
                FileFormat::Unset => {
                    // TODO: real error
                    panic!("query format is unset in call to align()");
                }
            };

            let hmms = parse_hmms_from_p7hmm_file(hmm_path)?;

            hmms.iter().map(Profile::new).collect()
        }
    };

    let seed_map = match seed_map {
        // if we happened to run the seed step before
        // this, the seeds will be passed in
        Some(seed_map) => seed_map,
        None => {
            let mut seeds_string = String::new();
            File::open(&args.seeds_path)
                .context(format!(
                    "failed to open alignment seeds file: {}",
                    &args.seeds_path.to_string_lossy(),
                ))?
                .read_to_string(&mut seeds_string)
                .context(format!(
                    "failed to read alignment seeds file: {}",
                    &args.seeds_path.to_string_lossy(),
                ))?;

            serde_json::from_str(&seeds_string).context(format!(
                "failed to parse alignment seeds file: {}",
                &args.seeds_path.to_string_lossy(),
            ))?
        }
    };

    let targets = Sequence::amino_from_fasta(&args.target_path)?;

    if args.num_threads == 1 {
        align_serial(args, profiles, targets, seed_map)?;
    } else {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.num_threads)
            .build_global()
            .unwrap();

        match args.thread_strategy {
            1 => align_threaded_a(args, profiles, targets, seed_map)?,
            2 => align_threaded_b(args, profiles, targets, seed_map)?,
            3 => align_threaded_c(args, profiles, targets, seed_map)?,
            4 => align_threaded_d(args, profiles, targets, seed_map)?,
            5 => align_threaded_e(args, profiles, targets, seed_map)?,
            6 => align_threaded_f(args, profiles, targets, seed_map)?,
            _ => {
                panic!("bad thread strategy")
            }
        }
    }

    Ok(())
}

pub fn align_serial(
    args: &AlignArgs,
    mut profiles: Vec<Profile>,
    targets: Vec<Sequence>,
    seed_map: SeedMap,
) -> anyhow::Result<()> {
    let mut score_params = ScoreParams::new(targets.len());

    let mut target_map: HashMap<String, Sequence> = HashMap::new();
    for target in targets {
        target_map.insert(target.name.clone(), target);
    }

    let max_profile_length = profiles
        .iter()
        .fold(0usize, |acc: usize, p: &Profile| acc.max(p.length));

    let max_target_length = target_map
        .values()
        .fold(0usize, |acc: usize, s: &Sequence| acc.max(s.length));

    let mut cloud_matrix = CloudMatrixLinear::new(max_profile_length);

    let mut forward_bounds = CloudBoundGroup::new(max_target_length, max_profile_length);
    let mut backward_bounds = CloudBoundGroup::new(max_target_length, max_profile_length);

    let mut forward_matrix =
        DpMatrixSparse::new(max_target_length, max_profile_length, &RowBounds::default());
    let mut backward_matrix =
        DpMatrixSparse::new(max_target_length, max_profile_length, &RowBounds::default());
    let mut posterior_matrix =
        DpMatrixSparse::new(max_target_length, max_profile_length, &RowBounds::default());
    let mut optimal_matrix =
        DpMatrixSparse::new(max_target_length, max_profile_length, &RowBounds::default());

    let mut results_writer = args.tsv_results_path.open(true)?;

    for profile in profiles.iter_mut() {
        let seeds = match seed_map.get(&profile.name) {
            Some(seeds) => seeds,
            None => {
                continue;
            }
        };

        for seed in seeds {
            let target =
                target_map
                    .get(&seed.target_name[..])
                    .ok_or_else(|| TargetNotFoundError {
                        target_name: seed.target_name.clone(),
                    })?;

            profile.configure_for_target_length(target.length);

            cloud_matrix.reuse(profile.length);
            forward_bounds.reuse(target.length, profile.length);
            backward_bounds.reuse(target.length, profile.length);

            cloud_search_forward(
                profile,
                target,
                seed,
                &mut cloud_matrix,
                &CloudSearchParams::default(),
                &mut forward_bounds,
            );

            cloud_search_backward(
                profile,
                target,
                seed,
                &mut cloud_matrix,
                &CloudSearchParams::default(),
                &mut backward_bounds,
            );

            CloudBoundGroup::join_bounds(&mut forward_bounds, &backward_bounds);

            if !forward_bounds.valid() {
                println!("cloud bound fail");
                continue;
            }

            forward_bounds.trim_wings();

            let row_bounds = RowBounds::new(&forward_bounds);

            if !row_bounds.valid() {
                println!("row bound fail");
                continue;
            }

            forward_matrix.reuse(target.length, profile.length, &row_bounds);
            backward_matrix.reuse(target.length, profile.length, &row_bounds);
            posterior_matrix.reuse(target.length, profile.length, &row_bounds);
            optimal_matrix.reuse(target.length, profile.length, &row_bounds);

            // we use the forward score to compute the final bit score (later)
            score_params.forward_score_nats =
                forward_bounded(profile, target, &mut forward_matrix, &row_bounds);

            backward_bounded(profile, target, &mut backward_matrix, &row_bounds);

            posterior_bounded(
                profile,
                &forward_matrix,
                &backward_matrix,
                &mut posterior_matrix,
                &row_bounds,
            );

            optimal_accuracy_bounded(profile, &posterior_matrix, &mut optimal_matrix, &row_bounds);

            let mut trace = Trace::new(target.length, profile.length);
            traceback_bounded(
                profile,
                &posterior_matrix,
                &optimal_matrix,
                &mut trace,
                row_bounds.target_end,
            );

            score_params.null_score_nats = null1_score(target.length);
            score_params.bias_correction_score_nats =
                null2_score(&posterior_matrix, profile, target, &row_bounds);

            let alignment = Alignment::from_trace(&trace, profile, target, &score_params);

            if alignment.evalue <= args.evalue_cutoff {
                writeln!(results_writer, "{}", alignment.tab_string());
            }
        }
    }
    Ok(())
}
