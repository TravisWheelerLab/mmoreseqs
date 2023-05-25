use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Write};

use crate::args::{Args, FileFormat};
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
use nale::structs::hmm::parse_hmms_from_p7hmm_file;
use nale::structs::{Alignment, Profile, Sequence, Trace};

use crate::pipeline::{
    align_threaded_a, align_threaded_b, align_threaded_c, align_threaded_d, align_threaded_e,
    align_threaded_f,
};
use anyhow::Context;
use nale::structs::alignment::ScoreParams;
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

pub fn align(
    args: &Args,
    profiles: Option<Vec<Profile>>,
    seed_map: Option<SeedMap>,
) -> anyhow::Result<()> {
    let profiles = match profiles {
        // if we happened to run the seed step before
        // this, the profiles will be passed in
        Some(profiles) => profiles,
        None => {
            let hmms = match args.query_format {
                FileFormat::Fasta => {
                    build_hmm_from_fasta(args)?;
                    parse_hmms_from_p7hmm_file(args.query_hmm().to_str().unwrap())?
                }
                FileFormat::Stockholm => {
                    build_hmm_from_stockholm(args)?;
                    parse_hmms_from_p7hmm_file(args.query_hmm().to_str().unwrap())?
                }
                FileFormat::Hmm => parse_hmms_from_p7hmm_file(args.paths.query.to_str().unwrap())?,
                FileFormat::Unset => {
                    panic!();
                }
            };

            hmms.iter().map(Profile::new).collect()
        }
    };

    let seed_map = match seed_map {
        // if we happened to run the seed step before
        // this, the seeds will be passed in
        Some(seed_map) => seed_map,
        None => {
            let mut seeds_string = String::new();
            File::open(&args.paths.seeds)
                .context(format!(
                    "failed to open alignment seeds file: {}",
                    &args.paths.seeds.to_string_lossy(),
                ))?
                .read_to_string(&mut seeds_string)
                .context(format!(
                    "failed to read alignment seeds file: {}",
                    &args.paths.seeds.to_string_lossy(),
                ))?;

            serde_json::from_str(&seeds_string).context(format!(
                "failed to parse alignment seeds file: {}",
                &args.paths.seeds.to_string_lossy(),
            ))?
        }
    };

    if args.threads == 1 {
        align_serial(args, profiles, seed_map)?;
    } else if args.threads == 2 {
        align_threaded_a(args, profiles, seed_map)?;
    } else if args.threads == 3 {
        align_threaded_b(args, profiles, seed_map)?;
    } else if args.threads == 4 {
        align_threaded_c(args, profiles, seed_map)?;
    } else if args.threads == 5 {
        align_threaded_d(args, profiles, seed_map)?;
    } else if args.threads == 6 {
        align_threaded_e(args, profiles, seed_map)?;
    } else if args.threads == 7 {
        align_threaded_f(args, profiles, seed_map)?;
    }
    Ok(())
}

pub fn align_serial(
    args: &Args,
    mut profiles: Vec<Profile>,
    seed_map: SeedMap,
) -> anyhow::Result<()> {
    let targets = Sequence::amino_from_fasta(&args.paths.target)?;

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

    let mut results_writer = args.paths.results.open(true)?;

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
