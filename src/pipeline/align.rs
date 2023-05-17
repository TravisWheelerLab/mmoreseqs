use std::collections::HashMap;
use std::fs::File;
use std::io::{stdout, Read};
use std::path::PathBuf;

use crate::args::{Args, FileFormat};
use crate::extension_traits::PathBufExt;
use crate::pipeline::prep::{build_hmm_from_fasta, build_hmm_from_stockholm};
use crate::pipeline::seed::SeedMap;

use nale::align::bounded::structs::{
    CloudBoundGroup, CloudMatrixLinear, CloudSearchParams, DpMatrixSparse, RowBounds, Seed,
};
use nale::align::bounded::{
    backward_bounded, cloud_search_backward, cloud_search_forward, forward_bounded, null1_score,
    null2_score, optimal_accuracy_bounded, posterior_bounded, traceback_bounded,
};
use nale::output::output_standard::write_standard_output;
use nale::output::output_tabular::write_tabular_output;
use nale::structs::dp_matrix::DpMatrix;
use nale::structs::hmm::parse_hmms_from_p7hmm_file;
use nale::structs::{Alignment, Profile, Sequence, Trace};

use anyhow::Context;
use nale::structs::alignment::ScoreParams;
use rayon::prelude::*;
use serde_json::json;
use thiserror::Error;

#[derive(Error, Debug)]
#[error("no profile with name: {profile_name}")]
pub struct ProfileNotFoundError {
    profile_name: String,
}

#[derive(Error, Debug)]
#[error("no seeds for profile: {profile_name}")]
pub struct SeedsNotFoundError {
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

    let mut profile_map: HashMap<String, Profile> = HashMap::new();
    for profile in profiles {
        profile_map.insert(profile.name.clone(), profile);
    }

    let targets = Sequence::amino_from_fasta(&args.paths.target)?;
    let target_count = targets.len();

    let mut target_map: HashMap<String, Sequence> = HashMap::new();
    for target in targets {
        target_map.insert(target.name.clone(), target);
    }

    let max_profile_length = profile_map
        .values()
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

    let mut alignments: Vec<Alignment> = vec![];

    let mut profile_names: Vec<&String> = seed_map.keys().collect();
    profile_names.sort();

    let mut score_params = ScoreParams {
        forward_score_nats: 0.0,
        null_score_nats: 0.0,
        bias_correction_score_nats: 0.0,
        target_count,
    };

    for profile_name in profile_names {
        let profile = profile_map
            .get_mut(profile_name)
            .ok_or_else(|| ProfileNotFoundError {
                profile_name: profile_name.clone(),
            })?;

        let seeds = seed_map
            .get(profile_name)
            .ok_or_else(|| SeedsNotFoundError {
                profile_name: profile_name.clone(),
            })?;

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
            )?;

            cloud_search_backward(
                profile,
                target,
                seed,
                &mut cloud_matrix,
                &CloudSearchParams::default(),
                &mut backward_bounds,
            )?;

            // let forward_json = forward_bounds.json("forward");
            // let backward_json = backward_bounds.json("backward");

            CloudBoundGroup::join_bounds(&mut forward_bounds, &backward_bounds)?;

            if !forward_bounds.valid() {
                println!("cloud bound fail");
                continue;
            }

            forward_bounds.trim_wings();
            // let joined_json = forward_bounds.json("joined");

            let row_bounds = RowBounds::new(&forward_bounds);

            // let row_json = row_bounds.json();

            // let mut out = PathBuf::from("bounds.json").open(true)?;
            // write!(
            //     out,
            //     "{}",
            //     json!({
            //         "start": 0,
            //         "end": profile.length + 1,
            //         "rowCount": target.length + 1,
            //         "forwardBounds": forward_json.0,
            //         "forwardDiagonal": forward_json.1,
            //         "backwardBounds": backward_json.0,
            //         "backwardDiagonal": backward_json.1,
            //         "joinedBounds": joined_json.0,
            //         "joinedDiagonal": joined_json.1,
            //         "rowBounds": row_json,
            //     })
            // )?;

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

            // let mut trace_out = PathBuf::from(format!("{}.trace", profile_name)).open(true)?;
            // trace.dump(&mut trace_out, profile, target)?;

            score_params.null_score_nats = null1_score(target.length);
            score_params.bias_correction_score_nats =
                null2_score(&posterior_matrix, profile, target, &row_bounds);

            // let mut post_out = PathBuf::from(format!("post.mtx")).open(true)?;
            // posterior_matrix.dump(&mut post_out)?;

            alignments.push(Alignment::from_trace(
                &trace,
                profile,
                target,
                &score_params,
            ));
        }
    }

    alignments = alignments
        .drain(..)
        .filter(|a| a.evalue <= args.evalue_cutoff)
        .collect::<Vec<Alignment>>();

    write_tabular_output(&alignments, &mut args.paths.results.open(true)?)?;
    write_standard_output(&alignments, &mut stdout())?;

    Ok(())
}

// fn threaded_align() {
//     #[derive(Clone)]
//     struct Step {
//         profile_name: String,
//         target_name: String,
//         target_start: usize,
//         target_end: usize,
//         profile_start: usize,
//         profile_end: usize,
//     }
//
//     let n_seeds: usize = seed_map.values().fold(0, |acc, s| acc + s.len());
//     let n_seeds_per_thread = n_seeds / args.threads + 1;
//     let mut thread_seeds: Vec<Vec<Step>> =
//         vec![Vec::with_capacity(n_seeds_per_thread); args.threads];
//
//     let mut current_vec = &mut thread_seeds[0];
//     let mut thread_idx: usize = 0;
//     let mut n_placed_in_thread: usize = 0;
//
//     for profile_name in profile_names {
//         let seeds = seed_map
//             .get(profile_name)
//             .ok_or_else(|| SeedsNotFoundError {
//                 profile_name: profile_name.clone(),
//             })?;
//
//         for seed in seeds {
//             if n_placed_in_thread >= n_seeds_per_thread {
//                 thread_idx += 1;
//                 current_vec = &mut thread_seeds[thread_idx];
//                 n_placed_in_thread = 0;
//             }
//             current_vec.push(Step {
//                 profile_name: profile_name.clone(),
//                 target_name: seed.target_name.clone(),
//                 target_start: seed.target_start,
//                 target_end: seed.target_end,
//                 profile_start: seed.profile_start,
//                 profile_end: seed.profile_end,
//             });
//             n_placed_in_thread += 1;
//         }
//     }
//
//     let time_at_spawn = Instant::now();
//     thread_seeds.into_par_iter().for_each_with(
//         (
//             time_at_spawn,
//             profile_map,
//             cloud_matrix,
//             forward_bounds,
//             backward_bounds,
//             forward_matrix,
//             backward_matrix,
//             posterior_matrix,
//             optimal_matrix,
//         ),
//         |data_structures, seeds| {
//             let (
//                 time_at_spawn,
//                 ref mut profile_map,
//                 ref mut cloud_matrix,
//                 ref mut forward_bounds,
//                 ref mut backward_bounds,
//                 ref mut forward_matrix,
//                 ref mut backward_matrix,
//                 ref mut posterior_matrix,
//                 ref mut optimal_matrix,
//             ) = data_structures;
//             println!(
//                 "thread copy time: {:3.2}",
//                 time_at_spawn.elapsed().as_secs_f32()
//             );
//             for seed in seeds {
//                 let profile = profile_map.get_mut(&seed.profile_name).unwrap();
//                 let target = target_map.get(&seed.target_name).unwrap();
//                 profile.configure_for_target_length(target.length);
//
//                 cloud_matrix.reuse(profile.length);
//                 forward_bounds.reuse(target.length, profile.length);
//                 backward_bounds.reuse(target.length, profile.length);
//
//                 let s = Seed {
//                     target_name: "".to_string(),
//                     target_start: seed.target_start,
//                     target_end: seed.target_end,
//                     profile_start: seed.profile_start,
//                     profile_end: seed.profile_end,
//                 };
//
//                 cloud_search_forward(
//                     profile,
//                     target,
//                     &s,
//                     cloud_matrix,
//                     &CloudSearchParams::default(),
//                     forward_bounds,
//                 )
//                     .unwrap();
//
//                 cloud_search_backward(
//                     profile,
//                     target,
//                     &s,
//                     cloud_matrix,
//                     &CloudSearchParams::default(),
//                     backward_bounds,
//                 )
//                     .unwrap();
//
//                 CloudBoundGroup::join_bounds(forward_bounds, backward_bounds).unwrap();
//
//                 forward_bounds.trim_wings();
//
//                 let row_bounds = RowBounds::new(&forward_bounds);
//
//                 forward_matrix.reuse(target.length, profile.length, &row_bounds);
//                 backward_matrix.reuse(target.length, profile.length, &row_bounds);
//                 posterior_matrix.reuse(target.length, profile.length, &row_bounds);
//                 optimal_matrix.reuse(target.length, profile.length, &row_bounds);
//
//                 forward_bounded(profile, target, forward_matrix, &row_bounds);
//
//                 backward_bounded(profile, target, backward_matrix, &row_bounds);
//
//                 posterior_bounded(
//                     profile,
//                     forward_matrix,
//                     backward_matrix,
//                     posterior_matrix,
//                     &row_bounds,
//                 );
//
//                 optimal_accuracy_bounded(profile, posterior_matrix, optimal_matrix, &row_bounds);
//
//                 let mut trace = Trace::new(target.length, profile.length);
//                 traceback_bounded(
//                     profile,
//                     posterior_matrix,
//                     optimal_matrix,
//                     &mut trace,
//                     row_bounds.target_end,
//                 );
//             }
//         },
//     );
//
//     Ok(())
// }
