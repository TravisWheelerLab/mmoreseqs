use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Mutex;

use crate::args::Args;
use crate::extension_traits::PathBufExt;
use crate::pipeline::seed::SeedMap;

use nale::align::bounded::structs::{
    CloudBoundGroup, CloudMatrixLinear, CloudSearchParams, DpMatrixSparse, RowBounds, Seed,
};
use nale::align::bounded::{
    backward_bounded, cloud_search_backward, cloud_search_forward, forward_bounded, null1_score,
    null2_score, optimal_accuracy_bounded, posterior_bounded, traceback_bounded,
};
use nale::structs::{Alignment, Profile, Sequence, Trace};

use nale::structs::alignment::ScoreParams;
use rayon::prelude::*;
use thread_local::ThreadLocal;

/// Each thread gets one profile and a seed list for that profile
///
/// DP structs are thread local, in RefCells
///
/// Mutex on single output file
pub fn align_threaded_c(
    args: &Args,
    mut profiles: Vec<Profile>,
    seed_map: SeedMap,
) -> anyhow::Result<()> {
    let targets = Sequence::amino_from_fasta(&args.paths.target)?;

    let results_writer: Mutex<BufWriter<File>> = Mutex::new(args.paths.results.open(true)?);

    let score_params = ScoreParams::new(targets.len());

    let mut target_map: HashMap<String, Sequence> = HashMap::new();
    for target in targets {
        target_map.insert(target.name.clone(), target);
    }

    let mut profile_seeds_pairs: Vec<(&mut Profile, &Vec<Seed>)> = vec![];

    for profile in profiles.iter_mut() {
        match seed_map.get(&profile.name) {
            Some(seeds) => profile_seeds_pairs.push((profile, seeds)),
            None => {
                continue;
            }
        }
    }

    let tl_cloud_matrix = ThreadLocal::new();
    let tl_forward_bounds = ThreadLocal::new();
    let tl_backward_bounds = ThreadLocal::new();
    let tl_forward_matrix = ThreadLocal::new();
    let tl_backward_matrix = ThreadLocal::new();
    let tl_posterior_matrix = ThreadLocal::new();
    let tl_optimal_matrix = ThreadLocal::new();

    profile_seeds_pairs.into_par_iter().for_each_with(
        score_params,
        |score_params, (profile, seeds)| {
            for seed in seeds {
                let mut cloud_matrix = tl_cloud_matrix
                    .get_or(|| RefCell::new(CloudMatrixLinear::default()))
                    .borrow_mut();

                let mut forward_bounds = tl_forward_bounds
                    .get_or(|| RefCell::new(CloudBoundGroup::default()))
                    .borrow_mut();

                let mut backward_bounds = tl_backward_bounds
                    .get_or(|| RefCell::new(CloudBoundGroup::default()))
                    .borrow_mut();

                let mut forward_matrix = tl_forward_matrix
                    .get_or(|| RefCell::new(DpMatrixSparse::default()))
                    .borrow_mut();

                let mut backward_matrix = tl_backward_matrix
                    .get_or(|| RefCell::new(DpMatrixSparse::default()))
                    .borrow_mut();

                let mut posterior_matrix = tl_posterior_matrix
                    .get_or(|| RefCell::new(DpMatrixSparse::default()))
                    .borrow_mut();

                let mut optimal_matrix = tl_optimal_matrix
                    .get_or(|| RefCell::new(DpMatrixSparse::default()))
                    .borrow_mut();

                let target = target_map.get(&seed.target_name).unwrap();
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

                score_params.null_score_nats = null1_score(target.length);
                score_params.bias_correction_score_nats =
                    null2_score(&posterior_matrix, profile, target, &row_bounds);

                optimal_accuracy_bounded(
                    profile,
                    &posterior_matrix,
                    &mut optimal_matrix,
                    &row_bounds,
                );

                let mut trace = Trace::new(target.length, profile.length);
                traceback_bounded(
                    profile,
                    &posterior_matrix,
                    &optimal_matrix,
                    &mut trace,
                    row_bounds.target_end,
                );

                let alignment = Alignment::from_trace(&trace, profile, target, score_params);

                if alignment.evalue <= args.evalue_cutoff {
                    let mut writer = results_writer.lock().unwrap();
                    writeln!(writer, "{}", alignment.tab_string());
                }
            }
        },
    );

    Ok(())
}
