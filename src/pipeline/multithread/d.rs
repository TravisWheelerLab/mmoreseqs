use crate::args::Args;
use crate::extension_traits::PathBufExt;
use crate::pipeline::seed::SeedMap;
use std::cell::RefCell;

use nale::align::bounded::structs::{
    CloudBoundGroup, CloudMatrixLinear, CloudSearchParams, DpMatrixSparse, RowBounds,
};
use nale::align::bounded::{
    backward_bounded, cloud_search_backward, cloud_search_forward, forward_bounded, null1_score,
    null2_score, optimal_accuracy_bounded, posterior_bounded, traceback_bounded,
};
use nale::structs::alignment::ScoreParams;
use nale::structs::{Alignment, Profile, Sequence, Trace};

use std::collections::HashMap;
use std::io::Write;
use std::sync::Mutex;

use rayon::iter::{IntoParallelIterator, ParallelIterator};
use thread_local::ThreadLocal;

/// Each thread gets a copy of all the profiles and a seed map
///
/// DP structs initialized once in a thread
///
/// Mutex on single output file
pub fn align_threaded_d(
    args: &Args,
    profiles: Vec<Profile>,
    mut seed_map: SeedMap,
) -> anyhow::Result<()> {
    let targets = Sequence::amino_from_fasta(&args.paths.target)?;

    let score_params = ScoreParams::new(targets.len());

    let mut target_map: HashMap<String, Sequence> = HashMap::new();
    for target in targets {
        target_map.insert(target.name.clone(), target);
    }

    // TODO: REMOVE HARDCODED THREAD COUNT HERE
    let mut thread_seed_maps: Vec<SeedMap> = vec![HashMap::new(); 8];
    let mut thread_idx: usize = 0;

    for profile in &profiles {
        let seeds = match seed_map.remove(&profile.name) {
            Some(seeds) => seeds,
            None => {
                continue;
            }
        };

        for seed in seeds {
            match thread_seed_maps[thread_idx].get_mut(&profile.name) {
                None => {
                    thread_seed_maps[thread_idx].insert(profile.name.clone(), vec![seed]);
                }
                Some(vec) => {
                    vec.push(seed);
                }
            }
            thread_idx += 1;
            // TODO: REMOVE HARDCODED THREAD COUNT HERE
            if thread_idx >= 8 {
                thread_idx = 0;
            }
        }
    }

    let thread_writer = ThreadLocal::new();
    let thread_count = Mutex::new(0);

    thread_seed_maps.into_par_iter().for_each_with(
        (profiles, &target_map, score_params),
        |(profiles, target_map, score_params), seed_map| {
            let mut cloud_matrix = CloudMatrixLinear::default();
            let mut forward_bounds = CloudBoundGroup::default();
            let mut backward_bounds = CloudBoundGroup::default();
            let mut forward_matrix = DpMatrixSparse::default();
            let mut backward_matrix = DpMatrixSparse::default();
            let mut posterior_matrix = DpMatrixSparse::default();
            let mut optimal_matrix = DpMatrixSparse::default();

            for profile in profiles.iter_mut() {
                let seeds = match seed_map.get(&profile.name) {
                    Some(seeds) => seeds,
                    None => {
                        continue;
                    }
                };

                for seed in seeds {
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
                        let mut writer = thread_writer
                            .get_or(|| {
                                let mut cnt = thread_count.lock().unwrap();
                                *cnt += 1;

                                RefCell::new(
                                    args.paths
                                        .results
                                        .with_extension(format!("{cnt}"))
                                        .open(true)
                                        .unwrap(),
                                )
                            })
                            .borrow_mut();

                        writeln!(writer, "{}", alignment.tab_string());
                    }
                }
            }
        },
    );

    Ok(())
}
