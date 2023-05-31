use std::cell::RefCell;
use std::collections::HashMap;
use std::io::Write;
use std::sync::Mutex;

use crate::extension_traits::PathBufExt;
use crate::pipeline::seed::SeedMap;

use nale::align::bounded::structs::{
    CloudBoundGroup, CloudMatrixLinear, CloudSearchParams, DpMatrixSparse, RowBounds, Seed,
};
use nale::align::bounded::{
    backward_bounded, cloud_search_backward, cloud_search_forward, forward_bounded, null1_score,
    null2_score, optimal_accuracy_bounded, posterior_bounded, traceback_bounded,
};
use nale::structs::alignment::ScoreParams;
use nale::structs::{Alignment, Profile, Sequence, Trace};

use crate::pipeline::align::AlignArgs;
use rayon::prelude::*;
use thread_local::ThreadLocal;

/// Each thread gets one profile and a seed list for that profile
///
/// DP structs copied for each thread
///
/// Mutex on single output file
pub fn align_threaded_e(
    args: &AlignArgs,
    mut profiles: Vec<Profile>,
    targets: Vec<Sequence>,
    mut seed_map: SeedMap,
) -> anyhow::Result<()> {
    let mut dp = AlignmentStructs::default();

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

    #[derive(Default, Clone)]
    struct AlignmentStructs {
        cloud_matrix: CloudMatrixLinear,
        forward_bounds: CloudBoundGroup,
        backward_bounds: CloudBoundGroup,
        forward_matrix: DpMatrixSparse,
        backward_matrix: DpMatrixSparse,
        posterior_matrix: DpMatrixSparse,
        optimal_matrix: DpMatrixSparse,
    }

    let thread_writer = ThreadLocal::new();
    let thread_count = Mutex::new(0);

    profile_seeds_pairs.into_par_iter().for_each_with(
        (dp, score_params),
        |(dp, score_params), (profile, seeds)| {
            for seed in seeds {
                let target = target_map.get(&seed.target_name).unwrap();
                profile.configure_for_target_length(target.length);

                dp.cloud_matrix.reuse(profile.length);
                dp.forward_bounds.reuse(target.length, profile.length);
                dp.backward_bounds.reuse(target.length, profile.length);

                cloud_search_forward(
                    profile,
                    target,
                    seed,
                    &mut dp.cloud_matrix,
                    &CloudSearchParams::default(),
                    &mut dp.forward_bounds,
                );

                cloud_search_backward(
                    profile,
                    target,
                    seed,
                    &mut dp.cloud_matrix,
                    &CloudSearchParams::default(),
                    &mut dp.backward_bounds,
                );

                CloudBoundGroup::join_bounds(&mut dp.forward_bounds, &dp.backward_bounds);

                if !dp.forward_bounds.valid() {
                    println!("cloud bound fail");
                    continue;
                }

                dp.forward_bounds.trim_wings();

                let row_bounds = RowBounds::new(&dp.forward_bounds);

                if !row_bounds.valid() {
                    println!("row bound fail");
                    continue;
                }

                dp.forward_matrix
                    .reuse(target.length, profile.length, &row_bounds);
                dp.backward_matrix
                    .reuse(target.length, profile.length, &row_bounds);
                dp.posterior_matrix
                    .reuse(target.length, profile.length, &row_bounds);
                dp.optimal_matrix
                    .reuse(target.length, profile.length, &row_bounds);

                // we use the forward score to compute the final bit score (later)
                score_params.forward_score_nats =
                    forward_bounded(profile, target, &mut dp.forward_matrix, &row_bounds);

                backward_bounded(profile, target, &mut dp.backward_matrix, &row_bounds);

                posterior_bounded(
                    profile,
                    &dp.forward_matrix,
                    &dp.backward_matrix,
                    &mut dp.posterior_matrix,
                    &row_bounds,
                );

                score_params.null_score_nats = null1_score(target.length);
                score_params.bias_correction_score_nats =
                    null2_score(&dp.posterior_matrix, profile, target, &row_bounds);

                optimal_accuracy_bounded(
                    profile,
                    &dp.posterior_matrix,
                    &mut dp.optimal_matrix,
                    &row_bounds,
                );

                let mut trace = Trace::new(target.length, profile.length);
                traceback_bounded(
                    profile,
                    &dp.posterior_matrix,
                    &dp.optimal_matrix,
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
                                args.tsv_results_path
                                    .with_extension(format!("{cnt}"))
                                    .open(true)
                                    .unwrap(),
                            )
                        })
                        .borrow_mut();

                    writeln!(writer, "{}", alignment.tab_string());
                }
            }
        },
    );

    Ok(())
}
