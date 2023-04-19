use crate::extension_traits::PathBufExt;
use crate::pipeline::seed::SeedMap;
use crate::Args;

use anyhow::Context;
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;

use nale::align::bounded::structs::{
    CloudBoundGroup, CloudMatrixLinear, CloudSearchParams, DpMatrixSparse, RowBounds,
};
use nale::align::bounded::{
    backward_bounded, cloud_search_backward, cloud_search_forward, forward_bounded,
    optimal_accuracy_bounded, posterior_bounded, traceback_bounded,
};
use nale::output::output_tabular::write_tabular_output;
use nale::structs::hmm::parse_hmms_from_p7hmm_file;
use nale::structs::{Alignment, Profile, Sequence, Trace};

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
            let hmms = parse_hmms_from_p7hmm_file(args.paths.query.to_str().unwrap())?;
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

    for profile_accession in profile_names {
        let profile = profile_map.get_mut(profile_accession).unwrap();
        let seeds = seed_map.get(profile_accession).unwrap();
        for seed in seeds {
            let target = target_map.get(&seed.target_name[..]).unwrap();
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

            CloudBoundGroup::join_bounds(&mut forward_bounds, &backward_bounds)?;

            forward_bounds.trim_wings();

            let row_bounds = RowBounds::new(&forward_bounds);

            forward_matrix.reuse(target.length, profile.length, &row_bounds);
            backward_matrix.reuse(target.length, profile.length, &row_bounds);
            posterior_matrix.reuse(target.length, profile.length, &row_bounds);
            optimal_matrix.reuse(target.length, profile.length, &row_bounds);

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

            alignments.push(Alignment::new(&trace, profile, target, target_count));
        }
    }

    alignments = alignments
        .drain(..)
        .filter(|a| a.evalue <= args.evalue_cutoff)
        .collect();

    write_tabular_output(&alignments, &mut args.paths.results.open(true)?)?;

    Ok(())
}
