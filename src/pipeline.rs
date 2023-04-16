use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::external_steps::{
    extract_mmseqs_profile_consensus_sequences, run_hmmbuild, run_mmseqs_align,
    run_mmseqs_convertalis, run_mmseqs_convertmsa, run_mmseqs_createdb, run_mmseqs_msa2profile,
    run_mmseqs_prefilter, ProfilesNotMappedError,
};
use crate::Args;

use nale::align::bounded::structs::{
    CloudBoundGroup, CloudMatrixLinear, CloudSearchParams, RowBounds, Seed,
};
use nale::align::bounded::{
    backward_bounded, cloud_search_backward, cloud_search_forward, forward_bounded,
    optimal_accuracy_bounded, posterior_bounded, traceback_bounded,
};
use nale::align::needleman_wunsch::{needleman_wunsch, SimpleTraceStep};
use nale::output::output_tabular::write_tabular_output;
use nale::output::path_buf_ext::PathBufExt;
use nale::structs::hmm::parse_hmms_from_p7hmm_file;
use nale::structs::{Alignment, DpMatrixFlat, Profile, Sequence, Trace};

use anyhow::Result;

fn map_p7_to_mmseqs_profiles(
    p7_profiles: &[Profile],
    args: &Args,
) -> Result<HashMap<String, Vec<usize>>> {
    let mmseqs_consensus_map = extract_mmseqs_profile_consensus_sequences(args)?;

    let mut profile_to_profile_idx_maps_by_accession: HashMap<String, Vec<usize>> = HashMap::new();

    for p7_profile in p7_profiles {
        let accession = &p7_profile.accession;
        let mmseqs_consensus = mmseqs_consensus_map.get(accession).unwrap();
        let p7_consensus = Sequence::from_utf8(&p7_profile.consensus_sequence[1..])?;
        let trace = needleman_wunsch(mmseqs_consensus, &p7_consensus);

        let mut mmseqs_to_p7: Vec<usize> = vec![0; mmseqs_consensus.length + 1];

        let mut mmseqs_idx: usize = 0;
        let mut p7_idx: usize = 0;
        for step in &trace {
            match step {
                SimpleTraceStep::Diagonal => {
                    mmseqs_idx += 1;
                    p7_idx += 1;
                }
                SimpleTraceStep::Up => {
                    mmseqs_idx += 1;
                }
                SimpleTraceStep::Left => {
                    p7_idx += 1;
                }
            }
            mmseqs_to_p7[mmseqs_idx] = p7_idx;
        }

        // this debug assert should guarantee that the NW
        // alignment fully covered both consensus sequences
        debug_assert_eq!(mmseqs_idx, mmseqs_consensus.length);
        debug_assert_eq!(p7_idx, p7_consensus.length);

        profile_to_profile_idx_maps_by_accession.insert(accession.clone(), mmseqs_to_p7);
    }

    Ok(profile_to_profile_idx_maps_by_accession)
}

pub fn build_alignment_seeds(
    profile_to_profile_idx_maps_by_accession: &HashMap<String, Vec<usize>>,
    args: &Args,
) -> Result<HashMap<String, Vec<Seed>>> {
    let mut profile_seeds_by_accession: HashMap<String, Vec<Seed>> = HashMap::new();

    let seeds_file = File::open(&args.paths.seeds)?;
    let seeds_buf_reader = BufReader::new(seeds_file);

    for line in seeds_buf_reader.lines().flatten() {
        let line_tokens: Vec<&str> = line.split_whitespace().collect();
        let accession = line_tokens[0];

        let seeds = match profile_seeds_by_accession.get_mut(accession) {
            Some(seeds) => seeds,
            None => {
                profile_seeds_by_accession.insert(accession.to_string(), vec![]);
                profile_seeds_by_accession.get_mut(accession).unwrap()
            }
        };

        let profile_idx_map = profile_to_profile_idx_maps_by_accession
            .get(accession)
            .ok_or(ProfilesNotMappedError)?;

        let target_name = line_tokens[1].to_string();
        let target_start = line_tokens[4].parse::<usize>()?;
        let target_end = line_tokens[5].parse::<usize>()?;
        let profile_start = line_tokens[2].parse::<usize>()?;
        let profile_end = line_tokens[3].parse::<usize>()?;

        seeds.push(Seed {
            target_name,
            target_start,
            target_end,
            profile_start: profile_idx_map[profile_start].max(1),
            profile_end: profile_idx_map[profile_end],
        })
    }
    Ok(profile_seeds_by_accession)
}

pub fn prep(args: &Args) -> Result<()> {
    run_mmseqs_convertmsa(args)?;
    run_mmseqs_msa2profile(args)?;
    run_mmseqs_createdb(args)?;
    run_hmmbuild(args)?;
    Ok(())
}

pub fn seed(args: &Args) -> Result<()> {
    run_mmseqs_prefilter(args)?;
    run_mmseqs_align(args)?;
    run_mmseqs_convertalis(args)?;
    Ok(())
}

pub fn align(args: &Args) -> Result<()> {
    let hmms = parse_hmms_from_p7hmm_file(args.paths.query_hmm.to_str().unwrap())?;
    let p7_profiles: Vec<Profile> = hmms.iter().map(Profile::new).collect();

    let profile_to_profile_idx_maps_by_accession = map_p7_to_mmseqs_profiles(&p7_profiles, args)?;

    let profile_seeds_by_accession =
        build_alignment_seeds(&profile_to_profile_idx_maps_by_accession, args)?;

    let mut profile_map: HashMap<String, Profile> = HashMap::new();
    for profile in p7_profiles {
        profile_map.insert(profile.accession.clone(), profile);
    }

    let targets = Sequence::amino_from_fasta(&args.paths.target_fasta)?;
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

    let mut forward_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);
    let mut backward_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);
    let mut posterior_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);
    let mut optimal_matrix = DpMatrixFlat::new(max_target_length, max_profile_length);

    let mut alignments: Vec<Alignment> = vec![];

    let mut profile_names: Vec<&String> = profile_seeds_by_accession.keys().collect();
    profile_names.sort();

    for profile_accession in profile_names {
        let profile = profile_map.get_mut(profile_accession).unwrap();
        let seeds = profile_seeds_by_accession.get(profile_accession).unwrap();
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

            forward_matrix.reuse(target.length, profile.length);
            backward_matrix.reuse(target.length, profile.length);
            posterior_matrix.reuse(target.length, profile.length);
            optimal_matrix.reuse(target.length, profile.length);

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

pub fn search(args: &Args) -> Result<()> {
    {
        // quickly make sure we can write the results
        args.paths.results.open(true)?;
    }
    prep(args)?;
    seed(args)?;
    align(args)?;
    Ok(())
}
