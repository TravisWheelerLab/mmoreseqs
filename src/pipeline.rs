use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::external_steps::{
    extract_mmseqs_profile_consensus_sequences, run_hmmbuild, run_mmseqs_align,
    run_mmseqs_convertmsa, run_mmseqs_createdb, run_mmseqs_msa2profile, run_mmseqs_prefilter,
    ProfilesNotMappedError,
};

use nale::align::needleman_wunsch::{needleman_wunsch, SimpleTraceStep};
use nale::pipelines::Seed;
use nale::structs::hmm::parse_hmms_from_p7hmm_file;
use nale::structs::{Profile, Sequence};

use crate::Args;
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
        profile_map.insert(profile.name.clone(), profile);
    }

    // let targets = Sequence::amino_from_fasta(&args.target)?;
    // let mut target_map: HashMap<String, Sequence> = HashMap::new();
    // for target in targets {
    //     target_map.insert(target.name.clone(), target);
    // }
    Ok(())
}

pub fn search(args: &Args) -> Result<()> {
    prep(args)?;
    seed(args)?;
    // align(args)?;
    Ok(())
}
