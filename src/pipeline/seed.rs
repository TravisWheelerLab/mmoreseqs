use crate::extension_traits::{CommandExt, PathBufExt};
use crate::Args;

use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::process::Command;

use nale::align::bounded::structs::Seed;
use nale::align::needleman_wunsch::{needleman_wunsch, SimpleTraceStep};
use nale::structs::hmm::parse_hmms_from_p7hmm_file;
use nale::structs::{Profile, Sequence};

use anyhow::Context;
use thiserror::Error;

pub fn seed(args: &Args) -> anyhow::Result<(Vec<Profile>, HashMap<String, Vec<Seed>>)> {
    Command::new("mmseqs")
        .arg("prefilter")
        .arg(&args.mmseqs_query_db())
        .arg(&args.mmseqs_target_db())
        .arg(&args.mmseqs_prefilter_db())
        .args(["--threads", &args.threads.to_string()])
        // -k INT                    k-mer length (0: automatically set to optimum) [0]
        // .args(["-k", "7"])
        // --k-score INT             k-mer threshold for generating similar k-mer lists [2147483647]
        .args(["--k-score", "80"])
        // --min-ungapped-score INT  Accept only matches with ungapped alignment score above
        //                             threshold [15]
        .args(["--min-ungapped-score", "15"])
        // --max-seqs INT            Maximum results per query sequence allowed to pass the
        //                             prefilter (affects sensitivity) [300]
        .args(["--max-seqs", "1000"])
        .run()?;

    Command::new("mmseqs")
        .arg("align")
        .arg(&args.mmseqs_query_db())
        .arg(&args.mmseqs_target_db())
        .arg(&args.mmseqs_prefilter_db())
        .arg(&args.mmseqs_align_db())
        .args(["--threads", &args.threads.to_string()])
        // -e DOUBLE      List matches below this E-value (range 0.0-inf) [1.000E-03]
        .args(["-e", "1e-2"])
        // --alt-ali INT  Show up to this many alternative alignments [0]
        .args(["--alt-ali", "0"])
        .args(["-a", "1"])
        .run()?;

    Command::new("mmseqs")
        .arg("convertalis")
        .arg(&args.mmseqs_query_db())
        .arg(&args.mmseqs_target_db())
        .arg(&args.mmseqs_align_db())
        .arg(&args.mmseqs_align_tsv())
        .args(["--threads", &args.threads.to_string()])
        .args([
            "--format-output",
            "query,target,qstart,qend,tstart,tend,evalue",
        ])
        .run()?;

    let hmms = parse_hmms_from_p7hmm_file(args.query_hmm().to_str().unwrap())?;
    let p7_profiles: Vec<Profile> = hmms.iter().map(Profile::new).collect();

    let profile_to_profile_idx_maps_by_accession =
        map_p7_to_mmseqs_profiles(&p7_profiles, args).context("failed to map profiles")?;

    let profile_seeds_by_accession =
        build_alignment_seeds(&profile_to_profile_idx_maps_by_accession, args)
            .context("failed to build alignment seeds")?;

    let mut seeds_out = args
        .paths
        .seeds
        .open(true)
        .context("failed to create alignment seeds file")?;

    write!(
        seeds_out,
        "{}",
        serde_json::to_string(&profile_seeds_by_accession)?
    )
    .context("failed to write alignment seeds")?;

    Ok((p7_profiles, profile_seeds_by_accession))
}

pub fn map_p7_to_mmseqs_profiles(
    p7_profiles: &[Profile],
    args: &Args,
) -> anyhow::Result<HashMap<String, Vec<usize>>> {
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

pub fn extract_mmseqs_profile_consensus_sequences(
    args: &Args,
) -> anyhow::Result<HashMap<String, Sequence>> {
    let mut offsets_and_lengths: Vec<(usize, usize)> = vec![];
    let mut accession_numbers: Vec<String> = vec![];

    let query_db_h_index_file =
        File::open(&args.mmseqs_query_db_h_index()).context("failed to open queryDB_h.index")?;

    let reader = BufReader::new(query_db_h_index_file);
    for line in reader.lines() {
        match line {
            Ok(line) => {
                let tokens: Vec<&str> = line.split_whitespace().collect();

                let offset = tokens[1].parse::<usize>()?;
                let length = tokens[2].parse::<usize>()?;
                offsets_and_lengths.push((offset, length));
            }
            Err(e) => {
                return Err(e).context("failed to parse line in queryDB_h.index");
            }
        }
    }

    let mut query_db_h_file =
        File::open(&args.mmseqs_query_db_h()).context("failed to open queryDB_h")?;

    for (offset, length) in &offsets_and_lengths {
        let mut buffer = vec![0; *length];
        query_db_h_file.seek(SeekFrom::Start(*offset as u64))?;
        query_db_h_file.read_exact(&mut buffer)?;

        let mut accession_string: Option<String> = None;
        for (buf_idx, byte) in buffer.iter().enumerate() {
            if byte.is_ascii_whitespace() {
                accession_string = Some(
                    std::str::from_utf8(&buffer[0..buf_idx])
                        .context("failed to create accession string")?
                        .to_string(),
                );
                break;
            }
        }

        match accession_string {
            Some(accession) => accession_numbers.push(accession),
            None => {
                panic!()
            }
        }
    }

    let query_db_index_file =
        File::open(&args.mmseqs_query_db_index()).context("failed to open queryDB.index")?;

    let reader = BufReader::new(query_db_index_file);
    for line in reader.lines() {
        match line {
            Ok(line) => {
                let tokens: Vec<&str> = line.split_whitespace().collect();

                let line_idx = tokens[0].parse::<usize>()?;
                let offset = tokens[1].parse::<usize>()?;
                let length = tokens[2].parse::<usize>()?;
                offsets_and_lengths[line_idx] = (offset, length);
            }
            Err(e) => {
                return Err(e).context("failed to parse line in queryDB.index");
            }
        }
    }

    let mut sequence_map: HashMap<String, Sequence> = HashMap::new();

    let mut query_db_file =
        File::open(&args.mmseqs_query_db()).context("failed to open queryDB")?;

    for (seq_idx, (offset, length)) in offsets_and_lengths.iter().enumerate() {
        let mut buffer = vec![0; *length];
        query_db_file.seek(SeekFrom::Start(*offset as u64))?;
        query_db_file.read_exact(&mut buffer)?;

        let mut consensus_digital_bytes: Vec<u8> = vec![];

        for byte_chunk in buffer.chunks(23) {
            if byte_chunk.len() == 23 {
                consensus_digital_bytes.push(byte_chunk[21]);
            }
        }

        sequence_map.insert(
            accession_numbers[seq_idx].clone(),
            Sequence::from_digital(&consensus_digital_bytes)?,
        );
    }

    Ok(sequence_map)
}

#[derive(Error, Debug)]
#[error("no profile to profile map")]
pub struct ProfilesNotMappedError;

pub fn build_alignment_seeds(
    profile_to_profile_idx_maps_by_accession: &HashMap<String, Vec<usize>>,
    args: &Args,
) -> anyhow::Result<HashMap<String, Vec<Seed>>> {
    let mut profile_seeds_by_accession: HashMap<String, Vec<Seed>> = HashMap::new();

    let mmseqs_align_file = File::open(&args.mmseqs_align_tsv()).context(format!(
        "couldn't open mmseqs align file at: {}",
        &args.mmseqs_align_tsv().to_string_lossy()
    ))?;
    let align_reader = BufReader::new(mmseqs_align_file);

    for line in align_reader.lines().flatten() {
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
