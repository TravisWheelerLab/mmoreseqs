use crate::command_ext::CommandExt;
use crate::Args;
use anyhow::{Context, Result};
use nale::structs::Sequence;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::process::Command;
use thiserror::Error;

#[derive(Error, Debug)]
#[error("no profile to profile map")]
pub struct ProfilesNotMappedError;

pub fn run_hmmbuild(args: &Args) -> Result<()> {
    Command::new("hmmbuild")
        .args(["--cpu", &args.threads.to_string()])
        .arg(&args.paths.query_hmm)
        .arg(&args.paths.query_msa)
        .run()
}

pub fn run_mmseqs_convertmsa(args: &Args) -> Result<()> {
    Command::new("mmseqs")
        .arg("convertmsa")
        .arg(&args.paths.query_msa)
        .arg(&args.paths.query_msa_db)
        .run()
}

pub fn run_mmseqs_msa2profile(args: &Args) -> Result<()> {
    Command::new("mmseqs")
        .arg("msa2profile")
        .arg(&args.paths.query_msa_db)
        .arg(&args.paths.query_db)
        .args(["--threads", &args.threads.to_string()])
        // --match-mode INT       0: Columns that have a residue in the first sequence are kept,
        //                        1: columns that have a residue in --match-ratio of all sequences
        //                           are kept [0]
        .args(["--match-mode", "1"])
        .run()
}

pub fn run_mmseqs_createdb(args: &Args) -> Result<()> {
    Command::new("mmseqs")
        .arg("createdb")
        .arg(&args.paths.target_fasta)
        .arg(&args.paths.target_db)
        .run()
}

pub fn run_mmseqs_prefilter(args: &Args) -> Result<()> {
    Command::new("mmseqs")
        .arg("prefilter")
        .arg(&args.paths.query_db)
        .arg(&args.paths.target_db)
        .arg(&args.paths.prefilter_db)
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
        .run()
}

pub fn run_mmseqs_align(args: &Args) -> Result<()> {
    Command::new("mmseqs")
        .arg("align")
        .arg(&args.paths.query_db)
        .arg(&args.paths.target_db)
        .arg(&args.paths.prefilter_db)
        .arg(&args.paths.align_db)
        .args(["--threads", &args.threads.to_string()])
        // -e DOUBLE      List matches below this E-value (range 0.0-inf) [1.000E-03]
        .args(["-e", "1e-2"])
        // --alt-ali INT  Show up to this many alternative alignments [0]
        .args(["--alt-ali", "0"])
        .args(["-a", "1"])
        .run()
}

pub fn run_mmseqs_convertalis(args: &Args) -> Result<()> {
    Command::new("mmseqs")
        .arg("convertalis")
        .arg(&args.paths.query_db)
        .arg(&args.paths.target_db)
        .arg(&args.paths.align_db)
        .arg(&args.paths.seeds)
        .args(["--threads", &args.threads.to_string()])
        .args([
            "--format-output",
            "query,target,qstart,qend,tstart,tend,evalue",
        ])
        .run()
}

pub fn extract_mmseqs_profile_consensus_sequences(
    args: &Args,
) -> Result<HashMap<String, Sequence>> {
    let mut offsets_and_lengths: Vec<(usize, usize)> = vec![];
    let mut accession_numbers: Vec<String> = vec![];

    let query_db_h_index_file =
        File::open(&args.paths.query_db_h_index).context("failed to open queryDB_h.index")?;

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
        File::open(&args.paths.query_db_h).context("failed to open queryDB_h")?;

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
        File::open(&args.paths.query_db_index).context("failed to open queryDB.index")?;

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

    let mut query_db_file = File::open(&args.paths.query_db).context("failed to open queryDB")?;

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
