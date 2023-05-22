use crate::args::{Args, FileFormat};
use crate::extension_traits::CommandExt;
use anyhow::{Context, Result};
use nale::structs::Sequence;
use std::process::Command;

pub fn prep(args: &Args) -> anyhow::Result<()> {
    match args.query_format {
        FileFormat::Fasta => {
            Command::new("mmseqs")
                .arg("createdb")
                .arg(&args.paths.query)
                .arg(args.mmseqs_query_db())
                .run()?;

            if args.build_hmm {
                build_hmm_from_fasta(args)?;
            }
        }
        FileFormat::Stockholm => {
            // the msa db is only used here, so it doesn't have an associated method on Args
            let msa_db_path = &args.paths.prep_dir.join("msaDB");
            Command::new("mmseqs")
                .arg("convertmsa")
                .arg(&args.paths.query)
                .arg(msa_db_path)
                .run()?;

            Command::new("mmseqs")
                .arg("msa2profile")
                .arg(msa_db_path)
                .arg(args.mmseqs_query_db())
                .args(["--threads", &args.threads.to_string()])
                // --match-mode INT       0: Columns that have a residue in the first sequence are kept,
                //                        1: columns that have a residue in --match-ratio of all sequences
                //                           are kept [0]
                .args(["--match-mode", "1"])
                .run()?;

            build_hmm_from_stockholm(args)?;
        }
        _ => {
            panic!()
        }
    }

    Command::new("mmseqs")
        .arg("createdb")
        .arg(&args.paths.target)
        .arg(&args.mmseqs_target_db())
        .run()?;

    Ok(())
}

pub fn build_hmm_from_stockholm(args: &Args) -> Result<()> {
    Command::new("hmmbuild")
        .args(["--cpu", &args.threads.to_string()])
        .arg(&args.query_hmm())
        .arg(&args.paths.query)
        .run()?;

    Ok(())
}

pub fn build_hmm_from_fasta(args: &Args) -> Result<()> {
    let query_seq = Sequence::amino_from_fasta(&args.paths.query).with_context(|| {
        format!(
            "failed to parse query fasta: {}",
            &args.paths.query.to_string_lossy()
        )
    })?;

    if query_seq.len() != 1 {
        panic!("multiple fasta queries are not supported at this time");
    }

    Command::new("hmmbuild")
        .args(["--cpu", &args.threads.to_string()])
        .args(["-n", &query_seq[0].name])
        .arg(args.query_hmm())
        .arg(&args.paths.query)
        .run()?;

    Ok(())
}
