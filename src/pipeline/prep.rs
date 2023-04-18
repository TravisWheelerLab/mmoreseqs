use crate::extension_traits::CommandExt;
use crate::{Args, FileFormat};
use std::process::Command;

pub fn prep(args: &Args) -> anyhow::Result<()> {
    match args.query_format {
        FileFormat::Fasta => {
            Command::new("mmseqs")
                .arg("createdb")
                .arg(&args.paths.query)
                .arg(&args.mmseqs_query_db())
                .run()?;
        }
        FileFormat::Stockholm => {
            // the msa db is only used here, so it doesn't have an associated method on Args
            let msa_db_path = &args.paths.prep_dir.join("msaDB");
            Command::new("mmseqs")
                .arg("convertmsa")
                .arg(&args.paths.query)
                .arg(&msa_db_path)
                .run()?;

            Command::new("mmseqs")
                .arg("msa2profile")
                .arg(&msa_db_path)
                .arg(&args.mmseqs_query_db())
                .args(["--threads", &args.threads.to_string()])
                // --match-mode INT       0: Columns that have a residue in the first sequence are kept,
                //                        1: columns that have a residue in --match-ratio of all sequences
                //                           are kept [0]
                .args(["--match-mode", "1"])
                .run()?;
        }
        _ => {
            panic!()
        }
    }

    // create the mmseqs2 target sequence database
    Command::new("mmseqs")
        .arg("createdb")
        .arg(&args.paths.target)
        .arg(&args.mmseqs_target_db())
        .run()?;

    // create a P7 HMM from the query
    Command::new("hmmbuild")
        .args(["--cpu", &args.threads.to_string()])
        .arg(&args.query_hmm())
        .arg(&args.paths.query)
        .run()?;

    Ok(())
}
