use crate::args::{guess_query_format_from_query_file, FileFormat};
use crate::extension_traits::CommandExt;
use std::fs::create_dir_all;

use std::path::{Path, PathBuf};
use std::process::Command;

use crate::pipeline::InvalidFileFormatError;
use anyhow::{Context, Result};
use clap::Args;
use nale::structs::Sequence;

#[derive(Args)]
pub struct PrepArgs {
    /// Query file
    #[arg(value_name = "QUERY.[fasta:sto]")]
    pub query_path: PathBuf,
    /// Target file
    #[arg(value_name = "TARGET.fasta")]
    pub target_path: PathBuf,
    /// Where to place the prepared files
    #[arg(short = 'p', long = "prep", default_value = "./prep/")]
    pub prep_dir_path: PathBuf,
    /// The number of threads to use
    #[arg(
        short = 't',
        long = "threads",
        default_value_t = 8usize,
        value_name = "n"
    )]
    pub num_threads: usize,
    /// Don't build a profile HMM with the input MSA
    #[arg(long, action)]
    pub skip_hmmbuild: bool,
}

pub trait PrepPaths {
    fn prep_dir_path(&self) -> &PathBuf;
    /// Produce a path to the query P7 HMM file
    fn prep_query_hmm_path(&self) -> PathBuf {
        self.prep_dir_path().join("query.hmm")
    }
    /// Produce a path to the MMseqs2 MSA database.
    ///
    /// This is created if a stockholm file query is provided.
    fn mmseqs_msa_db_path(&self) -> PathBuf {
        self.prep_dir_path().join("msaDB")
    }
    /// Produce a path to the MMseqs2 query database.
    ///
    /// If a fasta target was provided, this will be a sequence database.
    /// If a stockholm target was provided, this will be a profile database.
    fn mmseqs_query_db_path(&self) -> PathBuf {
        self.prep_dir_path().join("queryDB")
    }
    /// Produce a path to the MMseqs2 query database dbtype file.
    ///
    /// This file holds a byte that describes the original query file format.
    fn mmseqs_query_dbtype_path(&self) -> PathBuf {
        self.prep_dir_path().join("queryDB.dbtype")
    }
    /// Produce a path to the MMseqs2 query database index
    fn mmseqs_query_db_index_path(&self) -> PathBuf {
        self.prep_dir_path().join("queryDB.index")
    }
    /// Produce a path to the MMseqs2 query database h file
    fn mmseqs_query_db_h_path(&self) -> PathBuf {
        self.prep_dir_path().join("queryDB_h")
    }
    /// Produce a path to the MMseqs2 query database h file index
    fn mmseqs_query_db_h_index_path(&self) -> PathBuf {
        self.prep_dir_path().join("queryDB_h.index")
    }
    /// Produce a path to the MMseqs2 target database.
    ///
    /// This will always be a sequence database.
    fn mmseqs_target_db_path(&self) -> PathBuf {
        self.prep_dir_path().join("targetDB")
    }
    /// Produce a path to the MMseqs2 prefilter database.
    ///
    /// This is the result of running `mmseqs prefilter` on the query and target databases.
    fn mmseqs_prefilter_db_path(&self) -> PathBuf {
        self.prep_dir_path().join("prefilterDB")
    }
    /// Produce a path to the MMseqs2 alignment database.
    ///
    /// This is the result of running `mmseqs align` on the query, target, and prefilter databases.
    fn mmseqs_align_db_path(&self) -> PathBuf {
        self.prep_dir_path().join("alignDB")
    }
    /// Produce a path to the MMseqs2 alignment output.
    ///
    /// This is the result of running `mmseqs convertalis` on the query, target, and align databases.
    fn mmseqs_align_tsv_path(&self) -> PathBuf {
        self.prep_dir_path().join("align.tsv")
    }
}

impl PrepPaths for PrepArgs {
    fn prep_dir_path(&self) -> &PathBuf {
        &self.prep_dir_path
    }
}

pub fn prep(args: &PrepArgs) -> Result<()> {
    let query_format = guess_query_format_from_query_file(&args.query_path)?;

    create_dir_all(&args.prep_dir_path)?;

    match query_format {
        FileFormat::Fasta => {
            Command::new("mmseqs")
                .arg("createdb")
                .arg(&args.query_path)
                .arg(&args.mmseqs_query_db_path())
                .run()?;

            if !args.skip_hmmbuild {
                build_hmm_from_fasta(
                    &args.query_path,
                    &args.prep_query_hmm_path(),
                    args.num_threads,
                )?;
            }
        }
        FileFormat::Stockholm => {
            Command::new("mmseqs")
                .arg("convertmsa")
                .arg(&args.query_path)
                .arg(&args.mmseqs_msa_db_path())
                .run()?;

            Command::new("mmseqs")
                .arg("msa2profile")
                .arg(&args.mmseqs_msa_db_path())
                .arg(&args.mmseqs_query_db_path())
                .args(["--threads", &args.num_threads.to_string()])
                // --match-mode INT       0: Columns that have a residue in the first sequence are kept,
                //                        1: columns that have a residue in --match-ratio of all sequences
                //                           are kept [0]
                .args(["--match-mode", "1"])
                .run()?;

            if !args.skip_hmmbuild {
                build_hmm_from_stockholm(
                    &args.query_path,
                    &args.prep_query_hmm_path(),
                    args.num_threads,
                )?;
            }
        }
        ref format => {
            return Err(InvalidFileFormatError {
                format: format.clone(),
            })
            .context("invalid query file format in mmoreseqs prep")
        }
    }

    Command::new("mmseqs")
        .arg("createdb")
        .arg(&args.target_path)
        .arg(&args.mmseqs_target_db_path())
        .run()?;

    Ok(())
}

pub fn build_hmm_from_stockholm(
    stockholm_path: &impl AsRef<Path>,
    hmm_path: &impl AsRef<Path>,
    num_threads: usize,
) -> Result<()> {
    Command::new("hmmbuild")
        .args(["--cpu", &num_threads.to_string()])
        .arg(hmm_path.as_ref())
        .arg(stockholm_path.as_ref())
        .run()?;

    Ok(())
}

pub fn build_hmm_from_fasta(
    fasta_path: &impl AsRef<Path>,
    hmm_path: &impl AsRef<Path>,
    num_threads: usize,
) -> Result<()> {
    let fasta_path = fasta_path.as_ref();
    let hmm_path = hmm_path.as_ref();

    let query_seq = Sequence::amino_from_fasta(fasta_path).with_context(|| {
        format!(
            "failed to parse query fasta: {}",
            fasta_path.to_string_lossy()
        )
    })?;

    if query_seq.len() != 1 {
        panic!("multiple fasta queries are not supported at this time");
    }

    Command::new("hmmbuild")
        .args(["--cpu", &num_threads.to_string()])
        .args(["-n", &query_seq[0].name])
        .arg(hmm_path)
        .arg(fasta_path)
        .run()?;

    Ok(())
}
