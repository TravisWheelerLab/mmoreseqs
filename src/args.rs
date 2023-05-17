use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::PathBuf;

use anyhow::Context;
use thiserror::Error;

#[derive(Default)]
pub struct Paths {
    /// The query file provided at the command line
    pub query: PathBuf,
    /// The target file provided at the command line
    pub target: PathBuf,
    /// The directory under which prepared/intermediate files are placed
    pub prep_dir: PathBuf,
    /// The alignment seeds used for forward/backward
    pub seeds: PathBuf,
    /// The path that results will be written to
    pub results: PathBuf,
}

#[derive(Default)]
pub enum Command {
    Prep,
    Seed,
    Align,
    Search,
    #[default]
    NotSet,
}

#[derive(Default)]
/// The arguments that are passed throughout the pipeline
pub struct Args {
    pub command: Command,
    pub paths: Paths,
    pub query_format: FileFormat,
    pub threads: usize,
    pub evalue_cutoff: f64,
}

impl Args {
    pub fn guess_query_format(&mut self) -> anyhow::Result<()> {
        let file = File::open(&self.paths.query).context(format!(
            "failed to open query file: {}",
            &self.paths.query.to_string_lossy()
        ))?;

        let mut reader = BufReader::new(file);
        let mut first_line = String::new();
        reader.read_line(&mut first_line)?;

        if &first_line[0..1] == ">" {
            self.query_format = FileFormat::Fasta;
        } else if &first_line[0..11] == "# STOCKHOLM" {
            self.query_format = FileFormat::Stockholm;
        } else if &first_line[0..5] == "HMMER" {
            self.query_format = FileFormat::Hmm;
        } else {
            return Err(UnrecognizedFileFormatError).context(format!(
                "couldn't guess the format of query file: {}",
                &self.paths.query.to_string_lossy()
            ));
        };
        Ok(())
    }

    pub fn get_query_format_from_mmseqs_file(&mut self) -> anyhow::Result<()> {
        let mut file = File::open(self.paths.prep_dir.join("queryDB.dbtype"))?;
        let mut dbtype_buf: Vec<u8> = vec![];
        file.read_to_end(&mut dbtype_buf)?;
        //  from mmseqs2: commons/parameters.h
        //      DBTYPE_AMINO_ACIDS = 0;
        //      DBTYPE_NUCLEOTIDES = 1;
        //      DBTYPE_HMM_PROFILE = 2;
        //      //DBTYPE_PROFILE_STATE_SEQ = 3;
        //      //DBTYPE_PROFILE_STATE_PROFILE = 4;
        //      DBTYPE_ALIGNMENT_RES = 5;
        //      DBTYPE_CLUSTER_RES = 6;
        //      DBTYPE_PREFILTER_RES = 7;
        //      DBTYPE_TAXONOMICAL_RESULT = 8;
        //      DBTYPE_INDEX_DB = 9;
        //      DBTYPE_CA3M_DB = 10;
        //      DBTYPE_MSA_DB = 11;
        //      DBTYPE_GENERIC_DB = 12;
        //      DBTYPE_OMIT_FILE = 13;
        //      DBTYPE_PREFILTER_REV_RES = 14;
        //      DBTYPE_OFFSETDB = 15;
        //      DBTYPE_DIRECTORY = 16; // needed for verification
        //      DBTYPE_FLATFILE = 17; // needed for verification
        //      DBTYPE_SEQTAXDB = 18; // needed for verification
        //      DBTYPE_STDIN = 19; // needed for verification
        //      DBTYPE_URI = 20; // needed for verification
        self.query_format = match dbtype_buf[0] {
            0u8 => FileFormat::Fasta,
            2u8 => FileFormat::Stockholm,
            _ => {
                return Err(UnsupportedMmseqsDbError {
                    code: dbtype_buf[0],
                }
                .into())
            }
        };
        Ok(())
    }

    /// Produce a path to the query P7 HMM file
    pub fn query_hmm(&self) -> PathBuf {
        self.paths.prep_dir.join("query.hmm")
    }
    /// Produce a path to the MMseqs2 query database.
    ///
    /// If a fasta target was provided, this will be a sequence database.
    /// If a stockholm target was provided, this will be a profile database.
    pub fn mmseqs_query_db(&self) -> PathBuf {
        self.paths.prep_dir.join("queryDB")
    }
    /// Produce a path to the MMseqs2 query database index
    pub fn mmseqs_query_db_index(&self) -> PathBuf {
        self.paths.prep_dir.join("queryDB.index")
    }
    /// Produce a path to the MMseqs2 query database h file
    pub fn mmseqs_query_db_h(&self) -> PathBuf {
        self.paths.prep_dir.join("queryDB_h")
    }
    /// Produce a path to the MMseqs2 query database h file index
    pub fn mmseqs_query_db_h_index(&self) -> PathBuf {
        self.paths.prep_dir.join("queryDB_h.index")
    }
    /// Produce a path to the MMseqs2 target database.
    ///
    /// This will always be a sequence database.
    pub fn mmseqs_target_db(&self) -> PathBuf {
        self.paths.prep_dir.join("targetDB")
    }
    /// Produce a path to the MMseqs2 prefilter database.
    ///
    /// This is the result of running `mmseqs prefilter` on the query and target databases.
    pub fn mmseqs_prefilter_db(&self) -> PathBuf {
        self.paths.prep_dir.join("prefilterDB")
    }
    /// Produce a path to the MMseqs2 alignment database.
    ///
    /// This is the result of running `mmseqs align` on the query, target, and prefilter databases.
    pub fn mmseqs_align_db(&self) -> PathBuf {
        self.paths.prep_dir.join("alignDB")
    }
    /// Produce a path to the MMseqs2 alignment output.
    ///
    /// This is the result of running `mmseqs convertalis` on the query, target, and align databases.
    pub fn mmseqs_align_tsv(&self) -> PathBuf {
        self.paths.prep_dir.join("align.tsv")
    }
}

#[derive(Default)]
pub enum FileFormat {
    Fasta,
    Stockholm,
    Hmm,
    #[default]
    Unset,
}

#[derive(Error, Debug)]
#[error("can't guess file format")]
pub struct UnrecognizedFileFormatError;

#[derive(Error, Debug)]
#[error("unsupported MMseqs2 database type: {code}")]
pub struct UnsupportedMmseqsDbError {
    code: u8,
}
