# MMoreseqs

## About

MMoreseqs is a biological sequence alignment tool. Currently, only protein sequence alignment is supported.

MMoreseqs uses [MMseqs2](https://github.com/soedinglab/MMseqs2) to find rough alignment seeds to use as a starting point
for a highly sensitive, bounded sequence alignment algorithm.

This project is the continuation of pilot work done by David Rich. The original codebase can be
found [here](https://github.com/TravisWheelerLab/mmoreseqs-legacy).

## Installation

To build mmoreseqs from source, you'll first need to install Rust and Cargo.
The easiest way to do that is to use [rustup](https://rustup.rs/).

Once that's done, you can then build mmoreseqs:

    git clone https://github.com/TravisWheelerLab/mmoreseqs
    cd mmoreseqs/
    cargo build --release

You'll then find the compiled binary at: `target/release/mmoreseqs`

For example, try running:

    target/release/mmoreseqs -h

## Usage

The input to MMoreseqs is a query multiple sequence alignment (stockholm) file and a target sequence (fasta) file.
To run the MMoreseqs pipeline, use the `mmoreseqs search` command:

For example:

    $ mmoreseqs search query.sto target.fa

## License

MMoreseqs is licensed under the BSD-3-Clause license.

See `LICENSE` for details.

## Authors

Jack Roddy - jroddy@arizona.edu
