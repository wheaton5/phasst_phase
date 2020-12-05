#[macro_use]
extern crate clap;
extern crate fnv;
extern crate hashbrown;
extern crate flate2;
extern crate itertools;
extern crate rayon;
extern crate byteorder;
extern crate disjoint_set;
use rayon::prelude::*;
use flate2::read::GzDecoder;
use std::io::BufReader;
use std::io::BufRead;
use std::io::Read;
use itertools::izip;
use byteorder::{WriteBytesExt, LittleEndian};

use std::fs::File;
use std::io::{BufWriter, Write};

use hashbrown::{HashMap, HashSet};
use std::str;

use clap::{App};

fn main() {
    println!("Welcome to Phuzzy Phazer!");
    let params = load_params();
}


fn get_reader(filename: String) -> BufReader<Box<dyn Read>> {
    let filetype: Vec<&str> = filename.split(".").collect();
    let filetype = filetype[filetype.len()-1];
    let file = match File::open(filename.clone()) {
        Ok(file) => file,
        Err(error) => panic!("There was a problem opening the file: {:?}", error),
    };
    let reader: Box<dyn Read> = match filetype {
        "gz" => Box::new(GzDecoder::new(file)),
        _ => Box::new(file),
    };
    BufReader::new(reader)
}

#[derive(Clone)]
struct Params {
    het_kmers: String,
    txg_mols: Vec<String>,
    hic_mols: Vec<String>,
    output: String,
    assembly_kmers: String,
    assembly_fasta: String,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let het_kmers = params.value_of("het_kmers").unwrap();
    let output = params.value_of("output").unwrap();
    let txg_tmp = match params.values_of("txg_mols") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut txg_mols: Vec<String> = Vec::new();
    for x in txg_tmp { txg_mols.push(x.to_string()); }
    let hic_tmp =  match params.values_of("hic_mols") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut hic_mols: Vec<String> = Vec::new();
    for x in hic_tmp { hic_mols.push(x.to_string()); }
    let assembly_kmers = params.value_of("assembly_kmers").unwrap();
    let assembly_fasta = params.value_of("assembly_fasta").unwrap();
    Params{
        het_kmers: het_kmers.to_string(),
        output: output.to_string(),
        txg_mols: txg_mols,
        hic_mols: hic_mols,
        assembly_kmers: assembly_kmers.to_string(),
        assembly_fasta: assembly_fasta.to_string(),
    }
}

