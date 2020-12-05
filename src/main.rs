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
    let kmers = load_kmers(&params.het_kmers);

}

#[derive(PartialEq)]
enum KmerType {
    PairedHet,
    UnpairedHet,
    Homozygous,
}

struct Kmers {
    kmers: HashMap<i32, String>,
    kmer_counts: HashMap<i32, i32>,
    kmer_type: HashMap<i32, KmerType>,
}

fn load_kmers(kmerfile: &String) -> Kmers {
    let mut kmers: HashMap<i32, String> = HashMap::new();
    let mut kmer_id: i32 = 1; // no 0 kmer id as we are using sign for which of the pair
    let mut kmer_type: HashMap<i32, KmerType> = HashMap::new();
    let mut kmer_counts: HashMap<i32, i32> = HashMap::new();
    let reader = File::open(kmerfile).expect("cannot open kmer file");
    let mut reader = BufReader::new(reader);
    let mut buf1 = vec![]; let mut buf2 = vec![];
    let mut buf3 = vec![]; let mut buf4 = vec![];
    loop {
        buf1.clear(); buf2.clear(); buf3.clear(); buf4.clear();

        let bytes1 = reader.read_until(b'\t', &mut buf1).expect("cannot read file");
        if bytes1 == 0 { break; } 
   
        let bytes2 = reader.read_until(b'\t', &mut buf2).expect("cannot read file");
        if bytes2 == 0 { break; } 

        let bytes3 = reader.read_until(b'\t', &mut buf3).expect("cannot read file");
        if bytes3 == 0 { break; } 

        let bytes4 = reader.read_until(b'\n', &mut buf4).expect("cannot read file");
        if bytes4 == 0 { break; } 

        let switchme = std::str::from_utf8(&buf3).unwrap();
        match switchme {
            "HET\t" => {
                kmer_type.insert(kmer_id, KmerType::UnpairedHet);
                kmer_type.insert(-kmer_id, KmerType::UnpairedHet);
                kmers.insert(kmer_id, std::str::from_utf8(&buf1[0..(bytes1-1)]).unwrap().to_string());
                //println!("")
            },
            "HOM\t" => {
                kmer_type.insert(kmer_id, KmerType::Homozygous);
                kmer_type.insert(-kmer_id, KmerType::Homozygous);
            },
            _ => {
                kmer_type.insert(kmer_id, KmerType::PairedHet);
                kmer_type.insert(-kmer_id, KmerType::PairedHet);
                kmer_type.insert(kmer_id+1, KmerType::PairedHet);
                kmer_type.insert(-(kmer_id+1), KmerType::PairedHet);
                kmers.insert(kmer_id, std::str::from_utf8(&buf1[0..(bytes1-1)]).unwrap().to_string());
                kmers.insert(-kmer_id, std::str::from_utf8(&buf1[0..(bytes1-1)]).unwrap().to_string());
                let count = std::str::from_utf8(&buf2[0..(bytes2-1)]).unwrap().to_string().parse::<i32>().unwrap();
                kmer_counts.insert(kmer_id, count);
                kmer_counts.insert(-kmer_id, count);
                kmers.insert(kmer_id+1, std::str::from_utf8(&buf3[0..(bytes3-1)]).unwrap().to_string());
                kmers.insert(-(kmer_id+1), std::str::from_utf8(&buf3[0..(bytes3-1)]).unwrap().to_string());
                let count = std::str::from_utf8(&buf4[0..(bytes4-1)]).unwrap().to_string().parse::<i32>().unwrap();
                kmer_counts.insert(kmer_id+1, count);
                kmer_counts.insert(-(kmer_id+1), count);
            },
        }
        kmer_id += 2;
    }
    //(kmers, kmer_counts)
    Kmers{
        kmers: kmers,
        kmer_counts: kmer_counts,
        kmer_type: kmer_type,
    }
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

