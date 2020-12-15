#[macro_use]
extern crate clap;
extern crate hashbrown;
extern crate rayon;
extern crate phasst_lib;
extern crate rand;

use phasst_lib::{Kmers, load_molecule_kmers, load_assembly_kmers, Variants, Molecules, Assembly};
use rayon::prelude::*;

use rand::Rng;
use rand::rngs::StdRng;
use rand::SeedableRng;

use hashbrown::{HashMap};

use clap::{App};

fn main() {
    println!("Welcome to phasst phase!");
    let params = load_params();
    let kmers = Kmers::load_kmers(&params.het_kmers);
    let (variants, molecules) = load_molecule_kmers(&params.txg_mols, &params.hic_mols, &params.longread_mols, &kmers);
    let assembly = load_assembly_kmers(&params.assembly_kmers, &kmers);

    
    let variant_contig_order: Contig_Loci = good_assembly_loci(&assembly);
    let hic_links: HashMap<i32, Vec<HIC>> = gather_hic_links(&molecules, &variant_contig_order);
    phasst_phase_main(&params, &hic_links);
}

struct Contig_Loci {
    kmers: HashMap<i32, (i32, usize)>, // map from kmer id to contig id and position
    loci: HashMap<i32, usize>, // map from contig id to number of loci
}


fn good_assembly_loci(assembly: &Assembly) ->  Contig_Loci { // returning a map from kmer id to contig id and position
 
    let mut variant_contig_order: HashMap<i32, (i32, usize)> = HashMap::new();

    let mut contig_positions: HashMap<i32, Vec<(usize, i32, i32)>> = HashMap::new();
    for (kmer, (contig, num, _order, position)) in assembly.variants.iter() {
        if assembly.variants.contains_key(&Kmers::pair(*kmer)) { continue; } // we see both ref and alt in assembly, skip
        if *num > 1 { continue; } // we see this kmer multiple times in the assembly, skip
        let positions = contig_positions.entry(*contig).or_insert(Vec::new());
        positions.push((*position, *kmer, *contig));
    }

    for (_contig, positions) in contig_positions.iter() {
        let mut poses: Vec<(usize, i32, i32)> = Vec::new();
        for (pos, kmer, contig) in positions {
            poses.push((*pos, *kmer, *contig));
        }
        poses.sort();
        for (index, (_position, kmer, contig)) in poses.iter().enumerate() {
            variant_contig_order.insert(*kmer, (*contig, index));
            variant_contig_order.insert(Kmers::pair(*kmer), (*contig, index));
        }
    }

    let mut loci: HashMap<i32, usize> = HashMap::new();
    for (contig, positions) in contig_positions.iter() { 
        loci.insert(*contig, positions.len());
    }

    Contig_Loci{
        kmers: variant_contig_order,
        loci: loci,
    }
}

struct HIC {
    loci: Vec<usize>,
    alleles: Vec<f32>,
}

fn gather_hic_links(molecules: &Molecules, variant_contig_order: &Contig_Loci) -> HashMap<i32, Vec<HIC>> { // returns map from contig id to list of HIC data structures
    let mut hic_mols: HashMap<i32, Vec<HIC>> = HashMap::new();

    let mut contig_mols: HashMap<i32, HashMap<i32, Vec<(i32, usize)>>> = HashMap::new();
    for mol in molecules.get_hic_molecules() {
        for var in molecules.get_hic_variants(*mol) {
            if let Some((contig, order)) = variant_contig_order.kmers.get(&var.abs()) {
                let mols = contig_mols.entry(*contig).or_insert(HashMap::new());
                let vars = mols.entry(*mol).or_insert(Vec::new());
                vars.push((var.abs(), *order));
            }
        }
    }
    for (contig, mols) in contig_mols.iter() {
        let mut hic: Vec<HIC> = Vec::new();
        for (mol, vars) in mols.iter() {
            if vars.len() > 1 {
                let mut loci: Vec<usize> = Vec::new();
                let mut alleles: Vec<f32> = Vec::new();
                for (var, index) in vars {
                    loci.push(*index);
                    if var % 2 == 0 {
                        alleles.push(0.0);
                    } else {
                        alleles.push(1.0);
                    }
                }
                hic.push(HIC{alleles: alleles, loci: loci});
            }
        }
        hic_mols.insert(*contig, hic);
    }

    hic_mols
}


struct ThreadData {
    best_total_log_probability: HashMap<u32, f32>, // best log prob for each contig
    rng: StdRng,
    solves_per_thread: usize,
    thread_num: usize,
}

impl ThreadData {
    fn from_seed(seed: [u8; 32], solves_per_thread: usize, thread_num: usize) -> ThreadData {
        ThreadData {
            best_total_log_probability: HashMap::new(),
            rng: SeedableRng::from_seed(seed),
            solves_per_thread: solves_per_thread,
            thread_num: thread_num,
        }
    }
}

fn phasst_phase_main(params: &Params, hic_links: &HashMap<i32, Vec<HIC>>) {
    let seed = [params.seed; 32];
    let mut rng: StdRng = SeedableRng::from_seed(seed);
    let mut threads: Vec<ThreadData> = Vec::new();
    let solves_per_thread = ((params.restarts as f32)/(params.threads as f32)).ceil() as usize;
    for i in 0..params.threads {
        threads.push(ThreadData::from_seed(new_seed(&mut rng), solves_per_thread, i));
    }
    threads.par_iter_mut().for_each(|thread_data| {
        let cluster_centers: HashMap<i32, Vec<Vec<f32>>> = init_cluster_centers(loci_used, &cell_data, params, &mut thread_data.rng, &locus_to_index);
    });
}

fn init_cluster_centers(loci_used: usize, cell_data: &Vec<CellData>, params: &Params, rng: &mut StdRng, locus_to_index: &HashMap<usize, usize>) -> Vec<Vec<f32>> {

}

fn new_seed(rng: &mut StdRng) -> [u8; 32] {
    let mut seed = [0; 32];
    for i in 0..32 {
        seed[i] = rng.gen::<u8>();
    }
    seed
}


#[derive(Clone)]
struct Params {
    het_kmers: String,
    txg_mols: Option<Vec<String>>,
    hic_mols: Option<Vec<String>>,
    longread_mols: Option<Vec<String>>,
    output: String,
    assembly_kmers: String,
    assembly_fasta: String,
    threads: usize,
    seed: u8,
    ploidy: usize,
    restarts: u32,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    
    let het_kmers = params.value_of("het_kmers").unwrap();
    let output = params.value_of("output").unwrap();
    let txg_tmp: Option<Vec<&str>> = match params.values_of("linked_read_mols") {
        Some(x) => Some(x.collect()),
        None => None,
    };
    let mut txg_mols: Option<Vec<String>> = None;
    if let Some(mols) = txg_tmp {
        let mut tmp = Vec::new();
        for mol in mols {
            tmp.push(mol.to_string());
        }
        txg_mols = Some(tmp);
    }
    
    
    let hic_tmp: Option<Vec<&str>> =  match params.values_of("hic_mols") {
        Some(x) => Some(x.collect()),
        None => None,
    };
    let mut hic_mols = None;
    if let Some(mols) = hic_tmp {
        let mut tmp: Vec<String> = Vec::new();
        for x in mols { tmp.push(x.to_string()); }
        hic_mols = Some(tmp);
    }

    let long_tmp: Option<Vec<&str>> =  match params.values_of("longread_mols") {
        Some(x) => Some(x.collect()),
        None => None,
    };
    let mut long_mols = None;
    if let Some(mols) = long_tmp {
        let mut tmp: Vec<String> = Vec::new();
        for x in mols { tmp.push(x.to_string()); }
        long_mols = Some(tmp);
    }

    let threads = params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<usize>().unwrap();

    let seed = params.value_of("seed").unwrap_or("4"); // 4 is guarranteed random by dice roll https://xkcd.com/221/
    let seed = seed.to_string().parse::<u8>().unwrap();

    let restarts = params.value_of("restarts").unwrap_or("10");
    let restarts = restarts.to_string().parse::<u32>().unwrap();
    
    let assembly_kmers = params.value_of("assembly_kmers").unwrap();
    let assembly_fasta = params.value_of("assembly_fasta").unwrap();

    let ploidy = params.value_of("ploidy").unwrap_or("2");
    let ploidy = ploidy.to_string().parse::<usize>().unwrap();
    Params{
        het_kmers: het_kmers.to_string(),
        output: output.to_string(),
        txg_mols: txg_mols,
        hic_mols: hic_mols,
        longread_mols: long_mols,
        assembly_kmers: assembly_kmers.to_string(),
        assembly_fasta: assembly_fasta.to_string(),
        threads: threads,
        seed: seed,
        restarts: restarts,
        ploidy: ploidy,
    }
}

