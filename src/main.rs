#[macro_use]
extern crate clap;
extern crate hashbrown;
extern crate rayon;
extern crate phasst_lib;
extern crate rand;

use phasst_lib::{Kmers, load_assembly_kmers, Assembly, HicMols, load_hic};
use rayon::prelude::*;

use rand::Rng;
use rand::rngs::StdRng;
use rand::SeedableRng;

use hashbrown::{HashMap, HashSet};

use clap::{App};

fn main() {
    eprintln!("Welcome to phasst phase!");
    let params = load_params();
    eprintln!("loading kmers");
    let kmers = Kmers::load_kmers(&params.het_kmers);
    //let (_variants, molecules) = load_molecule_kmers(&params.txg_mols, &params.hic_mols, &params.longread_mols, &kmers);
    eprintln!("loading hic kmers");
    let hic_mols = load_hic(&params.hic_mols, &kmers);
    eprintln!("loading assembly kmers");
    let assembly = load_assembly_kmers(&params.assembly_kmers, &kmers);

    eprintln!("finding good loci");
    let variant_contig_order: ContigLoci = good_assembly_loci(&assembly);
    eprintln!("finding good hic reads");
    let hic_links: HashMap<i32, Vec<HIC>> = gather_hic_links(&hic_mols, &variant_contig_order);
    eprintln!("phasing");
    phasst_phase_main(&params, &hic_links, &variant_contig_order);
}

struct ContigLoci {
    kmers: HashMap<i32, (i32, usize)>, // map from kmer id to contig id and position
    loci: HashMap<i32, usize>, // map from contig id to number of loci
}


fn good_assembly_loci(assembly: &Assembly) ->  ContigLoci { // returning a map from kmer id to contig id and position
 
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

    ContigLoci{
        kmers: variant_contig_order,
        loci: loci,
    }
}

struct HIC {
    loci: Vec<usize>,
    alleles: Vec<bool>,
    //kmers: Vec<i32>,
}

fn gather_hic_links(hic_molecules: &HicMols, variant_contig_order: &ContigLoci) -> HashMap<i32, Vec<HIC>> { // returns map from contig id to list of HIC data structures
    let mut hic_mols: HashMap<i32, Vec<HIC>> = HashMap::new();
    let mut total = 0;
    for (contig, _) in variant_contig_order.loci.iter() {
        hic_mols.insert(*contig, Vec::new());
    }
    //let mut counts: [u32;3] = [0,0,0];
    let mut used_count = 0;
    let mut not_assembly = 0;
    let mut diff_contig = 0;
    for mol in hic_molecules.get_hic_molecules() {
        let mut the_contig: Option<i32> = None;
        let mut loci: Vec<usize> = Vec::new();
        let mut alleles: Vec<bool> = Vec::new();
        let mut used: HashSet<i32> = HashSet::new();
        //let mut total = 0;
        //let mut in_assembly = 0;

        for var in mol {
            if used.contains(&var.abs()) { used_count +=1 ; continue; }
            if let Some((contig, order)) = variant_contig_order.kmers.get(&var.abs()) {
                if let Some(chrom) = the_contig {
                    if *contig == chrom {
                        loci.push(*order);
                        if var.abs() % 2 == 0 {
                            alleles.push(false);
                        } else {
                            alleles.push(true);
                        }
                    } else { diff_contig += 1; }
                } else { 
                    the_contig = Some(*contig); 
                    loci.push(*order);
                    if var.abs() % 2 == 0 {
                        alleles.push(false);
                    } else {
                        alleles.push(true);
                    }
                }
            } else { not_assembly += 1; }
            used.insert(var.abs());
        }
        if loci.len() > 1 {
            let contig_mols = hic_mols.entry(the_contig.unwrap()).or_insert(Vec::new());
            contig_mols.push( HIC{loci: loci, alleles: alleles}); 
            total += 1;
        }
    }
    eprintln!("after culling we have {} hic molecules hitting >=2 distinct loci", total);
    eprintln!("why did we lose kmers? overlaps (same kmer twice) {}, no assembly locus {}, cross contig {}", used_count, not_assembly, diff_contig);
    /*
    for (contig, mols) in contig_mols.iter() {
        let mut hic: Vec<HIC> = Vec::new();
        for (_mol, vars) in mols.iter() {
            if vars.len() > 1 {
                let mut loci: Vec<usize> = Vec::new();
                let mut alleles: Vec<bool> = Vec::new();
                for (var, index) in vars {
                    loci.push(*index);
                    if var % 2 == 0 {
                        alleles.push(false);
                    } else {
                        alleles.push(true);
                    }
                }
                hic.push(HIC{alleles: alleles, loci: loci});
            }
        }
        hic_mols.insert(*contig, hic);
    }
    */
    hic_mols
}


struct ThreadData {
    best_total_log_probability: HashMap<i32, f32>, // best log prob for each contig
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

fn phasst_phase_main(params: &Params, hic_links: &HashMap<i32, Vec<HIC>>, contig_loci: &ContigLoci) {
    let seed = [params.seed; 32];
    let mut rng: StdRng = SeedableRng::from_seed(seed);
    let mut threads: Vec<ThreadData> = Vec::new();
    let solves_per_thread = ((params.restarts as f32)/(params.threads as f32)).ceil() as usize;
    for i in 0..params.threads {
        threads.push(ThreadData::from_seed(new_seed(&mut rng), solves_per_thread, i));
    }

    // First lets keep track of how many times each locus is hit by a hic molecule (that also hit something else)
    let mut locus_counts: HashMap<(i32, usize), u32> = HashMap::new();
    for (contig, _loci) in contig_loci.loci.iter() {
        let contig_hic_links = hic_links.get(contig).unwrap();
        for hic_read in contig_hic_links.iter() { 
            for locus in hic_read.loci.iter() {
                let count = locus_counts.entry((*contig, *locus)).or_insert(0);
                *count += 1;
            }
        } 
    }
    threads.par_iter_mut().for_each(|thread_data| {
        for iteration in 0..thread_data.solves_per_thread {
            for (contig, loci) in contig_loci.loci.iter() {
                let cluster_centers: Vec<Vec<f32>> = init_cluster_centers(*loci, params, &mut thread_data.rng);
                let contig_hic_links = hic_links.get(contig).unwrap();
                let (log_loss, cluster_centers) = expectation_maximization(*loci, cluster_centers, &contig_hic_links, params, iteration, thread_data.thread_num);
                let best_log_prob_so_far = thread_data.best_total_log_probability.entry(*contig).or_insert(f32::NEG_INFINITY);
                if &log_loss > best_log_prob_so_far {
                    thread_data.best_total_log_probability.insert(*contig, log_loss);
                }
                eprintln!("thread {} contig {} iteration {} done with {}, best so far {}", 
                    thread_data.thread_num, contig, iteration, log_loss, thread_data.best_total_log_probability.get(contig).unwrap());
                for index in 0..cluster_centers[0].len() {
                    let mut count = 0;
                    if let Some(x) = locus_counts.get(&(*contig, index)){
                        count = *x;
                    }
                    eprintln!("{}\t{}\t{}", cluster_centers[0][index], cluster_centers[1][index], count);
                }
            }
        }
    });
}

fn expectation_maximization(loci: usize, mut cluster_centers: Vec<Vec<f32>>, hic_links: &Vec<HIC>, 
    params: &Params, epoch: usize, thread_num: usize) -> (f32, Vec<Vec<f32>>) {
    if hic_links.len() == 0 {
        eprintln!("no hic links?");
    }
    let mut sums: Vec<Vec<f32>> = Vec::new();
    let mut denoms: Vec<Vec<f32>> = Vec::new();
    
    for cluster in 0..params.ploidy {
        sums.push(Vec::new());
        denoms.push(Vec::new());
        for _index in 0..loci {
            sums[cluster].push(0.5);
            denoms[cluster].push(1.0); // psuedocounts
        }
    }

     // now lets do some EM
    let log_prior: f32 = (1.0/(params.ploidy as f32)).ln();
    let mut iterations = 0;

    let mut total_log_loss = f32::NEG_INFINITY;

    let mut final_log_probabilities: Vec<Vec<f32>> = Vec::new();
    for _read in 0..hic_links.len() {
        final_log_probabilities.push(Vec::new());
    }

    let log_loss_change_limit = 0.000001 * (hic_links.len() as f32); // TODO fiddle with this in testing
    let mut last_log_loss = f32::NEG_INFINITY;
    let mut log_loss_change = 10000.0;
    while log_loss_change > log_loss_change_limit && iterations < 10000 {
        let mut log_bernoulli_loss = 0.0;
        reset_sums_denoms(loci, &mut sums, &mut denoms, params.ploidy);
        for (readdex, hic_read) in hic_links.iter().enumerate() { 
            let log_bernoullis = bernoulli_loss(hic_read, &cluster_centers, log_prior);
            log_bernoulli_loss += log_sum_exp(&log_bernoullis);
            let probabilities = normalize_in_log(&log_bernoullis);

            update_centers_average(&mut sums, &mut denoms, hic_read, &probabilities);
            /*
            eprintln!("hic read {} with probabilities {:?}", readdex, probabilities);
            for cluster in 0..2 {
                for i in 0..hic_read.loci.len() {
                    eprintln!("\tcluster {} locus {} cluster center {} allele {}", 
                        cluster, hic_read.loci[i], cluster_centers[cluster][hic_read.loci[i]], hic_read.alleles[i]);
                }
            }
            
            for cluster in 0..2 {
                for i in 0..hic_read.loci.len() {
                    eprintln!("\tcluster {} locus {} numerator {} denominator {}", 
                        cluster, hic_read.loci[i], sums[cluster][hic_read.loci[i]], denoms[cluster][hic_read.loci[i]]);
                }
            }
            */
            final_log_probabilities[readdex] = log_bernoullis;
        }
        total_log_loss = log_bernoulli_loss;
        log_loss_change = log_bernoulli_loss - last_log_loss;//log_loss - last_log_loss;
        last_log_loss = log_bernoulli_loss;//log_loss;

        update_final(loci, &sums, &denoms, &mut cluster_centers);
        iterations += 1;
        eprintln!("bernoulli\t{}\t{}\t{}\t{}\t{}", thread_num, epoch, iterations,  log_bernoulli_loss, log_loss_change);
    }
    (total_log_loss, cluster_centers)
}


fn update_final(loci: usize, sums: &Vec<Vec<f32>>, denoms: &Vec<Vec<f32>>, cluster_centers: &mut Vec<Vec<f32>>) {
    for locus in 0..loci {
        for cluster in 0..sums.len() {
            let update = sums[cluster][locus]/denoms[cluster][locus];
            cluster_centers[cluster][locus] = update.min(0.99).max(0.01);//max(0.0001, min(0.9999, update));
        }
    }
}

fn update_centers_average(sums: &mut Vec<Vec<f32>>, denoms: &mut Vec<Vec<f32>>, hic_read: &HIC, probabilities: &Vec<f32>) {
    for locus in 0..hic_read.loci.len() {
        for (cluster, probability) in probabilities.iter().enumerate() {
            match hic_read.alleles[locus] {
                false => sums[cluster][hic_read.loci[locus]] += 0.0,
                true => sums[cluster][hic_read.loci[locus]] += probability,
            }
            denoms[cluster][hic_read.loci[locus]] += probabilities[cluster];
        }
    }
}

fn normalize_in_log(log_probs: &Vec<f32>) -> Vec<f32> { // takes in a log_probability vector and converts it to a normalized probability
    let mut normalized_probabilities: Vec<f32> = Vec::new();
    let sum = log_sum_exp(log_probs);
    for i in 0..log_probs.len() {
        normalized_probabilities.push((log_probs[i]-sum).exp());
    }
    normalized_probabilities
}

fn bernoulli_loss(hic_read: &HIC, cluster_centers: &Vec<Vec<f32>>, log_prior: f32) -> Vec<f32> {
    let mut log_probabilities: Vec<f32> = Vec::new();
    for (cluster, center) in cluster_centers.iter().enumerate() {
        log_probabilities.push(log_prior);
        for (locus_index, locus) in hic_read.loci.iter().enumerate() {
            log_probabilities[cluster] += ln_bernoulli(hic_read.alleles[locus_index], center[*locus]);
        }
    }
    
    log_probabilities
}

fn ln_bernoulli(kmer: bool, p: f32) -> f32 {
    match kmer {
        false => (1.0 - p).ln(),
        true => p.ln(),
    }
}

fn log_sum_exp(p: &Vec<f32>) -> f32{
    let max_p: f32 = p.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    let sum_rst: f32 = p.iter().map(|x| (x - max_p).exp()).sum();
    max_p + sum_rst.ln()
}

fn reset_sums_denoms(loci: usize, sums: &mut Vec<Vec<f32>>, 
    denoms: &mut Vec<Vec<f32>>, num_clusters: usize) {
    for cluster in 0..num_clusters {
        for index in 0..loci {
            sums[cluster][index] = 0.50;
            denoms[cluster][index] = 1.0;
        }
    }
}

fn init_cluster_centers(loci: usize, params: &Params, rng: &mut StdRng) -> Vec<Vec<f32>> {
    let mut centers: Vec<Vec<f32>> = Vec::new();
    for cluster in 0..params.ploidy {
        centers.push(Vec::new());
        for _ in 0..loci {
            centers[cluster].push(rng.gen::<f32>().min(0.99).max(0.01));
        }
    }
    centers
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

