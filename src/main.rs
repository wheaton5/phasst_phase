#[macro_use]
extern crate clap;
extern crate hashbrown;
extern crate rayon;
extern crate phasst_lib;
extern crate rand;
extern crate disjoint_set;
extern crate statrs;

use phasst_lib::{Kmers, load_assembly_kmers, Assembly, HicMols, load_hic};
use rayon::prelude::*;

use rand::Rng;
use rand::rngs::StdRng;
use rand::SeedableRng;

use hashbrown::{HashMap, HashSet};
use disjoint_set::DisjointSet;
//use statrs::distribution::{Beta, Continuous};
use statrs::function::{beta};

use clap::{App};

const LONG_RANGE_HIC: usize = 15000;
const LONG_RANGE_HIC_WEIGHTING: f32 = 100.0;
const MIN_ALLELE_FRACTION_HIC: f32 = 0.15;

struct ContigLoci {
    kmers: HashMap<i32, (i32, usize, usize)>, // map from kmer id to contig id and position
    loci: HashMap<i32, usize>, // map from contig id to number of loci
}

struct HIC {
    loci: Vec<usize>,
    alleles: Vec<bool>,
    long_weighting: f32,
    //kmers: Vec<i32>,
}

struct ThreadData {
    best_total_log_probability: HashMap<i32, f32>, // best log prob for each contig
    rng: StdRng,
    solves_per_thread: usize,
    thread_num: usize,
}



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
    let allele_fractions = get_allele_fractions(&hic_mols); // MAYBE ADD LINKED READ AND LONG READ to this?
    let variant_contig_order: ContigLoci = good_assembly_loci(&assembly, &allele_fractions);
    eprintln!("finding good hic reads");
    let (hic_links, long_hic_links) = gather_hic_links(&hic_mols, &variant_contig_order);
    eprintln!("phasing");
    let mut connected_components = get_connected_components(&hic_links, &variant_contig_order, params.min_hic_links);
    phasst_phase_main(&params, &hic_links, &long_hic_links, &variant_contig_order);
}

fn get_connected_components(hic_links: &HashMap<i32, Vec<HIC>>, variant_contig_order: &ContigLoci, min_links: u32) -> 
    HashMap<i32, DisjointSet<usize>> {
    let mut components: HashMap<i32, DisjointSet<usize>> = HashMap::new();
    let mut edges: HashMap<i32, HashMap<(usize, usize), u32>> = HashMap::new();
    for (mol, hic_mols) in hic_links.iter() {
        let mut disjoint_set: DisjointSet<usize> = DisjointSet::new(); 
        for i in 0..*variant_contig_order.loci.get(mol).unwrap() {
            disjoint_set.make_set(i);
        }
        let mut contig_edges: HashMap<(usize, usize), u32> = HashMap::new();
        for hic_mol in hic_mols {
            for i in 0..hic_mol.loci.len() {
                for j in (i+1)..hic_mol.loci.len() {
                    let min = hic_mol.loci[i].min(hic_mol.loci[j]);
                    let max = hic_mol.loci[i].max(hic_mol.loci[j]);
                    let count = contig_edges.entry((min, max)).or_insert(0);
                    *count += 1;
                }
            }
        }
        components.insert(*mol, disjoint_set);
        edges.insert(*mol, contig_edges);
    }

    for (contig, contig_edges) in edges {
        let disjoint_set = components.get_mut(&contig).unwrap();
        for ((locus1, locus2), count) in contig_edges.iter() {
            if *count > min_links {
                disjoint_set.union(*locus1, *locus2).expect("cannot merge");
            }
        }
    }
    components
}

fn get_allele_fractions(hic_mols: &HicMols) -> HashMap<i32, f32> {
    let mut allele_fractions: HashMap<i32, f32> = HashMap::new();
    let mut allele_counts: HashMap<i32, [u32; 2]> = HashMap::new();
    for mol in hic_mols.get_hic_molecules() {
        for var in mol {
            let canonical = var.abs().min(Kmers::pair(var.abs()));
            let count = allele_counts.entry(canonical).or_insert([0;2]);
            if var.abs() % 2 == 0 { count[0] += 1; } else { count[1] += 1; }
        }
    }
    for (canonical, counts) in allele_counts {
        let fraction = (counts[0].min(counts[1]) as f32)/((counts[0]+ counts[1]) as f32);
        allele_fractions.insert(canonical, fraction);
        allele_fractions.insert(Kmers::pair(canonical), fraction);
    }
    allele_fractions
}

fn good_assembly_loci(assembly: &Assembly, allele_fractions: &HashMap<i32, f32>) ->  ContigLoci { // returning a map from kmer id to contig id and position
    let mut variant_contig_order: HashMap<i32, (i32, usize, usize)> = HashMap::new();

    let mut contig_positions: HashMap<i32, Vec<(usize, i32, i32)>> = HashMap::new();
    for (kmer, (contig, num, _order, position)) in assembly.variants.iter() {
        if assembly.variants.contains_key(&Kmers::pair(*kmer)) { continue; } // we see both ref and alt in assembly, skip
        if let Some(fraction) = allele_fractions.get(&kmer.abs()) {
            if *fraction < MIN_ALLELE_FRACTION_HIC { continue; }
        } else { continue; }

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
        for (index, (position, kmer, contig)) in poses.iter().enumerate() {
            variant_contig_order.insert(*kmer, (*contig, index, *position));
            variant_contig_order.insert(Kmers::pair(*kmer), (*contig, index, *position));
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

fn gather_hic_links(hic_molecules: &HicMols, variant_contig_order: &ContigLoci) -> 
        (HashMap<i32, Vec<HIC>>, HashMap<i32, Vec<HIC>>) { // returns map from contig id to list of HIC data structures
    let mut hic_mols: HashMap<i32, Vec<HIC>> = HashMap::new();
    let mut long_hic_mols: HashMap<i32, Vec<HIC>> = HashMap::new();
    let mut total = 0;
    let mut total_long_links = 0;
    for (contig, _) in variant_contig_order.loci.iter() {
        hic_mols.insert(*contig, Vec::new());
        long_hic_mols.insert(*contig, Vec::new());
    }
    
    let mut used_count = 0;
    let mut not_assembly = 0;
    let mut diff_contig = 0;
    for mol in hic_molecules.get_hic_molecules() {
        let mut the_contig: Option<i32> = None;
        let mut loci: Vec<usize> = Vec::new();
        let mut alleles: Vec<bool> = Vec::new();
        let mut used: HashSet<i32> = HashSet::new();
        let mut min: usize = std::usize::MAX;
        let mut max: usize = 0;

        for var in mol {
            if used.contains(&var.abs()) { used_count +=1 ; continue; }
            if let Some((contig, order, position)) = variant_contig_order.kmers.get(&var.abs()) {
                if let Some(chrom) = the_contig {
                    min = min.min(*position);
                    max = max.max(*position);
                    if *contig == chrom {
                        loci.push(*order);
                        if var.abs() % 2 == 0 {
                            alleles.push(false);
                        } else {
                            alleles.push(true);
                        }
                    } else { diff_contig += 1; }
                } else { 
                    min = min.min(*position);
                    max = max.max(*position);
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
            let mut long_range = 1.0;
            if max - min > LONG_RANGE_HIC {
                long_range = LONG_RANGE_HIC_WEIGHTING;
                let long_contig_mols = long_hic_mols.entry(the_contig.unwrap()).or_insert(Vec::new());
                let mut long_loci: Vec<usize> = Vec::new();
                let mut long_alleles: Vec<bool> = Vec::new();
                for index in 0..loci.len(){
                    long_loci.push(loci[index]);
                    long_alleles.push(alleles[index]);
                }
                long_contig_mols.push( HIC{loci: long_loci, alleles: long_alleles, long_weighting: LONG_RANGE_HIC_WEIGHTING} );
                total_long_links += 1;
            }
            contig_mols.push( HIC{loci: loci, alleles: alleles, long_weighting: long_range} ); 
            total += 1;
        }
    }
    eprintln!("after culling we have {} hic molecules hitting >=2 distinct loci", total);
    eprintln!("why did we lose kmers? overlaps (same kmer twice) {}, no assembly locus {}, cross contig {}", used_count, not_assembly, diff_contig);
    eprintln!("num long hic links {}", total_long_links);
    (hic_mols, long_hic_mols)
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

fn phasst_phase_main(params: &Params, hic_links: &HashMap<i32, Vec<HIC>>, long_hic_links: &HashMap<i32, Vec<HIC>>, 
        contig_loci: &ContigLoci) {
    let seed = [params.seed; 32];
    let mut rng: StdRng = SeedableRng::from_seed(seed);
    let mut threads: Vec<ThreadData> = Vec::new();
    let solves_per_thread = ((params.restarts as f32)/(params.threads as f32)).ceil() as usize;
    for i in 0..params.threads {
        threads.push(ThreadData::from_seed(new_seed(&mut rng), solves_per_thread, i));
    }

    // First lets keep track of how many times each locus is hit by a hic molecule (that also hit something else)
    let mut locus_counts: HashMap<(i32, usize), u32> = HashMap::new();
    let mut contig_locus_hic_mols: HashMap<i32, HashMap<usize, Vec<usize>>> = HashMap::new(); // TODO REMOVE DEBUG map from contig to a map from locus to vec of hic mols
    for (contig, loci) in contig_loci.loci.iter() {
        let contig_hic_links = hic_links.get(contig).unwrap();
        let locus_hic = contig_locus_hic_mols.entry(*contig).or_insert(HashMap::new()); //TODO REMOVE DEBUG
        for (read_index, hic_read) in contig_hic_links.iter().enumerate() {
            for locus in hic_read.loci.iter() {
                let hic_mols = locus_hic.entry(*locus).or_insert(Vec::new()); // TODO REMOVE DEBUG
                hic_mols.push(read_index); //TODO REMOVE DEBUG
                let count = locus_counts.entry((*contig, *locus)).or_insert(0);
                *count += 1;
            }
        } 
    }
    threads.par_iter_mut().for_each(|thread_data| {
        for iteration in 0..thread_data.solves_per_thread {
            for (contig, loci) in contig_loci.loci.iter() {
                let locus_hic_mols = contig_locus_hic_mols.get(contig).unwrap(); // TODO REMOVE DEBUG
                let cluster_centers: Vec<Vec<f32>> = init_cluster_centers(*loci, params, &mut thread_data.rng);
                let contig_hic_links = hic_links.get(contig).unwrap();
                let long_contig_hic_links = long_hic_links.get(contig).unwrap();
                eprintln!("solve with LONG LINKS ONLY");
                let (log_loss, cluster_centers, hic_probabilities) =  // first solve with long links only
                    expectation_maximization(*loci, cluster_centers, &long_contig_hic_links, params, iteration, thread_data.thread_num);
                eprintln!("ALL HIC LINKS");
                let (log_loss, cluster_centers, hic_probabilities) = 
                    expectation_maximization(*loci, cluster_centers, &contig_hic_links, params, iteration, thread_data.thread_num);
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
                        // now i have hic_probabilities which gives me the probabilities of each hic molecule to each cluster
                        // and i have contig_locus_hic_mols and now locus_hic_mols which gives me a map
                        // from locus to a vec of hic mol indexes
                    // TODO DEBUG HERE UNTIL NeXT DEBUG
                    if let Some(hic_moldexes) = locus_hic_mols.get(&index) {
                        for hicdex in hic_moldexes {
                            let mut centers: Vec<(f32, f32)> = Vec::new();
                            for locus in contig_hic_links[*hicdex].loci.iter() {
                                centers.push((cluster_centers[0][*locus], cluster_centers[1][*locus]));
                            }
                            let probs = &hic_probabilities[*hicdex];
                            eprintln!("\thicread\t{}\tloci\t{:?}\talleles\t{:?}\tclusters\t{:?}\tprobs{:?}",
                                hicdex, contig_hic_links[*hicdex].loci, contig_hic_links[*hicdex].alleles, centers, probs);
                        }
                    }
                    // END DEBUG
                }
            }
        }
    });
}

fn get_beta_priors(cluster_centers: &Vec<Vec<f32>>) -> Vec<Vec<f32>> {
    let mut prior: Vec<Vec<f32>> = Vec::new();
    for i in 0..cluster_centers[0].len() {
        let mut ones = Vec::new();
        for j in 0..cluster_centers.len() {
            ones.push(1.0);
        }
        prior.push(ones);
    }
    prior
}

fn fill_beta_priors(priors: &mut Vec<Vec<f32>>) {
    for i in 0..priors.len() {
        for j in 0..priors[i].len() {
            priors[i][j] = 1.0;
        }
    }
}

fn transfer_beta_priors(priors: &mut Vec<Vec<f32>>, next_priors: &Vec<Vec<f32>>) {
    for i in 0..priors.len() {
        for j in 0..priors[i].len() {
            priors[i][j] = next_priors[i][j];
        }
    }
}

fn expectation_maximization(loci: usize, mut cluster_centers: Vec<Vec<f32>>, hic_links: &Vec<HIC>, 
        params: &Params, epoch: usize, thread_num: usize) -> (f32, Vec<Vec<f32>>,
        Vec<Vec<f32>>) { // this is a vector of the hic molecules probabilities to each cluster
    if hic_links.len() == 0 {
        eprintln!("no hic links?");
    }
    let mut sums: Vec<Vec<f32>> = Vec::new();
    let mut denoms: Vec<Vec<f32>> = Vec::new();
    
    for cluster in 0..params.ploidy {
        sums.push(Vec::new());
        denoms.push(Vec::new());
        for _index in 0..loci {
            sums[cluster].push(0.1);
            denoms[cluster].push(0.2); // psuedocounts
        }
    }

     // now lets do some EM
    let log_prior: f32 = (1.0/(params.ploidy as f32)).ln();
    let mut iterations = 0;

    let mut total_log_loss = f32::NEG_INFINITY;

    let mut final_log_probabilities: Vec<Vec<f32>> = Vec::new();
    let mut alphas: Vec<Vec<f32>> = get_beta_priors(&cluster_centers);
    let mut betas: Vec<Vec<f32>> = get_beta_priors(&cluster_centers);
    let mut alphas_next: Vec<Vec<f32>> = get_beta_priors(&cluster_centers);
    let mut betas_next: Vec<Vec<f32>> = get_beta_priors(&cluster_centers);
    for _read in 0..hic_links.len() {
        final_log_probabilities.push(Vec::new());
    }

    let log_loss_change_limit = 0.000001 * (hic_links.len() as f32); // TODO fiddle with this in testing
    let mut last_log_loss = f32::NEG_INFINITY;
    let mut log_loss_change = 10000.0;
    let mut hic_probabilities: Vec<Vec<f32>> = Vec::new(); // TODO REMOVE DEBUG

    //while log_loss_change > log_loss_change_limit && iterations < 100 {
    while iterations < 150 { // TODO figure out something better here
        hic_probabilities.clear(); // TODO REMOVE DEBUG
        let mut log_likelihood = 0.0;
        fill_beta_priors(&mut alphas_next);
        fill_beta_priors(&mut betas_next);
        reset_sums_denoms(loci, &mut sums, &mut denoms, params.ploidy);
        for (readdex, hic_read) in hic_links.iter().enumerate() {
            let log_likelihoods;
            if iterations == 0 {
                log_likelihoods = bernoulli_likelihood(hic_read, &cluster_centers, log_prior);
            } else {
                log_likelihoods = beta_likelihood(hic_read, &cluster_centers, log_prior, &alphas, &betas);
            }
            log_likelihood += log_sum_exp(&log_likelihoods);
            let probabilities = normalize_in_log(&log_likelihoods);
            for (index, locus) in hic_read.loci.iter().enumerate() {
                for (cluster, probability) in probabilities.iter().enumerate() {
                    if hic_read.alleles[index] == true {
                        betas_next[*locus][cluster] += probability;
                    } else {
                        alphas_next[*locus][cluster] += probability;
                    }
                }
            }
            
            update_sums_denoms(&mut sums, &mut denoms, hic_read, &probabilities);
            hic_probabilities.push(probabilities);
            final_log_probabilities[readdex] = log_likelihoods;
        }
        total_log_loss = log_likelihood;
        log_loss_change = log_likelihood - last_log_loss;
        last_log_loss = log_likelihood;
        transfer_beta_priors(&mut alphas, &alphas_next);
        transfer_beta_priors(&mut betas, &betas_next);
        update_cluster_centers(loci, &sums, &denoms, &mut cluster_centers);
        iterations += 1;
        eprintln!("bernoulli\t{}\t{}\t{}\t{}\t{}", thread_num, epoch, iterations,  log_likelihood, log_loss_change);
    }
    (total_log_loss, cluster_centers, hic_probabilities)
}


fn update_cluster_centers(loci: usize, sums: &Vec<Vec<f32>>, denoms: &Vec<Vec<f32>>, cluster_centers: &mut Vec<Vec<f32>>) {
    for locus in 0..loci {
        for cluster in 0..sums.len() {
            if denoms[cluster][locus] > 1.0 {
                let update = sums[cluster][locus]/denoms[cluster][locus];
                cluster_centers[cluster][locus] = update.min(0.9999).max(0.0001);//max(0.0001, min(0.9999, update));
            }
        }
    }
}

fn update_sums_denoms(sums: &mut Vec<Vec<f32>>, denoms: &mut Vec<Vec<f32>>, hic_read: &HIC, probabilities: &Vec<f32>) {
    for locus in 0..hic_read.loci.len() {
        for (cluster, probability) in probabilities.iter().enumerate() {
            match hic_read.alleles[locus] {
                false => sums[cluster][hic_read.loci[locus]] += 0.0,
                true => sums[cluster][hic_read.loci[locus]] += probability,// * hic_read.long_weighting,
            }
            denoms[cluster][hic_read.loci[locus]] += probabilities[cluster];// * hic_read.long_weighting;
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

fn bernoulli_likelihood(hic_read: &HIC, cluster_centers: &Vec<Vec<f32>>, log_prior: f32) -> Vec<f32> {
    let mut log_probabilities: Vec<f32> = Vec::new();
    for (cluster, center) in cluster_centers.iter().enumerate() {
        log_probabilities.push(log_prior);
        for (locus_index, locus) in hic_read.loci.iter().enumerate() {
            log_probabilities[cluster] += ln_bernoulli(hic_read.alleles[locus_index], center[*locus]);
        }
    }
    
    log_probabilities
}

fn beta_likelihood(hic_read: &HIC, cluster_centers: &Vec<Vec<f32>>, log_prior: f32, alphas: &Vec<Vec<f32>>, betas: &Vec<Vec<f32>>) -> Vec<f32> {
    let mut log_probabilities: Vec<f32> = Vec::new();
    for (cluster, center) in cluster_centers.iter().enumerate() {
        log_probabilities.push(log_prior);
        for (locus_index, locus) in hic_read.loci.iter().enumerate() {
            log_probabilities[cluster] += ln_beta(hic_read.alleles[locus_index], center[*locus] as f64, 
                alphas[*locus][cluster] as f64, betas[*locus][cluster] as f64);
        }
    }
    log_probabilities
}

fn ln_beta(allele: bool, p: f64, alpha: f64, beta: f64) -> f32 {
    let p = match allele {
        false => 1.0 - p,
        true => p,
    };
    ((alpha - 1.0)*p + (beta - 1.0)*(1.0 - p) - beta::ln_beta(alpha, beta)) as f32
}

fn ln_bernoulli(allele: bool, p: f32) -> f32 {
    match allele {
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
            sums[cluster][index] = 0.10;
            denoms[cluster][index] = 0.20;
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
    min_hic_links: u32,
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

    let min_hic_links = params.value_of("min_hic_links").unwrap_or("4");
    let min_hic_links = min_hic_links.to_string().parse::<u32>().unwrap();

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
        min_hic_links: min_hic_links,
    }
}

//fn unwrap_param(name: &str, default: Option<&str>, )

