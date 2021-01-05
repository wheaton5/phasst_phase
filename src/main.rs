#[macro_use]
extern crate clap;
extern crate hashbrown;
extern crate rayon;
extern crate phasst_lib;
extern crate rand;
extern crate disjoint_set;

use phasst_lib::{Kmers, load_assembly_kmers, Assembly, HicMols, load_hic};
use rayon::prelude::*;

use rand::Rng;
use rand::rngs::StdRng;
use rand::SeedableRng;

use hashbrown::{HashMap, HashSet};
use disjoint_set::DisjointSet;

use clap::{App};

const LONG_RANGE_HIC: usize = 15000;
const LONG_RANGE_HIC_WEIGHTING: f32 = 100.0;
const MIN_ALLELE_FRACTION_HIC: f32 = 0.15;
const READ_DEBUG: bool = false;

struct ContigLoci {
    kmers: HashMap<i32, ContigLocus>, // map from kmer id to contig id and position and which allele assembly had
    loci: HashMap<i32, Vec<ContigLocus>>, // map from contig id to number of loci
}

#[derive(Debug, Clone, Copy)]
struct ContigLocus {
    contig_id: i32,
    index: usize,
    position: usize,
    allele: Allele,
}

#[derive(Debug, Clone, Copy)]
enum Allele {
    reference, alternate,
}

struct HIC {
    loci: Vec<usize>,
    alleles: Vec<Allele>,
    long_weighting: f32,
    //kmers: Vec<i32>,
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
    let bad_alleles = get_bad_alleles(&hic_mols);
    let variant_contig_order: ContigLoci = good_assembly_loci(&assembly, &allele_fractions, &bad_alleles);
    eprintln!("finding good hic reads");
    let (hic_links, long_hic_links) = gather_hic_links(&hic_mols, &variant_contig_order);
    eprintln!("phasing");
    let mut connected_components = get_connected_components(&hic_links, &variant_contig_order, params.min_hic_links);
    phasst_phase_main(&params, &hic_links, &long_hic_links, &variant_contig_order);
}

fn allele(kmer: i32) -> Allele {
    match kmer.abs() % 2 == 0 {
        true => Allele::reference,
        false => Allele::alternate,
    }
}

fn get_bad_alleles(hic_mols: &HicMols) -> HashSet<i32> {
    let mut bad: HashSet<i32> = HashSet::new();
    for mol in hic_mols.get_hic_molecules() {
        for i in 0..mol.len() {
            for j in (i+1)..mol.len() {
                if Kmers::pair(mol[i].abs()) == mol[j].abs() {
                    bad.insert(mol[i].abs());
                    bad.insert(mol[j].abs());
                }
            }
        }
    }
    bad
}

fn get_connected_components(hic_links: &HashMap<i32, Vec<HIC>>, variant_contig_order: &ContigLoci, min_links: u32) -> 
    HashMap<i32, DisjointSet<usize>> {
    let mut components: HashMap<i32, DisjointSet<usize>> = HashMap::new();
    let mut edges: HashMap<i32, HashMap<(usize, usize), u32>> = HashMap::new();
    for (mol, hic_mols) in hic_links.iter() {
        let mut disjoint_set: DisjointSet<usize> = DisjointSet::new(); 
        for i in 0..variant_contig_order.loci.get(mol).unwrap().len() {
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



fn good_assembly_loci(assembly: &Assembly, allele_fractions: &HashMap<i32, f32>, bad_alleles: &HashSet<i32>) ->  ContigLoci { // returning a map from kmer id to contig id and position
    let mut variant_contig_order: HashMap<i32, ContigLocus> = HashMap::new();

    let mut contig_positions: HashMap<i32, Vec<(usize, i32, i32)>> = HashMap::new();
    for (kmer, (contig, num, _order, position)) in assembly.variants.iter() {
        if assembly.variants.contains_key(&Kmers::pair(*kmer)) { continue; } // we see both ref and alt in assembly, skip
        if let Some(fraction) = allele_fractions.get(&kmer.abs()) {
            if *fraction < MIN_ALLELE_FRACTION_HIC { continue; }
            if bad_alleles.contains(kmer) { continue; }
        } else { continue; }

        if *num > 1 { continue; } // we see this kmer multiple times in the assembly, skip
        let positions = contig_positions.entry(*contig).or_insert(Vec::new());
        positions.push((*position, *kmer, *contig));
    }

    let mut loci: HashMap<i32, Vec<ContigLocus>> = HashMap::new();
    for (contig, positions) in contig_positions.iter() {
        let mut poses: Vec<(usize, i32, i32)> = Vec::new();
        for (pos, kmer, contig) in positions {
            poses.push((*pos, *kmer, *contig));
        }
        poses.sort();
        let contig_loci = loci.entry(*contig).or_insert(Vec::new());
        for (index, (position, kmer, contig)) in poses.iter().enumerate() {
            let locus = ContigLocus{
                contig_id: *contig,
                index: index,
                position: *position,
                allele: allele(*kmer),
            };
            variant_contig_order.insert(*kmer, locus);
            variant_contig_order.insert(Kmers::pair(*kmer), locus);
            contig_loci.push(locus);
        }

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
        let mut alleles: Vec<Allele> = Vec::new();
        let mut used: HashSet<i32> = HashSet::new();
        let mut min: usize = std::usize::MAX;
        let mut max: usize = 0;

        for var in mol {
            if used.contains(&var.abs()) { used_count +=1 ; continue; }
            if let Some(locus) = variant_contig_order.kmers.get(&var.abs()) {
                if let Some(chrom) = the_contig {
                    min = min.min(locus.position);
                    max = max.max(locus.position);
                    if locus.contig_id == chrom {
                        loci.push(locus.index);
                        alleles.push(allele(*var));
                    } else { diff_contig += 1; }
                } else { 
                    min = min.min(locus.position);
                    max = max.max(locus.position);
                    the_contig = Some(locus.contig_id); 
                    loci.push(locus.index);
                    alleles.push(allele(*var));
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
                let mut long_alleles: Vec<Allele> = Vec::new();
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

#[derive(Clone)]
struct ClusterCenters {
    clusters: Vec<ClusterCenter>,
}

#[derive(Clone)]
struct ClusterCenter {
    center: Vec<f32>,
}

struct ThreadData {
    best_total_log_probability: HashMap<i32, f32>, // best log prob for each contig
    cluster_centers: HashMap<i32, ClusterCenters>,
    rng: StdRng,
    solves_per_thread: usize,
    thread_num: usize,
} 

impl ThreadData {
    fn from_seed(seed: [u8; 32], solves_per_thread: usize, thread_num: usize) -> ThreadData {
        ThreadData {
            best_total_log_probability: HashMap::new(),
            rng: SeedableRng::from_seed(seed),
            cluster_centers: HashMap::new(),
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
        let mut best_centers: HashMap<i32, ClusterCenters> = HashMap::new();
        for iteration in 0..thread_data.solves_per_thread {
            for (contig, loci) in contig_loci.loci.iter() {
                
                let cluster_centers: ClusterCenters = init_cluster_centers(loci.len(), params, &mut thread_data.rng);
                let contig_hic_links = hic_links.get(contig).unwrap();
                let long_contig_hic_links = long_hic_links.get(contig).unwrap();
                eprintln!("solve with LONG LINKS ONLY");
                let (_log_loss, cluster_centers, _hic_probabilities, _) =  // first solve with long links only
                    expectation_maximization(loci.len(), cluster_centers, &long_contig_hic_links, params, iteration, thread_data.thread_num);
                eprintln!("ALL HIC LINKS");
                let (log_loss, cluster_centers, hic_probabilities, hic_likelihoods) = 
                    expectation_maximization(loci.len(), cluster_centers, &contig_hic_links, params, iteration, thread_data.thread_num);
                
                let best_log_prob_so_far = thread_data.best_total_log_probability.entry(*contig).or_insert(f32::NEG_INFINITY);
                if &log_loss > best_log_prob_so_far {
                    thread_data.best_total_log_probability.insert(*contig, log_loss);
                    best_centers.insert(*contig, cluster_centers.clone());
                }
                eprintln!("thread {} contig {} iteration {} done with {}, best so far {}", 
                    thread_data.thread_num, contig, iteration, log_loss, thread_data.best_total_log_probability.get(contig).unwrap());
                    
                if READ_DEBUG {
                    for index in 0..cluster_centers.clusters[0].center.len() {
                        let mut count = 0;
                        if let Some(x) = locus_counts.get(&(*contig, index)){
                            count = *x;
                        }
                        let locus_hic_mols = contig_locus_hic_mols.get(contig).unwrap(); 
                        if let Some(hic_moldexes) = locus_hic_mols.get(&index) {
                            for hicdex in hic_moldexes {
                                let mut centers: Vec<(f32, f32)> = Vec::new();
                                for locus in contig_hic_links[*hicdex].loci.iter() {
                                    centers.push((cluster_centers.clusters[0].center[*locus], cluster_centers.clusters[1].center[*locus]));
                                }
                                let probs = &hic_probabilities[*hicdex];
                                let likes = &hic_likelihoods[*hicdex];
                                eprintln!("\thicread {} loci {:?} alleles {:?} clusters {:?} probs{:?}likes {:?}",
                                    hicdex, contig_hic_links[*hicdex].loci, contig_hic_links[*hicdex].alleles, centers, probs, likes);
                            }
                        }
                    }
                }   // END READ DEBUG
            } // end contig loop
        } // end cluster center solve iteration
        thread_data.cluster_centers = best_centers;
    }); // end parallel iter across threads

    let mut best_centers: HashMap<i32, ClusterCenters> = HashMap::new();
    let mut best_center_log_probabilities: HashMap<i32, f32> = HashMap::new();
    let mut best_center_threads: HashMap<i32, usize> = HashMap::new();
    for (thread_index, thread) in threads.iter().enumerate() {
        for (contig, log_probability) in thread.best_total_log_probability.iter() {
            let best = best_center_log_probabilities.entry(*contig).or_insert(f32::NEG_INFINITY);
            if log_probability > best {
                best_center_log_probabilities.insert(*contig, *log_probability);
                best_center_threads.insert(*contig, thread_index);
            }
        }
    }

    for (contig, thread_id) in best_center_threads.iter() {
        best_centers.insert(*contig, threads[*thread_id].cluster_centers.get(contig).unwrap().clone());
    }
    eprintln!("contig\thap1\thap2\treads\tassembly_allele\tphase");
    for contig in 1..(best_centers.len()+1) {
        let cluster_centers = best_centers.get(&(contig as i32)).unwrap();
        let loci = contig_loci.loci.get(&(contig as i32)).unwrap();
        for index in 0..cluster_centers.clusters[0].center.len() {
            let mut count = 0;
            if let Some(x) = locus_counts.get(&((contig as i32), index)) {
                count = *x;
            }
            let phase = match loci[index].allele {
                Allele::reference => cluster_centers.clusters[0].center[index] - cluster_centers.clusters[1].center[index],
                Allele::alternate => cluster_centers.clusters[1].center[index] - cluster_centers.clusters[0].center[index],
            };
            eprintln!("{}\t{}\t{}\t{}\t{:?}\t{}", contig, cluster_centers.clusters[0].center[index], 
                        cluster_centers.clusters[1].center[index], count, loci[index].allele, phase);
        }
    }
     

}

fn expectation_maximization(loci: usize, mut cluster_centers: ClusterCenters, hic_links: &Vec<HIC>, 
        params: &Params, epoch: usize, thread_num: usize) -> (f32, ClusterCenters,
        Vec<Vec<f32>>, Vec<f32>) { // this is a vector of the hic molecules probabilities to each cluster
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
    for _read in 0..hic_links.len() {
        final_log_probabilities.push(Vec::new());
    }

    let log_loss_change_limit = 0.000001 * (hic_links.len() as f32); // TODO fiddle with this in testing
    let mut last_log_loss = f32::NEG_INFINITY;
    let mut log_loss_change = 10000.0;
    let mut hic_probabilities: Vec<Vec<f32>> = Vec::new(); // TODO REMOVE DEBUG
    let mut hic_likelihoods: Vec<f32> = Vec::new();

    //while log_loss_change > log_loss_change_limit && iterations < 100 {
    while iterations < 150 { // TODO figure out something better here
        hic_probabilities.clear(); // TODO REMOVE DEBUG
        hic_likelihoods.clear();
        let mut log_likelihood = 0.0;
        reset_sums_denoms(loci, &mut sums, &mut denoms, params.ploidy);
        for (readdex, hic_read) in hic_links.iter().enumerate() {
            
            let log_likelihoods = bernoulli_likelihood(hic_read, &cluster_centers, log_prior);
            let read_likelihood = log_sum_exp(&log_likelihoods);
            hic_likelihoods.push(read_likelihood);
            log_likelihood += read_likelihood;
            let probabilities = normalize_in_log(&log_likelihoods);
            
            update_sums_denoms(&mut sums, &mut denoms, hic_read, &probabilities); // &probabilities); // TESTING
            hic_probabilities.push(probabilities);
            final_log_probabilities[readdex] = log_likelihoods;
        }
        total_log_loss = log_likelihood;
        log_loss_change = log_likelihood - last_log_loss;
        last_log_loss = log_likelihood;
       
        update_cluster_centers(loci, &sums, &denoms, &mut cluster_centers);
        iterations += 1;
        eprintln!("bernoulli\t{}\t{}\t{}\t{}\t{}", thread_num, epoch, iterations,  log_likelihood, log_loss_change);
    }
    (total_log_loss, cluster_centers, hic_probabilities, hic_likelihoods)
}


fn update_cluster_centers(loci: usize, sums: &Vec<Vec<f32>>, denoms: &Vec<Vec<f32>>, cluster_centers: &mut ClusterCenters) {
    for locus in 0..loci {
        for cluster in 0..sums.len() {
            if denoms[cluster][locus] > 1.0 {
                let update = sums[cluster][locus]/denoms[cluster][locus];
                cluster_centers.clusters[cluster].center[locus] = update.min(0.9999).max(0.0001);//max(0.0001, min(0.9999, update));
            }
        }
    }
}

fn update_sums_denoms(sums: &mut Vec<Vec<f32>>, denoms: &mut Vec<Vec<f32>>, hic_read: &HIC, probabilities: &Vec<f32>) {
    for locus in 0..hic_read.loci.len() {
        for (cluster, probability) in probabilities.iter().enumerate() {
            match hic_read.alleles[locus] {
                reference => sums[cluster][hic_read.loci[locus]] += 0.0,
                alternate => sums[cluster][hic_read.loci[locus]] += probability,// * hic_read.long_weighting,
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

fn bernoulli_likelihood(hic_read: &HIC, cluster_centers: &ClusterCenters, log_prior: f32) -> Vec<f32> {
    let mut log_probabilities: Vec<f32> = Vec::new();
    for (cluster, center) in cluster_centers.clusters.iter().enumerate() {
        log_probabilities.push(log_prior);
        for (locus_index, locus) in hic_read.loci.iter().enumerate() {
            log_probabilities[cluster] += ln_bernoulli(hic_read.alleles[locus_index], center.center[*locus]);
        }
    }
    
    log_probabilities
}

fn ln_bernoulli(allele: Allele, p: f32) -> f32 {
    match allele {
        reference => (1.0 - p).ln(),
        alternate => p.ln(),
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

fn init_cluster_centers(loci: usize, params: &Params, rng: &mut StdRng) -> ClusterCenters {
    let mut centers: ClusterCenters = ClusterCenters{ clusters: Vec::new() };
    for _cluster in 0..params.ploidy {
        let mut center = ClusterCenter{ center: Vec::new() };
        for _ in 0..loci {
            center.center.push(rng.gen::<f32>().min(0.99).max(0.01));
        }
        centers.clusters.push(center);
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

