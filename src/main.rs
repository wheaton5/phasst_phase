#[macro_use]
extern crate clap;
extern crate disjoint_set;
extern crate hashbrown;
extern crate phasst_lib;
extern crate rand;
extern crate rayon;
extern crate bio;

use bio::io::fasta;
use bio::utils::TextSlice;
use bio::io::fasta::Record;
use std::path::Path;

use phasst_lib::{
    load_assembly_kmers, load_hic, load_hifi, load_linked_read_barcodes, Assembly, HicMols,
    HifiMols, Kmers, LinkedReadBarcodes,
};
use rayon::prelude::*;

use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;

use std::fs::File;
use std::io::{BufWriter, Write};

use disjoint_set::DisjointSet;
use hashbrown::{HashMap, HashSet};

use clap::App;

const LONG_RANGE_HIC: usize = 15000;
const LONG_RANGE_HIC_WEIGHTING: f32 = 100.0;
const MIN_ALLELE_FRACTION_HIC: f32 = 0.15;
const MOLECULE_DEBUG: bool = false;
const LINKED_READ_NEW_MOLECULE_GAP: usize = 30000;
const LINKED_READ_MIN_ALLELES: usize = 15;

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
    reference: i32,
    alternate: i32,
}

#[derive(Debug, Clone, Copy)]
enum Allele {
    Ref,
    Alt,
}

#[derive(Clone)]
struct Molecule {
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
    let hic_mols = load_hic(Some(&params.hic_mols), &kmers);
    eprintln!("loading long reads");
    let ccs = load_hifi(Some(&params.ccs_mols), &kmers);
    eprintln!("loading linked reads");
    let txg_barcodes = load_linked_read_barcodes(Some(&params.txg_mols), &kmers);
    eprintln!("loading assembly kmers");
    let assembly = load_assembly_kmers(&params.assembly_kmers, &params.assembly_fasta, &kmers);

    eprintln!("finding good loci");
    let allele_fractions = get_allele_fractions(&hic_mols); // MAYBE ADD LINKED Molecule AND LONG Molecule to this?
    let bad_alleles = get_bad_alleles(&hic_mols);
    let variant_contig_order: ContigLoci =
        good_assembly_loci(&assembly, &allele_fractions, &bad_alleles);
    eprintln!("finding good hic reads");
    let (hic_links, long_hic_links) = gather_hic_links(&hic_mols, &variant_contig_order);
    let ccs_mols = gather_hifi_mols(&ccs, &variant_contig_order);
    let txg_mols = gather_txg_mols(&txg_barcodes, &variant_contig_order);
    eprintln!("phasing");
    //MAYBE LATER let mut connected_components = get_connected_components(&hic_links, &variant_contig_order, params.min_hic_links);
    let (phasing, best_centers) = phasst_phase_main(
        &params,
        &hic_links,
        &long_hic_links,
        &ccs_mols,
        &txg_mols,
        &variant_contig_order,
        &assembly,
        &kmers,
    );
    let contig_chunks = assess_breakpoints(
        &hic_links,
        &ccs_mols,
        &txg_mols,
        &params,
        &variant_contig_order,
        &phasing,
        &assembly
    );
    output_phased_vcf(&kmers, &params, &best_centers, &variant_contig_order, &assembly, &contig_chunks);

}

fn output_phased_vcf(kmers: &Kmers, params: &Params, best_centers: &HashMap<i32, ClusterCenters>, contig_loci: &ContigLoci, assembly: &Assembly, contig_chunks: &HashMap<i32, Vec<(usize, usize)>>) {
    let mut output = params.output.to_string();
    output.push_str("/phasing_breaks.vcf");
    let f = File::create(output).expect("Unable to create file");
    let mut f = BufWriter::new(f);
    let mut phasing: HashMap<i32, Vec<Option<bool>>> = HashMap::new();
    //for (contig, center) in best_centers.iter() {
    for contig in 1..assembly.contig_names.len() {
        let contig = &(contig as i32);
        let empty: ClusterCenters = ClusterCenters{ clusters: Vec::new() };
        let center = best_centers.get(contig).unwrap_or(&empty);
        let contig_phasing = phasing.entry(*contig as i32).or_insert(Vec::new());
        let loci = contig_loci.loci.get(contig).unwrap();
        let contig_name = &assembly.contig_names[*contig as usize];

        let chunks = contig_chunks.get(contig).expect("why do you hate me");
        for (chunkdex, (left, right)) in chunks.iter().enumerate() {
            for ldex in *left..*right { //0..center.clusters[0].center.len() {
                // output is semi-vcf contig\tpos\t.\tREF\tALT\tqual\tfilter\tinfo\tformat\tsample
                let chunk_name = vec![contig_name.to_string(), (chunkdex+1).to_string(), left.to_string(), right.to_string()].join("_");
                let locus = loci[ldex];
                let pos = locus.position;
                let reference;
                let alternate;
                let flip;
                match locus.allele {
                    Allele::Ref => {
                        reference = kmers.kmers.get(&locus.reference).unwrap().to_string();
                        alternate = kmers.kmers.get(&locus.alternate).unwrap().to_string();
                        flip = false;
                    }
                    Allele::Alt => {
                        reference = kmers.kmers.get(&locus.alternate).unwrap().to_string();
                        alternate = kmers.kmers.get(&locus.reference).unwrap().to_string();
                        flip = true;
                    }
                }

                let mut genotype: Vec<String> = Vec::new();

                for cluster in center.clusters.iter() {
                    let value = cluster.center[ldex];
                    if value > 0.99 {
                        if !flip {
                            genotype.push("1".to_string());
                        } else {
                            genotype.push("0".to_string())
                        }
                    } else if value < 0.01 {
                        if !flip {
                            genotype.push("0".to_string())
                        } else {
                            genotype.push("1".to_string());
                        }
                    } else {
                        genotype.push(".".to_string());
                    }
                }
                if genotype[0] == "0" && genotype[1] == "1" {
                    if !flip {
                        contig_phasing.push(Some(true));
                    } else {
                        contig_phasing.push(Some(false));
                    }
                } else if genotype[0] == "1" && genotype[1] == "0" {
                    if !flip {
                        contig_phasing.push(Some(false));
                    } else {
                        contig_phasing.push(Some(true));
                    }
                } else {
                    contig_phasing.push(None);
                }

                let genotype = vec![genotype.join("|"), "60".to_string()].join(":");
                let mut line_vec: Vec<String> = vec![
                    chunk_name.to_string(),
                    pos.to_string(),
                    ".".to_string(),
                    reference,
                    alternate,
                    ".".to_string(),
                    ".".to_string(),
                    ".".to_string(),
                    "GT:PQ".to_string(),
                    genotype,
                ];
                let mut line = line_vec.join("\t");
                line.push_str("\n");
                f.write_all(line.as_bytes()).expect("Unable to write data");
            }
        }
        
    }
}

fn assess_breakpoints(
    hic_links: &HashMap<i32, Vec<Molecule>>,
    ccs_mols: &HashMap<i32, Vec<Molecule>>,
    txg_mols: &HashMap<i32, Vec<Molecule>>,
    params: &Params,
    contig_loci: &ContigLoci,
    phasing: &Phasing,
    assembly: &Assembly
) -> HashMap<i32, Vec<(usize, usize)>> {

    let mut chunks: HashMap<i32, Vec<(usize, usize)>> = HashMap::new(); // ranges for each contig

    //for (contig, hic) in hic_links.iter() {
    for contig in 1..assembly.contig_names.len() {
        let contig_name = &assembly.contig_names[contig];
        let contig = &(contig as i32);
        let empty: Vec<Molecule> = Vec::new();
        let hic = hic_links.get(contig).unwrap_or(&empty);//(&format!("contig {} {} has no hic links, contig size {}", contig, contig_name, assembly.contig_sizes.get(contig).unwrap()));
        let mut in_chunk = true;
        let contig_chunk = chunks.entry(*contig).or_insert(Vec::new());
        let mut current_chunk = (0,0);
        let ccs = ccs_mols
            .get(contig)
            .unwrap_or(&empty);//expect("cant find contig for hifi mols");
        let txg = txg_mols.get(contig).unwrap_or(&empty);//.expect("cant find contig for txg mols");
        let empty2: Vec<Option<bool>> = Vec::new();
        let contig_phasing = phasing.phasing.get(contig).unwrap_or(&empty2);//expect("contig not in phasing?");

        let mut locus_hic_mols: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut locus_ccs_mols: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut locus_txg_mols: HashMap<usize, Vec<usize>> = HashMap::new();
        for (read_index, hic_read) in hic.iter().enumerate() {
            for locus in hic_read.loci.iter() {
                let blah = locus_hic_mols.entry(*locus).or_insert(Vec::new());
                blah.push(read_index);
            }
        }
        for (read_index, ccs_read) in ccs.iter().enumerate() {
            for locus in ccs_read.loci.iter() {
                let blah = locus_ccs_mols.entry(*locus).or_insert(Vec::new());
                blah.push(read_index);
            }
        }
        for (read_index, txg_read) in txg.iter().enumerate() {
            for locus in txg_read.loci.iter() {
                let blah = locus_txg_mols.entry(*locus).or_insert(Vec::new());
                blah.push(read_index);
            }
        }
        let empty3: Vec<ContigLocus> = Vec::new();
        let loci = contig_loci
            .loci
            .get(contig)
            .unwrap_or(&empty3);

        let mut current_hic_mol_set: HashSet<usize> = HashSet::new();
        let mut current_txg_mol_set: HashSet<usize> = HashSet::new();
        let mut current_ccs_mol_set: HashSet<usize> = HashSet::new();
        for index in 0..params.break_window {
            let left = 0;
            let right = params.break_window;
            if let Some(hic_indexes) = locus_hic_mols.get(&index) {
                for hic_index in hic_indexes {
                    check_add(
                        &hic[*hic_index],
                        *hic_index,
                        &mut current_hic_mol_set,
                        left,
                        right,
                    );
                }
            }
            if let Some(ccs_indexes) = locus_ccs_mols.get(&index) {
                for ccs_index in ccs_indexes {
                    check_add(
                        &ccs[*ccs_index],
                        *ccs_index,
                        &mut current_ccs_mol_set,
                        left,
                        right,
                    );
                }
            }
            if let Some(txg_indexes) = locus_txg_mols.get(&index) {
                for txg_index in txg_indexes {
                    check_add(
                        &txg[*txg_index],
                        *txg_index,
                        &mut current_txg_mol_set,
                        left,
                        right,
                    );
                }
            }
        }
        for (mid, locus) in loci.iter().enumerate() {
            let mut left = mid - params.break_window;
            let right = mid + params.break_window;
            if mid > params.break_window {
                if let Some(hic_indexes) = locus_hic_mols.get(&(left - 1)) {
                    for hic_index in hic_indexes {
                        check_remove(
                            &hic[*hic_index],
                            *hic_index,
                            &mut current_hic_mol_set,
                            left,
                            right,
                        );
                    }
                }
                if let Some(ccs_indexes) = locus_ccs_mols.get(&(left - 1)) {
                    for ccs_index in ccs_indexes {
                        check_remove(
                            &ccs[*ccs_index],
                            *ccs_index,
                            &mut current_ccs_mol_set,
                            left,
                            right,
                        );
                    }
                }
                if let Some(txg_indexes) = locus_txg_mols.get(&(left - 1)) {
                    for txg_index in txg_indexes {
                        check_remove(
                            &txg[*txg_index],
                            *txg_index,
                            &mut current_txg_mol_set,
                            left,
                            right,
                        );
                    }
                }
            } else {
                left = 0;
            }
            if let Some(hic_indexes) = locus_hic_mols.get(&(right - 1)) {
                for hic_index in hic_indexes {
                    check_add(
                        &hic[*hic_index],
                        *hic_index,
                        &mut current_hic_mol_set,
                        left,
                        right,
                    );
                }
            }
            if let Some(ccs_indexes) = locus_ccs_mols.get(&(right - 1)) {
                for ccs_index in ccs_indexes {
                    check_add(
                        &ccs[*ccs_index],
                        *ccs_index,
                        &mut current_ccs_mol_set,
                        left,
                        right,
                    );
                }
            }
            if let Some(txg_indexes) = locus_txg_mols.get(&(right - 1)) {
                for txg_index in txg_indexes {
                    check_add(
                        &txg[*txg_index],
                        *txg_index,
                        &mut current_txg_mol_set,
                        left,
                        right,
                    );
                }
            }

            let mut counts: [u16; 4] = [0; 4];

            for hic_moldex in current_hic_mol_set.iter() {
                let hicmol = &hic[*hic_moldex];
                let mut molcounts: [u16; 4] = [0; 4];
                for index1 in 0..hicmol.loci.len() {
                    for index2 in (index1 + 1)..hicmol.loci.len() {
                        let locus1 = hicmol.loci[index1];
                        let locus2 = hicmol.loci[index2];
                        if locus1 < mid && locus1 >= left && locus2 >= mid && locus2 < right {
                            if let Some(phase_left) = contig_phasing[locus1] {
                                let mut phase_left = phase_left;
                                if let Some(phase_right) = contig_phasing[locus2] {
                                    let mut phase_right = phase_right;
                                    match hicmol.alleles[index1] {
                                        Allele::Alt => phase_left = !phase_left,
                                        Allele::Ref => (),
                                    }
                                    match hicmol.alleles[index2] {
                                        Allele::Alt => phase_right = !phase_right,
                                        Allele::Ref => (),
                                    }
                                    if phase_left && phase_right {
                                        molcounts[0] += 1;
                                    } else if !phase_left && !phase_right {
                                        molcounts[1] += 1;
                                    } else if phase_left && !phase_right {
                                        molcounts[2] += 1;
                                    } else if !phase_left && phase_right {
                                        molcounts[3] += 1;
                                    }
                                }
                            }
                        }
                    }
                }
                let mut best = 0;
                let mut bestdex = 0;
                let mut secondbest = 0;
                for index in 0..4 {
                    if molcounts[index] > best {
                        secondbest = best;
                        best = molcounts[index];
                        bestdex = index;
                    }
                }
                if best > 0 && secondbest == 0 {
                    counts[bestdex] += 1;
                }
            }
            eprintln!(
                "contig {}, mid {} position {}, break_counts {:?}",
                contig, locus.position, mid, counts
            );
            let cis = (counts[0]+counts[1]) as f64;
            let total = ((counts[2]+counts[3]) as f64) + cis;
            if mid > 250 && mid < loci.len() - 250 {
                if cis/(total + 1.0) > 0.75 && in_chunk {
                    current_chunk = (current_chunk.0, locus.position);
                } else if cis/(total + 1.0) < 0.75 && in_chunk {
                    in_chunk = false;
                    if current_chunk.1 > current_chunk.0 {
                        contig_chunk.push(current_chunk);
                        eprintln!("adding chunk for contig {}, chunk {:?}", contig_name, current_chunk);
                    }
                    current_chunk = (locus.position+1, locus.position+1);
                } else if cis/(total + 1.0) > 0.75 && !in_chunk {
                    in_chunk = true;
                    //current_chunk = (locus.position, locus.position);
                }
            } else if in_chunk {
                current_chunk = (current_chunk.0, locus.position);
            }
        }
        if in_chunk || contig_chunk.len() == 0 {
            current_chunk = (current_chunk.0, *assembly.contig_sizes.get(contig).unwrap());
            contig_chunk.push(current_chunk);
            eprintln!("adding chunk at finish for contig {}, chunk {:?}", contig_name, current_chunk);
        }
        if contig_chunk.len() > 1 {
            eprintln!("contig {} with size {} is split into {} chunks", contig, *assembly.contig_sizes.get(contig).unwrap(), contig_chunk.len());
            for (start, end) in contig_chunk.iter() {
                eprintln!("\t{} - {}", start, end);
            }
        } else {
            eprintln!("contig {} has no breaks", contig);
        }
        

    }
    let reader =  fasta::Reader::from_file(Path::new(&params.assembly_fasta.to_string())).expect("fasta not found");
    let mut writer = fasta::Writer::to_file(Path::new(&format!("{}/breaks.fa",params.output))).expect("cannot open fasta writer");
    for record in reader.records() {
        let record = record.unwrap();
        let contig_name = record.id().to_string();
        let contig_id = assembly.contig_ids.get(&contig_name).unwrap();
        
        if !chunks.contains_key(contig_id) {
            eprintln!("contig has no chunks??? {}", contig_id);
            let size = assembly.contig_sizes.get(contig_id).unwrap();
            let range = chunks.entry(*contig_id).or_insert(Vec::new());
            range.push((0,*size));
        } 
        let ranges = chunks.get(contig_id).unwrap();
        
        for (index, (start, stop)) in ranges.iter().enumerate() {
            let mut new_contig_name = contig_name.to_string();
            if ranges.len() > 0 {
                let list = vec![new_contig_name, (index+1).to_string(), start.to_string(), stop.to_string()];
                new_contig_name = list.join("_");
            }
            let seq: TextSlice = &record.seq()[*start..*stop];
            let record = Record::with_attrs(&new_contig_name, None, &seq);
            writer.write_record(&record).expect("could not write record");
        }
    }
    chunks
}


fn check_remove(
    mol: &Molecule,
    mol_index: usize,
    current_mol_set: &mut HashSet<usize>,
    left: usize,
    right: usize,
) {
    let mut count = 0;
    for locus in mol.loci.iter() {
        if *locus >= left && *locus < right {
            count += 1;
        }
    }
    if count < 2 {
        current_mol_set.remove(&mol_index);
    }
}

fn check_add(
    mol: &Molecule,
    mol_index: usize,
    current_mol_set: &mut HashSet<usize>,
    left: usize,
    right: usize,
) {
    let mut count = 0;
    for locus in mol.loci.iter() {
        if *locus >= left && *locus < right {
            count += 1;
        }
    }
    if count < 2 {
        current_mol_set.insert(mol_index);
    }
}

fn allele(kmer: i32) -> Allele {
    match kmer.abs() % 2 == 0 {
        true => Allele::Ref,
        false => Allele::Alt,
    }
}

fn get_bad_alleles(hic_mols: &HicMols) -> HashSet<i32> {
    let mut bad: HashSet<i32> = HashSet::new();
    for mol in hic_mols.get_hic_molecules() {
        for i in 0..mol.len() {
            for j in (i + 1)..mol.len() {
                if Kmers::pair(mol[i].abs()) == mol[j].abs() {
                    bad.insert(mol[i].abs());
                    bad.insert(mol[j].abs());
                }
            }
        }
    }
    bad
}

fn get_connected_components(
    hic_links: &HashMap<i32, Vec<Molecule>>,
    variant_contig_order: &ContigLoci,
    min_links: u32,
) -> HashMap<i32, DisjointSet<usize>> {
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
                for j in (i + 1)..hic_mol.loci.len() {
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
            let count = allele_counts.entry(canonical).or_insert([0; 2]);
            if var.abs() % 2 == 0 {
                count[0] += 1;
            } else {
                count[1] += 1;
            }
        }
    }
    for (canonical, counts) in allele_counts {
        let fraction = (counts[0].min(counts[1]) as f32) / ((counts[0] + counts[1]) as f32);
        allele_fractions.insert(canonical, fraction);
        allele_fractions.insert(Kmers::pair(canonical), fraction);
    }
    allele_fractions
}

fn good_assembly_loci(
    assembly: &Assembly,
    allele_fractions: &HashMap<i32, f32>,
    bad_alleles: &HashSet<i32>,
) -> ContigLoci {
    // returning a map from kmer id to contig id and position
    let mut variant_contig_order: HashMap<i32, ContigLocus> = HashMap::new();

    let mut contig_positions: HashMap<i32, Vec<(usize, i32, i32)>> = HashMap::new();
    for (kmer, (contig, num, _order, position)) in assembly.variants.iter() {
        // TODODODODODODODODODODo
        if *contig > 5 {
           continue;
        } // TODO remove

        if assembly.variants.contains_key(&Kmers::pair(*kmer)) {
            continue;
        } // we see both ref and alt in assembly, skip
        if let Some(fraction) = allele_fractions.get(&kmer.abs()) {
            if *fraction < MIN_ALLELE_FRACTION_HIC {
                continue;
            }
            if bad_alleles.contains(kmer) {
                continue;
            }
        } else {
            continue;
        }

        if *num > 1 {
            continue;
        } // we see this kmer multiple times in the assembly, skip
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
            let reference;
            let alternate;
            match allele(*kmer) {
                Allele::Ref => {
                    reference = kmer.abs();
                    alternate = reference - 1;
                }
                Allele::Alt => {
                    alternate = kmer.abs();
                    reference = alternate + 1;
                }
            }
            let locus = ContigLocus {
                contig_id: *contig,
                index: index,
                position: *position,
                allele: allele(*kmer),
                reference: reference,
                alternate: alternate,
            };
            variant_contig_order.insert(*kmer, locus);
            variant_contig_order.insert(Kmers::pair(*kmer), locus);
            contig_loci.push(locus);
        }
    }

    ContigLoci {
        kmers: variant_contig_order,
        loci: loci,
    }
}

fn gather_txg_mols(
    txg_barcodes: &LinkedReadBarcodes,
    variant_contig_order: &ContigLoci,
) -> HashMap<i32, Vec<Molecule>> {
    let mut to_return: HashMap<i32, Vec<Molecule>> = HashMap::new();
    let mut number = 0;
    for (contig, _) in variant_contig_order.loci.iter() {
        to_return.insert(*contig, Vec::new());
    }
    for (_bc, vars) in txg_barcodes.get_linked_read_barcodes().enumerate() {
        let mut contig_vars: HashMap<i32, Vec<(usize, usize, i32)>> = HashMap::new(); // map from contig to position and index and kmer
        for var in vars {
            if let Some(locus) = variant_contig_order.kmers.get(&var.abs()) {
                let contig = contig_vars.entry(locus.contig_id).or_insert(Vec::new());
                contig.push((locus.position, locus.index, *var));
            }
        }
        for (contig_id, vars) in contig_vars.iter_mut() {
            if vars.len() < LINKED_READ_MIN_ALLELES {
                continue;
            }
            vars.sort();

            let mut last_position = 0;
            let mut growing_molecule = Molecule {
                loci: Vec::new(),
                alleles: Vec::new(),
                long_weighting: 1.0,
            };
            let contig = to_return.entry(*contig_id).or_insert(Vec::new());
            for (position, index, var) in vars {
                if *position - last_position > LINKED_READ_NEW_MOLECULE_GAP {
                    if growing_molecule.loci.len() > LINKED_READ_MIN_ALLELES {
                        contig.push(growing_molecule);
                        number += 1;
                    }
                    growing_molecule = Molecule {
                        loci: Vec::new(),
                        alleles: Vec::new(),
                        long_weighting: 1.0,
                    };
                }
                growing_molecule.loci.push(*index);
                growing_molecule.alleles.push(allele(*var));
                last_position = *position;
            }
            if growing_molecule.loci.len() > LINKED_READ_MIN_ALLELES {
                contig.push(growing_molecule);
                number += 1;
            }
        }
    }
    eprintln!("linked read molecules detected {}", number);
    to_return
}

fn gather_hifi_mols(
    hifi: &HifiMols,
    variant_contig_order: &ContigLoci,
) -> HashMap<i32, Vec<Molecule>> {
    let mut to_return: HashMap<i32, Vec<Molecule>> = HashMap::new();
    for (contig, _) in variant_contig_order.loci.iter() {
        to_return.insert(*contig, Vec::new());
    }

    let mut total = 0;
    for mol in hifi.get_hifi_molecules() {
        total += 1;
        let mut the_contig: Option<i32> = None;
        let mut loci: Vec<usize> = Vec::new();
        let mut alleles: Vec<Allele> = Vec::new();
        for var in mol {
            if let Some(locus) = variant_contig_order.kmers.get(&var.abs()) {
                if let Some(chrom) = the_contig {
                    if locus.contig_id == chrom {
                        loci.push(locus.index);
                        alleles.push(allele(*var));
                    } else {
                        //do something here
                    }
                } else {
                    loci.push(locus.index);
                    alleles.push(allele(*var));
                    the_contig = Some(locus.contig_id);
                }
            }
        }
        if loci.len() > 1 {
            let contig = to_return.entry(the_contig.unwrap()).or_insert(Vec::new());
            contig.push(Molecule {
                alleles: alleles,
                loci: loci,
                long_weighting: 1.0,
            });
        }
    }

    eprintln!("hifi reads {}", total);
    to_return
}

fn gather_hic_links(
    hic_molecules: &HicMols,
    variant_contig_order: &ContigLoci,
) -> (HashMap<i32, Vec<Molecule>>, HashMap<i32, Vec<Molecule>>) {
    // returns map from contig id to list of HIC data structures
    let mut hic_mols: HashMap<i32, Vec<Molecule>> = HashMap::new();
    let mut long_hic_mols: HashMap<i32, Vec<Molecule>> = HashMap::new();
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
            if used.contains(&var.abs()) {
                used_count += 1;
                continue;
            }
            if let Some(locus) = variant_contig_order.kmers.get(&var.abs()) {
                if let Some(chrom) = the_contig {
                    min = min.min(locus.position);
                    max = max.max(locus.position);
                    if locus.contig_id == chrom {
                        loci.push(locus.index);
                        alleles.push(allele(*var));
                    } else {
                        diff_contig += 1;
                    }
                } else {
                    min = min.min(locus.position);
                    max = max.max(locus.position);
                    the_contig = Some(locus.contig_id);
                    loci.push(locus.index);
                    alleles.push(allele(*var));
                }
            } else {
                not_assembly += 1;
            }
            used.insert(var.abs());
        }
        if loci.len() > 1 {
            let contig_mols = hic_mols.entry(the_contig.unwrap()).or_insert(Vec::new());
            let mut long_range = 1.0;
            if max - min > LONG_RANGE_HIC {
                long_range = LONG_RANGE_HIC_WEIGHTING;
                let long_contig_mols = long_hic_mols
                    .entry(the_contig.unwrap())
                    .or_insert(Vec::new());
                let mut long_loci: Vec<usize> = Vec::new();
                let mut long_alleles: Vec<Allele> = Vec::new();
                for index in 0..loci.len() {
                    long_loci.push(loci[index]);
                    long_alleles.push(alleles[index]);
                }
                long_contig_mols.push(Molecule {
                    loci: long_loci,
                    alleles: long_alleles,
                    long_weighting: LONG_RANGE_HIC_WEIGHTING,
                });
                total_long_links += 1;
            }
            contig_mols.push(Molecule {
                loci: loci,
                alleles: alleles,
                long_weighting: long_range,
            });
            total += 1;
        }
    }
    eprintln!(
        "after culling we have {} hic molecules hitting >=2 distinct loci",
        total
    );
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

struct Phasing {
    phasing: HashMap<i32, Vec<Option<bool>>>,
}

fn phasst_phase_main(
    params: &Params,
    hic_links: &HashMap<i32, Vec<Molecule>>,
    long_hic_links: &HashMap<i32, Vec<Molecule>>,
    hifi_mols: &HashMap<i32, Vec<Molecule>>,
    linked_read_mols: &HashMap<i32, Vec<Molecule>>,
    contig_loci: &ContigLoci,
    assembly: &Assembly,
    kmers: &Kmers,
) -> (Phasing, HashMap<i32, ClusterCenters>) {
    let seed = [params.seed; 32];
    let mut rng: StdRng = SeedableRng::from_seed(seed);
    let mut threads: Vec<ThreadData> = Vec::new();
    let solves_per_thread = ((params.restarts as f32) / (params.threads as f32)).ceil() as usize;
    for i in 0..params.threads {
        threads.push(ThreadData::from_seed(
            new_seed(&mut rng),
            solves_per_thread,
            i,
        ));
    }

    // First lets keep track of how many times each locus is hit by a hic molecule (that also hit something else)
    let mut locus_counts: HashMap<(i32, usize), u32> = HashMap::new();
    let mut contig_locus_hic_mols: HashMap<i32, HashMap<usize, Vec<usize>>> = HashMap::new(); // TODO REMOVE DEBUG map from contig to a map from locus to vec of hic mols
    for (contig, _loci) in contig_loci.loci.iter() {
        let contig_hic_links = hic_links.get(contig).unwrap();
        let locus_hic = contig_locus_hic_mols
            .entry(*contig)
            .or_insert(HashMap::new()); //TODO REMOVE DEBUG
        for (read_index, hic_read) in contig_hic_links.iter().enumerate() {
            for locus in hic_read.loci.iter() {
                let hic_mols = locus_hic.entry(*locus).or_insert(Vec::new()); // TODO REMOVE DEBUG
                hic_mols.push(read_index); //TODO REMOVE DEBUG
                let count = locus_counts.entry((*contig, *locus)).or_insert(0);
                *count += 1;
            }
        }
    }
    let empty_vec: Vec<Molecule> = Vec::new();
    threads.par_iter_mut().for_each(|thread_data| {
        let mut best_centers: HashMap<i32, ClusterCenters> = HashMap::new();
        for iteration in 0..thread_data.solves_per_thread {
            for (contig, loci) in contig_loci.loci.iter() {

                let cluster_centers: ClusterCenters = init_cluster_centers(loci.len(), params, &mut thread_data.rng);
                let contig_hic_links = hic_links.get(contig).unwrap();
                let contig_hifi_mols = hifi_mols.get(contig).unwrap();
                let contig_txg_mols = linked_read_mols.get(contig).unwrap();
                let long_contig_hic_links = long_hic_links.get(contig).unwrap();
                eprintln!("solve with LONG LINKS ONLY");
                let (_log_loss, cluster_centers, _hic_probabilities, _) =  // first solve with long links only
                    expectation_maximization(loci.len(), cluster_centers, &long_contig_hic_links, &empty_vec, &empty_vec, params, iteration, thread_data.thread_num);
                eprintln!("ALL MoleculeS: {} hic and {} hifi", contig_hic_links.len(), contig_hifi_mols.len());
                let (log_loss, cluster_centers, hic_probabilities, hic_likelihoods) =
                    expectation_maximization(loci.len(), cluster_centers, &contig_hic_links, &contig_hifi_mols, &contig_txg_mols, params, iteration, thread_data.thread_num);
                let best_log_prob_so_far = thread_data.best_total_log_probability.entry(*contig).or_insert(f32::NEG_INFINITY);
                if &log_loss > best_log_prob_so_far {
                    thread_data.best_total_log_probability.insert(*contig, log_loss);
                    best_centers.insert(*contig, cluster_centers.clone());
                }
                eprintln!("thread {} contig {} iteration {} done with {}, best so far {}", 
                    thread_data.thread_num, contig, iteration, log_loss, thread_data.best_total_log_probability.get(contig).unwrap());

                if MOLECULE_DEBUG {
                    for index in 0..cluster_centers.clusters[0].center.len() {
                        let mut count = 0;
                        if let Some(x) = locus_counts.get(&(*contig, index)){
                            count = *x;
                        }
                        let locus_hic_mols = contig_locus_hic_mols.get(contig).unwrap();

                        let phase = match loci[index].allele {
                            Allele::Ref => cluster_centers.clusters[0].center[index] - cluster_centers.clusters[1].center[index],
                            Allele::Alt => cluster_centers.clusters[1].center[index] - cluster_centers.clusters[0].center[index],
                        };
                        eprintln!("{}\t{}\t{}\t{}\t{:?}\t{}", contig, cluster_centers.clusters[0].center[index], 
                        cluster_centers.clusters[1].center[index], count, loci[index].allele, phase);
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
                } // END Molecule DEBUG
            } // end contig loop
        } // end cluster center solve iteration
        thread_data.cluster_centers = best_centers;
    }); // end parallel iter across threads

    let mut best_centers: HashMap<i32, ClusterCenters> = HashMap::new();
    let mut best_center_log_probabilities: HashMap<i32, f32> = HashMap::new();
    let mut best_center_threads: HashMap<i32, usize> = HashMap::new();
    for (thread_index, thread) in threads.iter().enumerate() {
        for (contig, log_probability) in thread.best_total_log_probability.iter() {
            let best = best_center_log_probabilities
                .entry(*contig)
                .or_insert(f32::NEG_INFINITY);
            if log_probability > best {
                best_center_log_probabilities.insert(*contig, *log_probability);
                best_center_threads.insert(*contig, thread_index);
            }
        }
    }

    for (contig, thread_id) in best_center_threads.iter() {
        best_centers.insert(
            *contig,
            threads[*thread_id]
                .cluster_centers
                .get(contig)
                .unwrap()
                .clone(),
        );
    }
    println!("contig\tposition\thap1\thap2\treads\tassembly_allele\tphase");

    for contig in 1..(best_centers.len() + 2) {
        if !contig_loci.loci.contains_key(&(contig as i32)) {
            eprintln!(
                "contig loci doesnt contain contig {}, {}",
                contig, assembly.contig_names[contig]
            );
            continue;
        }
        if !best_centers.contains_key(&(contig as i32)) {
            eprintln!(
                "best centers doesnt contain contig {}, {}",
                contig, assembly.contig_names[contig]
            );
            continue;
        }
        let cluster_centers = best_centers.get(&(contig as i32)).unwrap();
        let loci = contig_loci.loci.get(&(contig as i32)).unwrap();
        for index in 0..cluster_centers.clusters[0].center.len() {
            let mut count = 0;
            if let Some(x) = locus_counts.get(&((contig as i32), index)) {
                count = *x;
            }
            let position = loci[index].position;
            let phase = match loci[index].allele {
                Allele::Ref => {
                    cluster_centers.clusters[0].center[index]
                        - cluster_centers.clusters[1].center[index]
                }
                Allele::Alt => {
                    cluster_centers.clusters[1].center[index]
                        - cluster_centers.clusters[0].center[index]
                }
            };
            let alt_frac = match loci[index].allele {
                Allele::Ref => cluster_centers.clusters[1].center[index],
                Allele::Alt => cluster_centers.clusters[0].center[index],
            };
            let ref_frac = match loci[index].allele {
                Allele::Ref => cluster_centers.clusters[0].center[index],
                Allele::Alt => cluster_centers.clusters[1].center[index],
            };
            println!(
                "{}\t{}\t{}\t{}\t{}\t{:?}\t{}",
                contig, position, ref_frac, alt_frac, count, loci[index].allele, phase
            );
        }
    }
    let mut output = params.output.to_string();
    output.push_str("/phasing.vcf");
    let f = File::create(output).expect("Unable to create file");
    let mut f = BufWriter::new(f);
    let mut phasing: HashMap<i32, Vec<Option<bool>>> = HashMap::new();
    for (contig, center) in best_centers.iter() {
        let contig_phasing = phasing.entry(*contig as i32).or_insert(Vec::new());
        let loci = contig_loci.loci.get(contig).unwrap();
        let contig_name = &assembly.contig_names[*contig as usize];

        for ldex in 0..center.clusters[0].center.len() {
            // output is semi-vcf contig\tpos\t.\tREF\tALT\tqual\tfilter\tinfo\tformat\tsample
            let locus = loci[ldex];
            let pos = locus.position;
            let reference;
            let alternate;
            let flip;
            match locus.allele {
                Allele::Ref => {
                    reference = kmers.kmers.get(&locus.reference).unwrap().to_string();
                    alternate = kmers.kmers.get(&locus.alternate).unwrap().to_string();
                    flip = false;
                }
                Allele::Alt => {
                    reference = kmers.kmers.get(&locus.alternate).unwrap().to_string();
                    alternate = kmers.kmers.get(&locus.reference).unwrap().to_string();
                    flip = true;
                }
            }

            let mut genotype: Vec<String> = Vec::new();

            for cluster in center.clusters.iter() {
                let value = cluster.center[ldex];
                if value > 0.99 {
                    if !flip {
                        genotype.push("1".to_string());
                    } else {
                        genotype.push("0".to_string())
                    }
                } else if value < 0.01 {
                    if !flip {
                        genotype.push("0".to_string())
                    } else {
                        genotype.push("1".to_string());
                    }
                } else {
                    genotype.push(".".to_string());
                }
            }
            if genotype[0] == "0" && genotype[1] == "1" {
                if !flip {
                    contig_phasing.push(Some(true));
                } else {
                    contig_phasing.push(Some(false));
                }
            } else if genotype[0] == "1" && genotype[1] == "0" {
                if !flip {
                    contig_phasing.push(Some(false));
                } else {
                    contig_phasing.push(Some(true));
                }
            } else {
                contig_phasing.push(None);
            }

            let genotype = vec![genotype.join("|"), "60".to_string()].join(":");
            let mut line_vec: Vec<String> = vec![
                contig_name.to_string(),
                pos.to_string(),
                ".".to_string(),
                reference,
                alternate,
                ".".to_string(),
                ".".to_string(),
                ".".to_string(),
                "GT:PQ".to_string(),
                genotype,
            ];
            let mut line = line_vec.join("\t");
            line.push_str("\n");
            f.write_all(line.as_bytes()).expect("Unable to write data");
        }
    }
    (Phasing { phasing: phasing }, best_centers)
}

fn expectation_maximization(
    loci: usize,
    mut cluster_centers: ClusterCenters,
    hic_links: &Vec<Molecule>,
    hifi: &Vec<Molecule>,
    txg: &Vec<Molecule>,
    params: &Params,
    epoch: usize,
    thread_num: usize,
) -> (f32, ClusterCenters, Vec<Vec<f32>>, Vec<f32>) {
    // this is a vector of the hic molecules probabilities to each cluster
    if hic_links.len() == 0 {
        eprintln!("no hic links?");
    }
    let mut sums: Vec<Vec<f32>> = Vec::new();
    let mut denoms: Vec<Vec<f32>> = Vec::new();

    let mut variant_hic_reads: Vec<Vec<Molecule>> = Vec::new();
    for _ in 0..cluster_centers.clusters[0].center.len() {
        variant_hic_reads.push(Vec::new());
    }

    for hic in hic_links {
        for locus in hic.loci.iter() {
            variant_hic_reads[*locus].push(hic.clone());
        }
    }

    for cluster in 0..params.ploidy {
        sums.push(Vec::new());
        denoms.push(Vec::new());
        for _index in 0..loci {
            sums[cluster].push(0.1);
            denoms[cluster].push(0.2); // psuedocounts
        }
    }

    // now lets do some EM
    let log_prior: f32 = (1.0 / (params.ploidy as f32)).ln();
    let mut iterations = 0;

    let mut total_log_loss = f32::NEG_INFINITY;

    let mut final_log_probabilities: Vec<Vec<f32>> = Vec::new();
    for _read in hic_links.iter().chain(hifi.iter()).chain(txg.iter()) {
        final_log_probabilities.push(Vec::new());
    }

    let log_loss_change_limit = 0.000001 * (hic_links.len() as f32); // TODO fiddle with this in testing
    let mut log_loss_change;
    let mut last_log_loss = f32::NEG_INFINITY;
    let mut hic_probabilities: Vec<Vec<f32>> = Vec::new(); // TODO REMOVE DEBUG
    let mut hic_likelihoods: Vec<f32> = Vec::new();

    //while log_loss_change > log_loss_change_limit && iterations < 100 {
    while iterations < 150 {
        // TODO figure out something better here
        hic_probabilities.clear(); // TODO REMOVE DEBUG
        hic_likelihoods.clear();
        let mut log_likelihood = 0.0;
        reset_sums_denoms(loci, &mut sums, &mut denoms, params.ploidy);
        for (readdex, read) in hic_links
            .iter()
            .chain(hifi.iter())
            .chain(txg.iter())
            .enumerate()
        {
            let log_likelihoods = bernoulli_likelihood(read, &cluster_centers, log_prior);
            let read_likelihood = log_sum_exp(&log_likelihoods);
            hic_likelihoods.push(read_likelihood);
            //if iterations > 140 && read_likelihood < -1.0 { continue; } // TODO DEBUG NOT SURE
            log_likelihood += read_likelihood;
            let probabilities = normalize_in_log(&log_likelihoods);

            update_sums_denoms(&mut sums, &mut denoms, read, &probabilities);
            hic_probabilities.push(probabilities);
            final_log_probabilities[readdex] = log_likelihoods;
        }

        /*
        if iterations == 120 {

            for loci in 0..cluster_centers.clusters[0].center.len() {
                let mut total_likelihood = 0.0;
                let reads = variant_hic_reads[loci].len() as f32;
                for read in variant_hic_reads[loci].iter() {
                    let log_likelihoods = bernoulli_likelihood(read, &cluster_centers, log_prior, &locus_filter);
                    let read_likelihood = log_sum_exp(&log_likelihoods);
                    total_likelihood += read_likelihood;
                }
                //println!("{}\t{}", total_likelihood, reads);
                if reads != 0.0 && total_likelihood / reads < -0.75 {
                    locus_filter[loci] = false;
                }
            }
        }
        */
        total_log_loss = log_likelihood;
        log_loss_change = log_likelihood - last_log_loss;
        last_log_loss = log_likelihood;

        update_cluster_centers(loci, &sums, &denoms, &mut cluster_centers);
        iterations += 1;
        eprintln!(
            "bernoulli\t{}\t{}\t{}\t{}\t{}",
            thread_num, epoch, iterations, log_likelihood, log_loss_change
        );
    }
    (
        total_log_loss,
        cluster_centers,
        hic_probabilities,
        hic_likelihoods,
    )
}

fn update_cluster_centers(
    loci: usize,
    sums: &Vec<Vec<f32>>,
    denoms: &Vec<Vec<f32>>,
    cluster_centers: &mut ClusterCenters,
) {
    for locus in 0..loci {
        for cluster in 0..sums.len() {
            if denoms[cluster][locus] > 1.0 {
                let update = sums[cluster][locus] / denoms[cluster][locus];
                cluster_centers.clusters[cluster].center[locus] = update.min(0.9999).max(0.0001);
                //max(0.0001, min(0.9999, update));
            }
        }
    }
}

fn update_sums_denoms(
    sums: &mut Vec<Vec<f32>>,
    denoms: &mut Vec<Vec<f32>>,
    hic_read: &Molecule,
    probabilities: &Vec<f32>,
) {
    for locus in 0..hic_read.loci.len() {
        for (cluster, probability) in probabilities.iter().enumerate() {
            match hic_read.alleles[locus] {
                Allele::Ref => sums[cluster][hic_read.loci[locus]] += 0.0,
                Allele::Alt => sums[cluster][hic_read.loci[locus]] += probability, // * hic_read.long_weighting,
            }
            denoms[cluster][hic_read.loci[locus]] += probabilities[cluster]; // * hic_read.long_weighting;
        }
    }
}

fn normalize_in_log(log_probs: &Vec<f32>) -> Vec<f32> {
    // takes in a log_probability vector and converts it to a normalized probability
    let mut normalized_probabilities: Vec<f32> = Vec::new();
    let sum = log_sum_exp(log_probs);
    for i in 0..log_probs.len() {
        normalized_probabilities.push((log_probs[i] - sum).exp());
    }
    normalized_probabilities
}

fn bernoulli_likelihood(
    read: &Molecule,
    cluster_centers: &ClusterCenters,
    log_prior: f32,
) -> Vec<f32> {
    let mut log_probabilities: Vec<f32> = Vec::new();
    for (cluster, center) in cluster_centers.clusters.iter().enumerate() {
        log_probabilities.push(log_prior);
        for (locus_index, locus) in read.loci.iter().enumerate() {
            log_probabilities[cluster] +=
                ln_bernoulli(read.alleles[locus_index], center.center[*locus]);
        }
    }

    log_probabilities
}

fn ln_bernoulli(allele: Allele, p: f32) -> f32 {
    match allele {
        Allele::Ref => (1.0 - p).ln(),
        Allele::Alt => p.ln(),
    }
}

fn log_sum_exp(p: &Vec<f32>) -> f32 {
    let max_p: f32 = p.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    let sum_rst: f32 = p.iter().map(|x| (x - max_p).exp()).sum();
    max_p + sum_rst.ln()
}

fn reset_sums_denoms(
    loci: usize,
    sums: &mut Vec<Vec<f32>>,
    denoms: &mut Vec<Vec<f32>>,
    num_clusters: usize,
) {
    for cluster in 0..num_clusters {
        for index in 0..loci {
            sums[cluster][index] = 0.10;
            denoms[cluster][index] = 0.20;
        }
    }
}

fn init_cluster_centers(loci: usize, params: &Params, rng: &mut StdRng) -> ClusterCenters {
    let mut centers: ClusterCenters = ClusterCenters {
        clusters: Vec::new(),
    };
    for _cluster in 0..params.ploidy {
        let mut center = ClusterCenter { center: Vec::new() };
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
    txg_mols: Vec<String>,
    hic_mols: Vec<String>,
    ccs_mols: Vec<String>,
    output: String,
    assembly_kmers: String,
    assembly_fasta: String,
    threads: usize,
    seed: u8,
    ploidy: usize,
    restarts: u32,
    min_hic_links: u32,
    break_window: usize,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();

    let het_kmers = params.value_of("het_kmers").unwrap();
    let output = params.value_of("output").unwrap();
    let txg_tmp = match params.values_of("linked_read_mols") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut txg_mols: Vec<String> = Vec::new();
    for x in txg_tmp {
        txg_mols.push(x.to_string());
    }
    let hic_tmp = match params.values_of("hic_mols") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut hic_mols: Vec<String> = Vec::new();
    for x in hic_tmp {
        hic_mols.push(x.to_string());
    }

    let long_tmp = match params.values_of("long_read_mols") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut ccs_mols: Vec<String> = Vec::new();
    for x in long_tmp {
        ccs_mols.push(x.to_string());
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

    let break_window = params.value_of("break_window").unwrap_or("500");
    let break_window = break_window.to_string().parse::<usize>().unwrap();
    eprintln!("break window {}", break_window);

    Params {
        het_kmers: het_kmers.to_string(),
        output: output.to_string(),
        txg_mols: txg_mols,
        hic_mols: hic_mols,
        ccs_mols: ccs_mols,
        assembly_kmers: assembly_kmers.to_string(),
        assembly_fasta: assembly_fasta.to_string(),
        threads: threads,
        seed: seed,
        restarts: restarts,
        ploidy: ploidy,
        min_hic_links: min_hic_links,
        break_window: break_window,
    }
}

//fn unwrap_param(name: &str, default: Option<&str>, )
