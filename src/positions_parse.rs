use crate::vcf_structs::SiteRow;
use crate::vcf_structs::VCFData;
use std::collections::HashSet;

pub fn get_intersection(first: &Vec<SiteRow>, second: &Vec<SiteRow>)-> Vec<SiteRow> {
    let mut main_set = HashSet::new();

    for row in first {
        let chr = row.chromosome.clone();
        let pos = row.position;
        main_set.insert((chr,pos));
    }

    let mut final_vals : Vec<SiteRow> = Vec::new();

    for row in second {
        let chr = row.chromosome.clone();
        let pos = row.position;

        if main_set.contains(&(chr,pos)) {
            final_vals.push(row.clone());
        }
    }

    return final_vals;



}

pub fn keep_sites(sites_to_keep: &Vec<SiteRow>, vcf_data: &VCFData) -> VCFData {
    let mut new_vcf = vcf_data.clone();

    let mut main_set = HashSet::new();

    for row in sites_to_keep {
        let chr = row.chromosome.clone();
        let pos = row.position;
        main_set.insert((chr,pos));
    }

    for j in (0..vcf_data.positions.len()) {
        let position = vcf_data.positions[j];
        let chromosome = vcf_data.chromosomes[j].clone();

        if !main_set.contains(&(chromosome,position)) {
            for i in (0..vcf_data.haplotype_names.len()) {
                new_vcf.vcf_data[i][j] = 255;
            }
        }
    }

    return new_vcf;
}