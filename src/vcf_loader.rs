use vcf::{VCFReader, U8Vec, VCFHeaderFilterAlt, VCFError, VCFRecord};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::BufReader;
use std;

use crate::vcf_structs::VCFData;

pub fn read(filename: &str) -> Result<VCFData, VCFError> {


    let mut reader = VCFReader::new(BufReader::new(MultiGzDecoder::new(File::open(
        filename,
    )?)))?;

    
    let mut sample_names = Vec::new();
    let mut haplotype_names = Vec::new();


    for i in reader.header().samples() {
        let cur_sample = String::from_utf8(i.to_vec()).unwrap();

        let mut first = cur_sample.clone();
        let mut second = cur_sample.clone();

        first.push_str("[0]");
        second.push_str("[1]");

        sample_names.push(cur_sample);
        haplotype_names.push(first);
        haplotype_names.push(second);
    }


    let mut panel_matrix: Vec<Vec<u8>> = Vec::with_capacity(haplotype_names.len());

    for i in 0..haplotype_names.len() {
        panel_matrix.push(Vec::new());
    }

    let mut vcf_record = reader.empty_record();

    let mut last_position : u64 = 0;
    let mut cur_position : u64 = 0;

    let mut started = false;
    let mut finished = false;

    let mut position_list: Vec<u32> = Vec::new();

    let mut chromosome_list = Vec::new();


    while !started |!finished {

        last_position = cur_position;


        reader.next_record(&mut vcf_record);

        cur_position = vcf_record.position.clone();
        


        if started & (cur_position == last_position) {
            finished = true;
            break;
        }
        started = true;

        //let chr_u64: u64 = String::from_utf8(vcf_record.chromosome.clone()).unwrap().parse().unwrap();
        


        let chr_str: String = String::from_utf8(vcf_record.chromosome.clone()).unwrap();
        

        chromosome_list.push(chr_str);
        position_list.push(cur_position as u32);

        let mut current = 0;



        
        for sam in &sample_names {

            let byteval = sam.as_bytes();
            let vals = vcf_record.genotype(byteval,b"GT").unwrap();
            let value = String::from_utf8(vals[0].clone()).unwrap();
            let first_value = char::to_digit(value.as_bytes()[0] as char,10).unwrap() as u8;
            let second_value = char::to_digit(value.as_bytes()[2] as char,10).unwrap() as u8;
            let first_push = 2*current;
            let second_push = 2*current+1;

            panel_matrix[first_push].push(first_value);
            panel_matrix[second_push].push(second_value);

            current += 1;

        }

    }

    let mut mafs: Vec<f64> = Vec::with_capacity(position_list.len());

    for j in (0..mafs.len()) {
        let mut ones : f64 = 0.0;
        let mut tot : f64 = 0.0;

        for i in 0..haplotype_names.len() {
            if panel_matrix[i][j] == 1 {
                ones += 1.0;
            }
            tot += 1.0;
        }

        let mut rat = ones/tot;

        mafs.push(rat);
    }

    let loaded = VCFData { vcf_data: panel_matrix, chromosomes: chromosome_list,
        positions: position_list, sample_names: sample_names,
         haplotype_names: haplotype_names, alt_allele_freq : mafs};

    return Ok(loaded);
}

pub fn reverse_vcf(vcf: &VCFData) -> VCFData {
    let mut rev_vcf_data = vcf.vcf_data.clone();
    for i in 0..rev_vcf_data.len() {
        rev_vcf_data[i].reverse();
    }

    let mut rev_chromosomes = vcf.chromosomes.clone();
    rev_chromosomes.reverse();

    let mut rev_positions = vcf.positions.clone();
    rev_positions.reverse();

    let rev_sample_names = vcf.sample_names.clone();

    let rev_haplotype_names = vcf.haplotype_names.clone();

    let mut rev_alt_allele_freq = vcf.alt_allele_freq.clone();
    rev_alt_allele_freq.reverse();

    return VCFData {
        vcf_data: rev_vcf_data,
        chromosomes: rev_chromosomes,
        positions: rev_positions,
        sample_names: rev_sample_names,
        haplotype_names: rev_haplotype_names,
        alt_allele_freq: rev_alt_allele_freq,
    }

}