use std::io;
use vcf::*;
use indexmap::indexmap;
use std::str;
use std::fs::File;

use crate::vcf_structs;


pub fn write_vcf(vcf_data: vcf_structs::VCFData, sites_data: Vec<vcf_structs::SiteRow>, filename: &str) -> () {

    let mut write_buf: Vec<u8> = Vec::new();
    let sample_names = vcf_data.sample_names;
    let mut u8samples: Vec<Vec<u8>> = Vec::with_capacity(sample_names.len());

    for name in sample_names {
        let bytename = name.into_bytes();
        u8samples.push(bytename);
    }

    let u8copy = u8samples.clone();

    let mut header_line = VCFHeaderLine::from_bytes(b"##Imputed VCF File \n",0).unwrap();
    
    let header = VCFHeader::new(vec![header_line],u8samples);

    let mut writer = VCFWriter::new(File::create(filename).unwrap(),&header).unwrap();

    let panel = &vcf_data.vcf_data;
    let n = panel.len();
    let m = panel[0].len();

    let mut j = 0;
    
    for row in sites_data {

        let genotypes = n/2;


        let chr = row.chromosome;
        let position = row.position;
        let reference = row.reference;
        let alternative = row.alternate;
        let format = "GT".to_string();
        let an = n;
        let mut ac = 0;

        let mut new_rec = VCFRecord::new(header.clone());
        new_rec.chromosome = chr.to_string().into_bytes();
        new_rec.position = position as u64;
        new_rec.reference = reference.into_bytes();
        new_rec.alternative = vec![alternative.into_bytes()];
        new_rec.filter = vec!["PASS".to_string().into_bytes()];
        new_rec.format = vec![format.into_bytes()];

        let mut full_geno = vec![];

        for i in (0..genotypes) {
            let first_i = i*2;
            let second_i = i*2+1;

            let first_val = panel[first_i][j];
            let second_val = panel[second_i][j];

            if first_val == 1 {
                ac += 1;
            }
            if second_val == 1 {
                ac += 1;
            }

            let sam_name = &u8copy[i]; 

            let full_vals = vec![vec![first_val.to_string().into_bytes()[0],124,second_val.to_string().into_bytes()[0]]];
            full_geno.push(vec![full_vals]);
            
        }

        let info = vec![(vec![65, 67], vec![ac.to_string().into_bytes()]), (vec![65, 78], vec![an.to_string().into_bytes()])];
        
        new_rec.genotype = full_geno;
        new_rec.info = info;

        writer.write_record(&new_rec);

        j += 1

    }

    return ()
}