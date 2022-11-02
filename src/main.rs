#![allow(warnings, unused)]

use rand::Rng;
use rand::thread_rng;
use rand::seq::SliceRandom;

use std::collections::HashMap;

pub mod fasta_loader;

pub mod pbwt;

pub mod imputer;

pub mod vcf_loader;

pub mod sites_loader;

pub mod vcf_writer;

pub mod vcf_structs;
fn main() {

/*     let mut reference_data = fasta_loader::remove_all_zeros(&fasta_loader::bool_leveling(&fasta_loader::get_variations(&vec!["testing.fasta"],true)));
    let mut test_sequence = reference_data.remove(199);
    //let mut test_sequence = reference_data[665].clone();
    let test_base = test_sequence.clone();

    println!("Loaded Data");

    let pbwt_dat = pbwt::pbwt(&reference_data,100);

    println!("Computed PBWT");

    let mut tot_changed = 0;
    let mut correct = 0;
    let mut zero_correct = 0;
    let mut one_correct = 0;
    let mut zero_imp = 0;
    let mut one_imp = 0;
    let mut zero_hit = 0;
    let mut one_hit = 0;

    let mut rand_muts: Vec<usize> = Vec::new();
    let mut old_vals: HashMap<usize,u8> = HashMap::new();
    
    let ratio = 0.7;
    let limit = ((ratio*(test_sequence.len() as f64)).ceil() as usize);

    println!("{}",limit);
    let mut shuf: Vec<usize> = (0..test_sequence.len()).collect();
    shuf.shuffle(&mut thread_rng());


    for i in 0..limit {     
        //let r_val: usize = rand::thread_rng().gen_range(0..test_sequence.len());
        let r_val: usize = shuf[i];


        if !old_vals.contains_key(&r_val) {
            rand_muts.push(r_val);
            old_vals.insert(r_val,test_sequence[r_val]);
            test_sequence[r_val] = 255;
            tot_changed += 1;
        }
    }

    let imp = imputer::impute(&reference_data, &pbwt_dat, &test_sequence);


    for i in 0..rand_muts.len() {
        let cur_point = rand_muts[i];
        let imputed_value = imp[cur_point];
        let old_value = *old_vals.get(&cur_point).unwrap();

        if imputed_value == old_value {
            correct += 1;
        } 
        if imputed_value == 0 && old_value == 0 {
            zero_hit += 1;
        }
        if imputed_value == 1 && old_value == 1 {
            one_hit += 1;
        }
        if old_value == 0 {
            zero_correct += 1;
        }
        if old_value == 1 {
            one_correct += 1;
        }
        if imputed_value == 0 {
            zero_imp += 1;
        }
        if imputed_value == 1 {
            one_imp += 1;
        }
    }

    //for i in 0..imp.len() {
    //    println!("{},{}",test_base[i],imp[i]);
    //}

    println!("");
    println!("Length: {}",test_sequence.len());
    println!("Percentage Imputed: {}",(tot_changed as f64)/(test_sequence.len() as f64));
    println!("");
    println!("Correct: {} out of total: {}",correct,tot_changed);
    println!("Ratio: {:.2}%",100.0*(correct as f64)/(tot_changed as f64));

    println!("");
    println!("Percent Zero: {}",100.0*(zero_correct as f64)/(tot_changed as f64));
    println!("Percent One: {}",100.0*(one_correct as f64)/(tot_changed as f64));
    //println!("{:?}",imp);

    println!("");
    println!("Imp Percent Zero: {}",100.0*(zero_imp as f64)/(tot_changed as f64));
    println!("Imp Percent One: {}",100.0*(one_imp as f64)/(tot_changed as f64));

    println!("");
    println!("Correct Zero Percentage: {}",100.0*(zero_hit as f64)/(zero_correct as f64));
    println!("Correct One Percentage: {}",100.0*(one_hit as f64)/(one_correct as f64));

    println!("");
    println!("Precision: {}",(one_hit as f64)/(one_imp as f64));
    println!("Recall: {}",(one_hit as f64)/(one_correct as f64)); */

    let read_vcf = vcf_loader::read("./vcf_data/omni10.vcf.gz").unwrap();

    let read_sites = sites_loader::read_sites("./vcf_data/omni10.sites").unwrap();
    
    vcf_writer::write_vcf(read_vcf,read_sites,"hkhan.vcf");

}




