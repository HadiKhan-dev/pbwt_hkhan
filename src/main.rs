#![allow(warnings, unused)]

use rand::Rng;
use rand::thread_rng;
use rand::seq::SliceRandom;
use vcf_structs::SiteRow;
use vcf_structs::VCFData;

use std::collections::HashMap;
use std::thread::available_parallelism;
use std::sync::{Arc,Mutex,mpsc};
use std;

pub mod fasta_loader;

pub mod pbwt;

pub mod spaced_pbwt;

pub mod pbwt_structs;

pub mod imputer;

pub mod vcf_loader;

pub mod sites_loader;

pub mod vcf_writer;

pub mod vcf_structs;

pub mod positions_parse;

pub mod storer;

pub mod comparison;

use crate::storer::Basic;


fn main() {

    let panel_vcf = vcf_loader::read("./vcf_data/omni4k-10.vcf.gz").unwrap();
    
    let test_vcf = vcf_loader::read("./vcf_data/omni10.vcf.gz").unwrap();
    
    let read_sites = sites_loader::read_sites("./vcf_data/omni10.sites").unwrap();
    
    let illu_sites = sites_loader::read_sites("./vcf_data/illu1M.sites").unwrap();
    

    let kept_sites = positions_parse::get_intersection(&read_sites,&illu_sites);


    let new_test_vcf = positions_parse::keep_sites(&kept_sites,&test_vcf);

    let arc_sites = Arc::new(kept_sites);

    let arc_panel = Arc::new(panel_vcf);

    let now = std::time::Instant::now();

    let dual_panel_pbwt = spaced_pbwt::dual_pbwt(arc_panel,arc_sites,100);

    let elapsed = now.elapsed();

    println!("PBWT Time: {:.4?}", elapsed);

    let mut imputed_col_flags = dual_panel_pbwt.forward_pbwt.pbwt_col_flags.clone();

    for i in 0..imputed_col_flags.len() {
        if imputed_col_flags[i] == 0 {
            imputed_col_flags[i] = 1;
        } else {
            imputed_col_flags[i] = 0;
        }
    }

    //storer::write_pbwt(&panel_pbwt,"outputs/test_pbwt.pbwt");

    // //println!("{:?}",new_test_vcf.vcf_data);

    // let panel_pbwt = storer::read_pbwt("outputs/test_pbwt.pbwt");

    // let now = std::time::Instant::now();

    let mut freqs : Vec<f64> = Vec::with_capacity(dual_panel_pbwt.forward_pbwt.num_total_sites as usize);

    for i in 0..dual_panel_pbwt.forward_pbwt.num_total_sites {
        let zero_freq = (dual_panel_pbwt.forward_pbwt.count[i as usize] as f64)/(dual_panel_pbwt.forward_pbwt.num_samples as f64);
        freqs.push(100.0*(1.0-zero_freq));
    }

    let mut a = new_test_vcf.vcf_data.clone();
    let mut b = test_vcf.vcf_data.clone();

    // for i in 0..99 {
    //     let mut m = new_test_vcf.vcf_data.clone();
    //     let mut n = test_vcf.vcf_data.clone();
    //     a.append(&mut m);
    //     b.append(&mut n);
    // }

    let am = Arc::new(dual_panel_pbwt);
    let bm = Arc::clone(&am);

    let now = std::time::Instant::now();


    let imputed = imputer::new_impute(am,a,8);

    let elapsed = now.elapsed();
    
    let dual_panel_pbwt = Arc::<pbwt_structs::DualPbwt>::try_unwrap(bm).unwrap();

    let bucket_bounds = vec![0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0, 50.0, 70.0, 90.0, 100.0];

    comparison::compare_results(&b,&imputed,&imputed_col_flags, &freqs,&bucket_bounds);


    println!("Time total: {:.4?}", elapsed);

    // let  X =
    // vec![
    // vec![0, 1, 0, 1, 0, 1],
    // vec![1, 1, 0, 0, 0, 1],
    // vec![1, 1, 1, 1, 1, 1],
    // vec![0, 1, 1, 1, 1, 0],
    // vec![0, 0, 0, 0, 1, 1],
    // vec![1, 0, 0, 0, 1, 0],
    // vec![1, 1, 0, 0, 0, 1],
    // vec![0, 1, 0, 1, 1, 0]];

    // let Y = X.clone();

    // let l = X.len();
    // let m = X[0].len();

    // let vcf = VCFData {
    //     vcf_data: X,
    //     chromosomes : vec![String::from("1"); m],
    //     positions: (1..((m+1) as u32)).collect(),

    //     sample_names: vec![String::from("A"),String::from("B"),
    //     String::from("C"),String::from("D")],

    //     haplotype_names: vec![String::from("A[0]"),String::from("A[1]"),String::from("B[0]"),
    //     String::from("B[1]"),String::from("C[0]"),String::from("C[1]"),String::from("D[0]"),
    //     String::from("D[1]")],

    //     alt_allele_freq: vec![0.0;m],
    // };

    // let test = vec![0,0,0,255,255,255];

    // let mut keep_rows = Vec::new();

    // for i in vec![0,1,2] {
    //     let new_sr = SiteRow {
    //         chromosome: String::from("1"),
    //         position: (i+1) as u32,
    //         reference: String::from("A"),
    //         alternate: String::from("C"),
    //     };
    //     keep_rows.push(new_sr);
    // }

    // let p = spaced_pbwt::spaced_pbwt(&vcf,&keep_rows,10);
    // let l = spaced_pbwt::dual_pbwt(&vcf,&keep_rows,1);

    // let ins = spaced_pbwt::spaced_insert_place(&p, &test);

    // let impd = imputer::new_impute_single(&l,&test);
    // //let recd = spaced_pbwt::spaced_recover_sequence(&p,7);
    // println!("{:?}",impd);


}



