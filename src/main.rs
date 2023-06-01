#![allow(warnings, unused)]


use xgboost_rs;

use rand::Rng;
use rand::thread_rng;
use rand::seq::SliceRandom;
use vcf_structs::SiteRow;
use vcf_structs::VCFData;


use std::collections::HashMap;
use std::thread::available_parallelism;
use std::sync::{Arc,Mutex,mpsc};
use std;

use crate::storer::Basic;

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

pub mod pca_weights;

pub mod helper_structs;

fn main() {

    // let mut data: Vec<Vec<f32>> = 
    // vec![vec![0.0,0.1,0.12,0.2,5.0,8.0,10.0,3.0],
    // vec![0.0,0.1,0.12,0.2,5.0,8.0,10.0,3.0]];

    // let num_rows = data.len();

    // let mut flat_data = flatten_matrix(data);

    // fix_zeros(&mut flat_data);


    // let mut d_matrix = xgboost_rs::DMatrix::from_dense(&mut flat_data,num_rows).unwrap();

    let model_0dot1 = xgboost_rs::Booster::load(
        "../pbwt_python/xgboost_models/model_2.json").unwrap();
    let model_0dot2 = xgboost_rs::Booster::load(
        "../pbwt_python/xgboost_models/model_2.json").unwrap();
    let model_0dot3 = xgboost_rs::Booster::load(
        "../pbwt_python/xgboost_models/model_2.json").unwrap();
    let model_0dot5 = xgboost_rs::Booster::load(
        "../pbwt_python/xgboost_models/model_2.json").unwrap();
    let model_0dot7 = xgboost_rs::Booster::load(
        "../pbwt_python/xgboost_models/model_2.json").unwrap();
    let model_1 = xgboost_rs::Booster::load(
        "../pbwt_python/xgboost_models/model_2.json").unwrap();
    let model_2 = xgboost_rs::Booster::load(
        "../pbwt_python/xgboost_models/model_2.json").unwrap();
    let model_3 = xgboost_rs::Booster::load(
        "../pbwt_python/xgboost_models/model_3.json").unwrap();
    let model_5 = xgboost_rs::Booster::load(
        "../pbwt_python/xgboost_models/model_5.json").unwrap();
    let model_7 = xgboost_rs::Booster::load(
        "../pbwt_python/xgboost_models/model_7.json").unwrap();
    let model_10 = xgboost_rs::Booster::load(
        "../pbwt_python/xgboost_models/model_10.json").unwrap();
    let model_20 = xgboost_rs::Booster::load(
        "../pbwt_python/xgboost_models/model_20.json").unwrap();
    let model_30 = xgboost_rs::Booster::load(
        "../pbwt_python/xgboost_models/model_30.json").unwrap();
    let model_50 = xgboost_rs::Booster::load(
        "../pbwt_python/xgboost_models/model_50.json").unwrap();
    let model_70 = xgboost_rs::Booster::load(
        "../pbwt_python/xgboost_models/model_70.json").unwrap();
    let model_90 = xgboost_rs::Booster::load(
        "../pbwt_python/xgboost_models/model_90.json").unwrap();
    let model_100 = xgboost_rs::Booster::load(
        "../pbwt_python/xgboost_models/model_100.json").unwrap();
    let basic_model = xgboost_rs::Booster::load(
        "../pbwt_python/xgboost_models/basic_model.json").unwrap();

    let xgboost_models = vec![model_0dot1,model_0dot2,model_0dot3,model_0dot5,
                                        model_0dot7,model_1,model_2,model_3,model_5,model_7,
                                        model_10,model_20,model_30,model_50,model_70,
                                        model_90,model_100];

    // let xgboost_models = vec![basic_model];

    let arc_models = Arc::new(xgboost_models);
    
    // println!("Predictions: {:?}",model.predict(&d_matrix));


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

    let mut freqs : Vec<f64> = Vec::with_capacity(dual_panel_pbwt.forward_pbwt.num_total_sites as usize);

    for i in 0..dual_panel_pbwt.forward_pbwt.num_total_sites {
        let zero_freq = (dual_panel_pbwt.forward_pbwt.count[i as usize] as f64)/(dual_panel_pbwt.forward_pbwt.num_samples as f64);
        freqs.push(100.0*(1.0-zero_freq));
    }

    let mut a = new_test_vcf.vcf_data.clone();
    let mut b = test_vcf.vcf_data.clone();

    for i in 0..0 {
        let mut m = new_test_vcf.vcf_data.clone();
        let mut n = test_vcf.vcf_data.clone();
        a.append(&mut m);
        b.append(&mut n);
    }

    let am = Arc::new(dual_panel_pbwt);
    let bm = Arc::clone(&am);

    let now = std::time::Instant::now();

    // let imp = imputer::impute_single_xgboost(&am,&arc_models,&a[0]);
    let imputed = imputer::impute_xgboost(am,arc_models,a,8);

    let elapsed = now.elapsed();
    
    // let dual_panel_pbwt = Arc::<pbwt_structs::DualPbwt>::try_unwrap(bm).unwrap();

    let bucket_bounds = vec![0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0, 50.0, 70.0, 90.0, 100.0];

    comparison::compare_results(&b,&imputed,&imputed_col_flags, &freqs,&bucket_bounds);


    println!("Time total: {:.4?}", elapsed);

}



