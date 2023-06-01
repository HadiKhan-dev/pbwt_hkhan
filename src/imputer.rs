use crate::pbwt;
use crate::pbwt_structs;
use crate::pbwt_structs::{SpacedPbwt, DualPbwt};
use crate::spaced_pbwt::spaced_insert_place;
use crate::pca_weights;
use core::slice::Iter;
use std::cmp;
use std::cmp::max;
use std::cmp::min;
use std::thread;
use std::collections::HashMap;
use std::sync::{Arc,Mutex,mpsc};
use closure::closure;
use threadpool::ThreadPool;
use xgboost_rs;


struct RowWPos{
    pub position: u32,
    pub row: Vec<u8>,
}

fn fix_zeros(float_data: &mut [f32]) {
    /*
    Replace zeros with epsilon to stop XGBoost from considering them as NaN values
     */
    let ZERO: f32 = 1e-15;
    for x in float_data {
        if *x == 0.0 {
            *x = ZERO;
        }
    }
}

fn flatten_matrix(data: Vec<Vec<f32>>) -> Vec<f32> {
    /*
    Flatten a 2D matrix to 1D
     */
    return data.into_iter().flatten().collect();
}


pub fn impute_single_pca(pbwt_data: &DualPbwt, test_sequence: &Vec<u8>) -> Vec<u8> {

    let buckets : Vec<f64> = vec![0.1,0.2,0.3,0.5,0.7,1.0,2.0,3.0,
    5.0,7.0,10.0,20.0,30.0,50.0,70.0,90.0,100.0];

    let prob_cutoff : Vec<f64> = vec![0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
    0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5];

    let side_distance = 10;
    let mut final_sequence: Vec<u8> = Vec::new();
    let m = test_sequence.len();
    let panel_size = pbwt_data.forward_pbwt.num_samples as usize;

    let mut test_rev = test_sequence.clone();
    test_rev.reverse();

    let mut forward_places = spaced_insert_place(&pbwt_data.forward_pbwt, test_sequence,10,0);
    let mut forward_positions = forward_places.insert_positions;
    let mut forward_side_data = forward_places.side_values;
    let mut forward_divergence = forward_places.divergence_values;
    
    let mut reverse_places = spaced_insert_place(&pbwt_data.reverse_pbwt, &test_rev,10,0);
    let mut reverse_positions = reverse_places.insert_positions;
    let mut reverse_side_data = reverse_places.side_values;
    let mut reverse_divergence = reverse_places.divergence_values;
    
    reverse_positions.reverse();
    reverse_side_data.reverse();
    reverse_divergence.reverse();


    for i in 0..m {
        if test_sequence[i] != 255 {
            final_sequence.push(test_sequence[i]);
        } else {

            let mut loc;

            let one_freq = 100.0*(1.0-(pbwt_data.forward_pbwt.count[i] as f64)/(pbwt_data.forward_pbwt.num_samples as f64));
            
            let bucket_position = buckets.binary_search_by(|v| {
                v.partial_cmp(&one_freq).expect("Couldn't compare values")
            });

            match bucket_position {
                Ok(i) => {
                    loc = i;
                }
                Err(i) => {
                    loc = i;
                }
            }

            let bucket_midpoint: f64;

            if loc == 0 {
                bucket_midpoint = buckets[0]/2.0;
            } else {
                bucket_midpoint = (buckets[loc]+buckets[loc-1])/2.0;
            }

            let mut made_coeffs: Vec<f64> = vec![0.0; side_distance as usize];

            let place_forward = forward_positions[i];
            let place_reverse = reverse_positions[i];

            let forward_vals_lower = &forward_side_data[i][0];
            let forward_vals_upper = &forward_side_data[i][1];

            let reverse_vals_lower = &reverse_side_data[i][0];
            let reverse_vals_upper = &reverse_side_data[i][1];



            for i in 0..(side_distance as usize) {
                let forward_lower_val: f64;
                let forward_upper_val: f64;
                let reverse_lower_val: f64;
                let reverse_upper_val: f64;

                if i < forward_vals_lower.len() {
                    forward_lower_val = forward_vals_lower[i] as f64; 
                } else {
                    forward_lower_val = bucket_midpoint;
                }

                if i < forward_vals_upper.len() {
                    forward_upper_val = forward_vals_upper[i] as f64; 
                } else {
                    forward_upper_val = bucket_midpoint;
                }

                if i < reverse_vals_lower.len() {
                    reverse_lower_val = reverse_vals_lower[i] as f64; 
                } else {
                    reverse_lower_val = bucket_midpoint;
                }

                if i < reverse_vals_upper.len() {
                    reverse_upper_val = reverse_vals_upper[i] as f64; 
                } else {
                    reverse_upper_val = bucket_midpoint;
                }

                if i < made_coeffs.len() {
                    made_coeffs[i] += 0.25*forward_lower_val;
                    made_coeffs[i] += 0.25*forward_upper_val;
                    made_coeffs[i] += 0.25*reverse_lower_val;
                    made_coeffs[i] += 0.25*reverse_upper_val;
                    
                }

            }

            let prob: f64 = pca_weights::get_probability(&made_coeffs,loc);

            if prob > prob_cutoff[loc] {
                final_sequence.push(1);
            } else {
                final_sequence.push(0);
            }

        }
    }

    return final_sequence;
}

pub fn impute_single_xgboost(pbwt_data: &DualPbwt, xgboost_models: &Vec<xgboost_rs::Booster>,
    test_sequence: &Vec<u8>) -> Vec<u8> {


    fn small_test(x: f32,y: f32) -> u8 {
        if x <= y {
            return 0;
        }
        return 1;
    }

    let side_length: usize = 10;
    let divergence_length: usize = 10;

    let buckets : Vec<f32> = vec![0.1,0.2,0.3,0.5,0.7,1.0,2.0,3.0,
    5.0,7.0,10.0,20.0,30.0,50.0,70.0,90.0,100.0];

    let prob_cutoff : Vec<f32> = vec![0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
    0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5];



    // let xgb_models = vec![model.clone(); buckets.len()];

    let mut impute_inputs: Vec<Vec<f32>> = vec![Vec::new(); buckets.len()];

    let mut impute_counts: Vec<usize> = vec![0; buckets.len()];

    let mut middle_sequence: Vec<u8> = Vec::new();

    let mut final_sequence: Vec<u8> = Vec::new();

    let m = test_sequence.len();
    let panel_size = pbwt_data.forward_pbwt.num_samples as usize;

    let mut test_rev = test_sequence.clone();
    test_rev.reverse();

    let mut forward_places = spaced_insert_place(&pbwt_data.forward_pbwt, test_sequence,10,10);
    let mut reverse_places = spaced_insert_place(&pbwt_data.reverse_pbwt, &test_rev,10,10);
    
    let mut forward_positions = forward_places.insert_positions;
    let mut forward_side_data = forward_places.side_values;
    let mut forward_divergence = forward_places.divergence_values;
    
    let mut reverse_positions = reverse_places.insert_positions;
    let mut reverse_side_data = reverse_places.side_values;
    let mut reverse_divergence = reverse_places.divergence_values;
    
    reverse_positions.reverse();
    reverse_side_data.reverse();
    reverse_divergence.reverse();


    for i in 0..m {
        if test_sequence[i] != 255 {
            middle_sequence.push(test_sequence[i]);
        } else {
            let mut loc;

            let one_freq = 100.0*(1.0-(pbwt_data.forward_pbwt.count[i] as f32)/(pbwt_data.forward_pbwt.num_samples as f32));
            
            let bucket_position = buckets.binary_search_by(|v| {
                v.partial_cmp(&one_freq).expect("Couldn't compare values")
            });

            match bucket_position {
                Ok(i) => {
                    loc = i;
                }
                Err(i) => {
                    loc = i;
                }
            }

            let bucket_midpoint: f32;

            if loc == 0 {
                bucket_midpoint = buckets[0]/2.0;
            } else {
                bucket_midpoint = (buckets[loc]+buckets[loc-1])/2.0;
            }

            let place_forward = forward_positions[i];
            let place_reverse = reverse_positions[i];

            let mut forward_side_lower: Vec<f32> = forward_side_data[i][0].clone().iter().map(|x| *x as f32).collect();
            let mut forward_side_upper: Vec<f32> = forward_side_data[i][1].clone().iter().map(|x| *x as f32).collect();

            let mut reverse_side_lower: Vec<f32> = reverse_side_data[i][0].clone().iter().map(|x| *x as f32).collect();
            let mut reverse_side_upper: Vec<f32> = reverse_side_data[i][1].clone().iter().map(|x| *x as f32).collect();


            let forward_divergence_lower = forward_divergence[i][0].clone();
            let forward_divergence_upper = forward_divergence[i][1].clone();

            let reverse_divergence_lower = reverse_divergence[i][0].clone();
            let reverse_divergence_upper = reverse_divergence[i][1].clone();

            forward_side_lower.append(&mut vec![bucket_midpoint; side_length-forward_side_lower.len()]);
            forward_side_upper.append(&mut vec![bucket_midpoint; side_length-forward_side_upper.len()]);

            reverse_side_lower.append(&mut vec![bucket_midpoint; side_length-reverse_side_lower.len()]);
            reverse_side_upper.append(&mut vec![bucket_midpoint; side_length-reverse_side_upper.len()]);


            let mut full_side_data = [forward_side_lower,forward_side_upper,reverse_side_lower,reverse_side_upper].concat();
            let mut full_divergence_data = [forward_divergence_lower,forward_divergence_upper,reverse_divergence_lower,reverse_divergence_upper].concat().iter().map(|x| *x as f32).collect();
            
            let mut full_data = [full_side_data,full_divergence_data].concat();
            
            fix_zeros(&mut full_data);

            impute_inputs[loc].append(&mut full_data);
            impute_counts[loc] += 1;

            middle_sequence.push(100+loc as u8);
        }
    }

    let mut imputed_values: Vec<Vec<u8>> = Vec::new();

    for j in 0..buckets.len() {
        let d_matrix = xgboost_rs::DMatrix::from_dense(&mut impute_inputs[j],impute_counts[j]).unwrap();
        let impute_numbers = xgboost_models[j].predict(&d_matrix).unwrap();

        imputed_values.push(impute_numbers.iter().map(
            |x| small_test(*x,prob_cutoff[j])).collect());
    }

    let mut imputed_iter: Vec<Iter<u8>> = imputed_values.iter().map(|x| x.iter()).collect();


    for i in 0..m {
        let value = middle_sequence[i];
        if value < 100 {
            final_sequence.push(value)
        } else {
            let location = (value-100) as usize;

            final_sequence.push(*imputed_iter[location].next().unwrap());
        }
    }

    return final_sequence;
}

pub fn impute_pca(pbwt_data : Arc<pbwt_structs::DualPbwt>, test_sequences : Vec<Vec<u8>>,num_threads: usize) -> Vec<Vec<u8>> {

    let N = test_sequences.len();

    let mut fin_imputed: Vec<Vec<u8>> = vec![Vec::new();N];
    let safe_pbwt = Arc::new(pbwt_data);

    //  let xgboost_model: xgboost_rs::Booster = xgboost_rs::Booster::load(
    //      "../pbwt_python/model_dump.json").unwrap();

    //  let xgb_arc = Arc::new(xgboost_model);

    let (tx,rx) = mpsc::channel();

    let pool = ThreadPool::new(num_threads);

    let mut i = 0;
    for tst_seq in test_sequences.into_iter(){

        let txr = tx.clone();
        let safe_test_seq = Arc::new(tst_seq);
        let refer = Arc::clone(&safe_pbwt);
        //let refer_model = Arc::clone(&xgb_arc);

        let clo = closure!(move txr, move refer, move safe_test_seq, || {

            let val = impute_single_pca(&refer,&safe_test_seq);

            let ret = RowWPos { position: i as u32, row: val};

            txr.send(ret).unwrap();
        });

        pool.execute(clo);
        i += 1;

    }

    drop(tx);
    for recd in rx.iter() {

        fin_imputed[recd.position as usize] = recd.row;
    }

    return fin_imputed;
}

pub fn impute_xgboost(pbwt_data : Arc<pbwt_structs::DualPbwt>,
    xgboost_models : Arc<Vec<xgboost_rs::Booster>>, test_sequences : Vec<Vec<u8>>,
    num_threads: usize) -> Vec<Vec<u8>> {

    let N = test_sequences.len();

    let mut fin_imputed: Vec<Vec<u8>> = vec![Vec::new();N];
    let safe_pbwt = Arc::new(pbwt_data);

    let xgb_arc = Arc::new(xgboost_models);

    let (tx,rx) = mpsc::channel();

    let pool = ThreadPool::new(num_threads);

    let mut i = 0;
    for tst_seq in test_sequences.into_iter(){

        let txr = tx.clone();
        let safe_test_seq = Arc::new(tst_seq);
        let refer = Arc::clone(&safe_pbwt);
        let refer_models = Arc::clone(&xgb_arc);

        let clo = closure!(move txr, move refer, move safe_test_seq, move refer_models, || {

            let val = impute_single_xgboost(&refer,&refer_models,&safe_test_seq);

            let ret = RowWPos { position: i as u32, row: val};

            txr.send(ret).unwrap();
        });

        pool.execute(clo);
        i += 1;

    }

    drop(tx);
    for recd in rx.iter() {

        fin_imputed[recd.position as usize] = recd.row;
    }

    return fin_imputed;
}
