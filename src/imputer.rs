use crate::pbwt;
use crate::pbwt_structs;
use crate::pbwt_structs::{SpacedPbwt, DualPbwt};
use crate::spaced_pbwt::spaced_insert_place;
use crate::pca_weights;
use std::cmp;
use std::cmp::max;
use std::cmp::min;
use std::thread;
use std::collections::HashMap;
use std::sync::{Arc,Mutex,mpsc};
use closure::closure;
use threadpool::ThreadPool;


struct RowWPos{
    pub position: u32,
    pub row: Vec<u8>,
}

pub fn impute_single(pbwt_data: &DualPbwt, test_sequence: &Vec<u8>) -> Vec<u8> {

    let side_distance = 10;
    let mut final_sequence: Vec<u8> = Vec::new();
    let m = test_sequence.len();
    let panel_size = pbwt_data.forward_pbwt.num_samples as usize;

    let mut test_rev = test_sequence.clone();
    test_rev.reverse();

    let forward_places = spaced_insert_place(&pbwt_data.forward_pbwt, test_sequence);
    let mut reverse_places = spaced_insert_place(&pbwt_data.reverse_pbwt, &test_rev);
    reverse_places.reverse();

    for i in 0..m {
        if test_sequence[i] != 255 {
            final_sequence.push(test_sequence[i]);
        } else {
            let place_forward = forward_places[i];
            let place_reverse = reverse_places[i];

            let forward_min = max(place_forward-side_distance+1,0) as usize;
            let forward_max = min(place_forward+side_distance+1,panel_size as i32) as usize;
            let reverse_min = max(place_reverse-side_distance+1,0) as usize;
            let reverse_max = min(place_reverse+side_distance+1,panel_size as i32) as usize;


            let forward_vals_lower_ref = &pbwt_data.forward_pbwt.bin_pbwt[i][forward_min..(place_forward+1) as usize];
            let forward_vals_upper_ref = &pbwt_data.forward_pbwt.bin_pbwt[i][(place_forward+1) as usize..forward_max];

            let reverse_vals_lower_ref = &pbwt_data.reverse_pbwt.bin_pbwt[m-i-1][reverse_min..(place_reverse+1) as usize];
            let reverse_vals_upper_ref = &pbwt_data.reverse_pbwt.bin_pbwt[m-i-1][(place_reverse+1) as usize..reverse_max];

            let mut forward_vals_lower = Vec::new();
            let mut forward_vals_upper = Vec::new();

            let mut reverse_vals_lower = Vec::new();
            let mut reverse_vals_upper = Vec::new();

            for i in 0..forward_vals_lower_ref.len() {
                forward_vals_lower.push(forward_vals_lower_ref[forward_vals_lower_ref.len()-1-i]);
            }

            for i in 0..forward_vals_upper_ref.len() {
                forward_vals_upper.push(forward_vals_upper_ref[i]);
            }

            for i in 0..reverse_vals_lower_ref.len() {
                reverse_vals_lower.push(reverse_vals_lower_ref[reverse_vals_lower_ref.len()-1-i]);
            }

            for i in 0..reverse_vals_upper_ref.len() {
                reverse_vals_upper.push(reverse_vals_upper_ref[i]);
            }




            let mut zero_tot = 0;
            let mut one_tot = 0;

            let mut expo: f64 = 0.0;

            let intercept = -7.986;

            let coeffs = vec![16.175,0.957,0.116];

            expo += intercept;


            for i in 0..(side_distance as usize) {
                let forward_lower_val: f64;
                let forward_upper_val: f64;
                let reverse_lower_val: f64;
                let reverse_upper_val: f64;

                if i < forward_vals_lower.len() {
                    forward_lower_val = forward_vals_lower[i] as f64; 
                } else {
                    forward_lower_val = 0.5;
                }

                if i < forward_vals_upper.len() {
                    forward_upper_val = forward_vals_upper[i] as f64; 
                } else {
                    forward_upper_val = 0.5;
                }

                if i < reverse_vals_lower.len() {
                    reverse_lower_val = reverse_vals_lower[i] as f64; 
                } else {
                    reverse_lower_val = 0.5;
                }

                if i < reverse_vals_upper.len() {
                    reverse_upper_val = reverse_vals_upper[i] as f64; 
                } else {
                    reverse_upper_val = 0.5;
                }

                if i < coeffs.len() {
                    expo += 0.25*forward_lower_val*coeffs[i];
                    expo += 0.25*forward_upper_val*coeffs[i];
                    expo += 0.25*reverse_lower_val*coeffs[i];
                    expo += 0.25*reverse_upper_val*coeffs[i];
                }


                if forward_lower_val == 0.0 {
                    zero_tot += 1;
                } else if forward_lower_val == 1.0 {
                    one_tot += 1;
                }

                if forward_upper_val == 0.0 {
                    zero_tot += 1;
                } else if forward_upper_val == 1.0 {
                    one_tot += 1;
                }

                if reverse_lower_val == 0.0 {
                    zero_tot += 1;
                } else if reverse_lower_val == 1.0 {
                    one_tot += 1;
                }

                if reverse_upper_val == 0.0 {
                    zero_tot += 1;
                } else if reverse_upper_val == 1.0 {
                    one_tot += 1;
                }
            }

            let prob: f64 = 1.0/(1.0+(-expo).exp());

            if prob > 0.5 {
                final_sequence.push(1);
            } else {
                final_sequence.push(0);
            }

            // if one_tot >= zero_tot {
            //     final_sequence.push(1);
            // } else {
            //     final_sequence.push(0);
            // }


        }
    }

    return final_sequence;
}

pub fn new_impute_single(pbwt_data: &DualPbwt, test_sequence: &Vec<u8>) -> Vec<u8> {

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

    let forward_places = spaced_insert_place(&pbwt_data.forward_pbwt, test_sequence);
    let mut reverse_places = spaced_insert_place(&pbwt_data.reverse_pbwt, &test_rev);
    reverse_places.reverse();

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

            // println!("{} and {} and {} and {}",i,one_freq,loc,buckets[loc]);
            let bucket_midpoint: f64;

            if loc == 0 {
                bucket_midpoint = buckets[0]/2.0;
            } else {
                bucket_midpoint = (buckets[loc]+buckets[loc-1])/2.0;
            }

            let mut made_coeffs: Vec<f64> = vec![0.0; side_distance as usize];

            let place_forward = forward_places[i];
            let place_reverse = reverse_places[i];

            let forward_min = max(place_forward-side_distance+1,0) as usize;
            let forward_max = min(place_forward+side_distance+1,panel_size as i32) as usize;
            let reverse_min = max(place_reverse-side_distance+1,0) as usize;
            let reverse_max = min(place_reverse+side_distance+1,panel_size as i32) as usize;


            let forward_vals_lower_ref = &pbwt_data.forward_pbwt.bin_pbwt[i][forward_min..(place_forward+1) as usize];
            let forward_vals_upper_ref = &pbwt_data.forward_pbwt.bin_pbwt[i][(place_forward+1) as usize..forward_max];

            let reverse_vals_lower_ref = &pbwt_data.reverse_pbwt.bin_pbwt[m-i-1][reverse_min..(place_reverse+1) as usize];
            let reverse_vals_upper_ref = &pbwt_data.reverse_pbwt.bin_pbwt[m-i-1][(place_reverse+1) as usize..reverse_max];

            let mut forward_vals_lower = Vec::new();
            let mut forward_vals_upper = Vec::new();

            let mut reverse_vals_lower = Vec::new();
            let mut reverse_vals_upper = Vec::new();

            for i in 0..forward_vals_lower_ref.len() {
                forward_vals_lower.push(forward_vals_lower_ref[forward_vals_lower_ref.len()-1-i]);
            }

            for i in 0..forward_vals_upper_ref.len() {
                forward_vals_upper.push(forward_vals_upper_ref[i]);
            }

            for i in 0..reverse_vals_lower_ref.len() {
                reverse_vals_lower.push(reverse_vals_lower_ref[reverse_vals_lower_ref.len()-1-i]);
            }

            for i in 0..reverse_vals_upper_ref.len() {
                reverse_vals_upper.push(reverse_vals_upper_ref[i]);
            }




            let mut zero_tot = 0;
            let mut one_tot = 0;

            let mut expo: f64 = 0.0;

            let mut pretransform = vec![0; side_distance as usize];


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
            // let prob: f64 = 1.0/(1.0+(-expo).exp());

            if prob > prob_cutoff[loc] {
                final_sequence.push(1);
            } else {
                final_sequence.push(0);
            }

        }
    }

    return final_sequence;
}

pub fn impute(pbwt_data : Arc<pbwt_structs::DualPbwt>, test_sequences : Vec<Vec<u8>>,num_threads: usize) -> Vec<Vec<u8>> {

    let N = test_sequences.len();

    let mut fin_imputed: Vec<Vec<u8>> = vec![Vec::new();N];
    let safe_pbwt = Arc::new(pbwt_data);

    let (tx,rx) = mpsc::channel();

    let pool = ThreadPool::new(num_threads);

    let mut i = 0;
    for tst_seq in test_sequences.into_iter(){

        let txr = tx.clone();
        let safe_test_seq = Arc::new(tst_seq);
        let refer = Arc::clone(&safe_pbwt);

        let clo = closure!(move txr, move refer, move safe_test_seq, || {

            let val = new_impute_single(&refer,&safe_test_seq);

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
