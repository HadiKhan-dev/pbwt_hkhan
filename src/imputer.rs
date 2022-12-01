use crate::pbwt;
use crate::pbwt_structs;
use crate::pbwt_structs::{SpacedPbwt, DualPbwt};
use crate::spaced_pbwt::spaced_insert_place;
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

fn get_score(first : u8,second : u8,p : f64) -> f64 {
    if first == second {
        return p.ln();
    } else {
        return (1.0-p).ln();
    }
}

pub fn impute_single(pbwt_data : &pbwt_structs::PbwtInfo, test_sequence : &Vec<u8>) -> Vec<u8> {
    let p: f64 = 0.999;
    let alpha: f64 = 0.9;

    let mut side_distance = 3;
    
    let mut final_sequence: Vec<u8> = Vec::new();
    let mut position : i32 = (pbwt_data.bin_pbwt[0].len()-1) as i32;

    let min_pos: i32 = 0;
    let max_pos: i32 = (pbwt_data.bin_pbwt[0].len()) as i32;

    let mut cur_value: u8 = 0;
    for i in 0..test_sequence.len() {
        cur_value = test_sequence[i];


        if cur_value != 255 {
            final_sequence.push(cur_value);
            position = pbwt::get_position(pbwt_data,i,position,cur_value);
        } else {
            let lowest = cmp::max(min_pos,position-2) as usize;
            let highest = cmp::min(max_pos,position+4) as usize;

            let mut zero_score : f64 = 0.0;
            let mut one_score : f64 = 0.0;

            let mut total_zeros = 0;
            let mut total_ones = 0;

            let low_rows = &pbwt_data.bin_pbwt[i][lowest..(position+1) as usize];
            let high_rows = &pbwt_data.bin_pbwt[i][(position+1) as usize..highest];

            let low_rows_len = low_rows.len();
            
            for j in 0..low_rows.len() {
                //let cur_row = low_rows[j];
                //let cur_row_val = haplotypes[cur_row as usize][i];
                let cur_row_val = low_rows[j];
                
                if cur_row_val == 0 {
                    total_zeros += 1;
                } else {
                    total_ones += 1;
                }
                let zero_add = get_score(0,cur_row_val,p);
                let one_add = get_score(1,cur_row_val,p);

                zero_score += alpha.powf((low_rows_len-1-j) as f64)*get_score(0,cur_row_val,p);
                one_score += alpha.powf((low_rows_len-1-j) as f64)*get_score(1,cur_row_val,p);
            }

            for j in 0..high_rows.len() {
                //let cur_row = high_rows[j];
                //let cur_row_val = haplotypes[cur_row as usize][i];
                let cur_row_val = high_rows[j];

                if cur_row_val == 0 {
                    total_zeros += 1;
                } else {
                    total_ones += 1;
                }

                let zero_add = get_score(0,cur_row_val,p);
                let one_add = get_score(1,cur_row_val,p);

                zero_score += alpha.powf(j as f64)*get_score(0,cur_row_val,p);
                one_score += alpha.powf(j as f64)*get_score(1,cur_row_val,p);
            }

            if zero_score > one_score {
                cur_value = 0;
            } else {
                cur_value = 1;
            }

            // if total_zeros > total_ones {
            //     cur_value = 0;
            // } else {
            //     cur_value = 1;
            // }

            final_sequence.push(cur_value);
            position = pbwt::get_position(pbwt_data,i,position,cur_value);
        }
    }

    return final_sequence;
}

pub fn new_impute_single(pbwt_data: &DualPbwt, test_sequence: &Vec<u8>) -> Vec<u8> {

    let side_distance = 4;
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

            let intercept = -7.96707151;

            let coeffs = vec![16.39305966,0.97798448];

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

pub fn impute(pbwt_data : Arc<pbwt_structs::PbwtInfo>, test_sequences : Vec<Vec<u8>>,num_threads: usize) -> Vec<Vec<u8>> {

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

            let val = impute_single(&refer,&safe_test_seq);

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


pub fn new_impute(pbwt_data : Arc<pbwt_structs::DualPbwt>, test_sequences : Vec<Vec<u8>>,num_threads: usize) -> Vec<Vec<u8>> {

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