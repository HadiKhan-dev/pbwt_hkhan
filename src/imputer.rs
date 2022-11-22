use crate::pbwt;
use crate::pbwt_structs;
use std::cmp;
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