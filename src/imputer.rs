// use crate::pbwt;
// use std::cmp;

// fn get_score(first : u8,second : u8,p : f64) -> f64 {
//     if first == second {
//         return p.ln();
//     } else {
//         return (1.0-p).ln();
//     }
// }

// pub fn impute(haplotypes: &Vec<Vec<u8>>, pbwt_data : &pbwt::PbwtInfo, test_sequence : &Vec<u8>) -> Vec<u8> {
//     let p: f64 = 0.999;
//     let alpha: f64 = 0.9;
    
//     let mut final_sequence: Vec<u8> = Vec::new();
//     let mut position : i32 = (pbwt_data.pbwt_data[0].len()-1) as i32;

//     let min_pos: i32 = 0;
//     let max_pos: i32 = (pbwt_data.pbwt_data[0].len()) as i32;

//     let mut cur_value: u8 = 0;
//     for i in 0..test_sequence.len() {
//         cur_value = test_sequence[i];


//         if cur_value != 255 {
//             final_sequence.push(cur_value);
//             position = pbwt::get_position(pbwt_data,i,position,cur_value);
//         } else {
//             let lowest = cmp::max(min_pos,position-4) as usize;
//             let highest = cmp::min(max_pos,position+6) as usize;

//             let mut zero_score : f64 = 0.0;
//             let mut one_score : f64 = 0.0;

//             let mut total_zeros = 0;
//             let mut total_ones = 0;

//             let low_rows = &pbwt_data.pbwt_data[i][lowest..(position+1) as usize];
//             let high_rows = &pbwt_data.pbwt_data[i][(position+1) as usize..highest];

//             let low_rows_len = low_rows.len();
            
//             for j in 0..low_rows.len() {
//                 let cur_row = low_rows[j];
//                 let cur_row_val = haplotypes[cur_row as usize][i];

//                 if cur_row_val == 0 {
//                     total_zeros += 1;
//                 } else {
//                     total_ones += 1;
//                 }


//                 let zero_add = get_score(0,cur_row_val,p);
//                 let one_add = get_score(1,cur_row_val,p);

//                 zero_score += alpha.powf((low_rows_len-1-j) as f64)*get_score(0,cur_row_val,p);
//                 one_score += alpha.powf((low_rows_len-1-j) as f64)*get_score(1,cur_row_val,p);
//             }

//             for j in 0..high_rows.len() {
//                 let cur_row = high_rows[j];
//                 let cur_row_val = haplotypes[cur_row as usize][i];

//                 if cur_row_val == 0 {
//                     total_zeros += 1;
//                 } else {
//                     total_ones += 1;
//                 }

//                 let zero_add = get_score(0,cur_row_val,p);
//                 let one_add = get_score(1,cur_row_val,p);

//                 zero_score += alpha.powf(j as f64)*get_score(0,cur_row_val,p);
//                 one_score += alpha.powf(j as f64)*get_score(1,cur_row_val,p);
//             }

//             //if zero_score > one_score {
//             //    cur_value = 0;
//             //} else {
//             //    cur_value = 1;
//             //}

//             if total_zeros >= total_ones {
//                 cur_value = 0;
//             } else {
//                 cur_value = 1;
//             }

//             final_sequence.push(cur_value);
//             position = pbwt::get_position(pbwt_data,i,position,cur_value);
//         }
//     }

//     return final_sequence;

// }