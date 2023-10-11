use threadpool::ThreadPool;

use crate::vcf_loader;
use crate::vcf_structs::{VCFData, SiteRow};
use crate::pbwt_structs::{SpacedPbwt, DualPbwt};
use crate::helper_structs::{PositionData,SideData,FullInsertData};
use closure::closure;
use std::collections::HashSet;
use std::sync::{Arc, mpsc};
use std::cmp;

pub fn spaced_pbwt(vcf: &VCFData, pbwt_cols: &Vec<SiteRow>, fm_gap: u32) -> SpacedPbwt {

    let data_positions: Vec<u32> = vcf.positions.clone();
    let mut pbwt_positions: Vec<u32> = Vec::new();
    let mut insert_positions: Vec<u32> = Vec::new();
    let data: &Vec<Vec<u8>> = &vcf.vcf_data;

    let mut allele_freqs: Vec<f32> = Vec::new();

    let mut col_set: HashSet<u32> = HashSet::new();

    let mut n: usize = 0;

    for col in pbwt_cols {
        let pos = col.position;
        col_set.insert(pos);
        n += 1;
    }
    

    let m = data.len();
    let n_full = data[0].len();

    let n_other = n_full-n;


    let mut is_pbwt_col :Vec<u8> = Vec::with_capacity(n_full+1);
    let mut pbwt_positions: Vec<u32> = Vec::new();
    let mut inserted_positions: Vec<u32> = Vec::new();
    let mut prefixes : Vec<Vec<u32>> = Vec::with_capacity(n+1);
    let mut divergences : Vec<Vec<u32>> = Vec::with_capacity(n+1);
    let mut binaries: Vec<Vec<u8>> = Vec::with_capacity(n_full+1);


    let cur_prefix : Vec<u32> = Vec::from_iter(0..m as u32);
    let cur_divergence : Vec<u32> = vec![0; m];
    let mut j: usize = 0;
    let mut j_pbwt = 0;

    let mut count_vec: Vec<u32> = Vec::new();
    let mut occ_vec : Vec<Vec<Vec<u32>>> = Vec::new();

    prefixes.push(cur_prefix);
    divergences.push(cur_divergence);

    let mut cur_prefix_ref: &Vec<u32> = &(prefixes[prefixes.len()-1]);
    let mut cur_divergence_ref: &Vec<u32> = &divergences[divergences.len()-1];

    let mut ct: i32 = 0;
    let mut ct_extra: u32 = 0;
    let mut zero_tot: u32 = 0;
    let mut one_tot: u32 = 0;
    let mut occ_positions: Vec<Vec<u32>> = vec![Vec::new(),Vec::new()];
    let mut new_add: Vec<u8> = Vec::with_capacity(m);

    let mut a: Vec<u32> = Vec::with_capacity(m);
    let mut b: Vec<u32> = Vec::with_capacity(m);
    let mut d: Vec<u32> = Vec::with_capacity(m);
    let mut e: Vec<u32> = Vec::with_capacity(m);

    let mut bin_values: Vec<u8> = Vec::with_capacity(m);


    for col in &vcf.positions {
        if !col_set.contains(&col) {


            ct = 0;
            ct_extra = 0;
            zero_tot = 0;
            one_tot = 0;
            occ_positions = vec![Vec::new(),Vec::new()];

            new_add = Vec::with_capacity(m);

            let mut val: u8;
            for idx in cur_prefix_ref.iter() {


                unsafe {
                val = *data.get_unchecked(*idx as usize).get_unchecked(j);
                }


                new_add.push(val);

                if val == 0 {
                    zero_tot += 1;
                } else if val == 1 {
                    one_tot += 1;
                }

                if (ct == 0) || (ct_extra == fm_gap) {
                    occ_positions[0].push(zero_tot);
                    occ_positions[1].push(one_tot);
                    ct_extra = 0;
    
                }

                ct += 1;
                ct_extra += 1;


            }

            binaries.push(new_add);
            is_pbwt_col.push(0);
            inserted_positions.push(*col);
            count_vec.push(zero_tot);
            occ_vec.push(occ_positions);


        } else {


            a = Vec::with_capacity(m);
            b = Vec::with_capacity(m);
            d = Vec::with_capacity(m);
            e = Vec::with_capacity(m);

            bin_values = Vec::with_capacity(m);

            let mut p: u32 = j_pbwt+1;
            let mut q: u32 = j_pbwt+1;

            occ_positions = vec![Vec::new(),Vec::new()];

            ct = 0;
            ct_extra = 0;
            zero_tot = 0;
            one_tot = 0;

            let mut cur_allele: u8;
            for (idx,start_point) in
            cur_prefix_ref.iter().zip(cur_divergence_ref.iter()) {

                let idx_val = *idx;

                unsafe {
                cur_allele = *data.get_unchecked(idx_val as usize).get_unchecked(j);
                }

                bin_values.push(cur_allele);

                let st = *start_point;

                if st > p {
                p = st;
                }

                if st > q {
                q = st;
                }

                if cur_allele == 0 {
                    a.push(idx_val);
                    d.push(p);
                    p = 0;

                    zero_tot += 1;
                }
    
                if cur_allele == 1 {
                    b.push(idx_val);
                    e.push(q);
                    q = 0;
    
                    one_tot += 1;
                }

                if (ct == 0) || (ct_extra == fm_gap) {
                    occ_positions[0].push(zero_tot);
                    occ_positions[1].push(one_tot);
                    ct_extra = 0;
    
                }
                ct += 1;
                ct_extra += 1;
            
            }  

            
            let mut new_prefix = a;
            new_prefix.append(&mut b);
            let mut new_divergence = d;
            new_divergence.append(&mut e);


            prefixes.push(new_prefix);
            divergences.push(new_divergence);
            binaries.push(bin_values);

            cur_prefix_ref = &(prefixes[prefixes.len()-1]);
            cur_divergence_ref = &divergences[divergences.len()-1];

            count_vec.push(zero_tot);
            occ_vec.push(occ_positions);

            is_pbwt_col.push(1);
            pbwt_positions.push(*col);

        }
        j += 1;

    }

    for i in 0..n_full {
        let one_freq = (1.0-(count_vec[i] as f32)/(m as f32));

        allele_freqs.push(one_freq);
    
    }


    return SpacedPbwt {
        num_samples: m as u32,
        num_pbwt_sites: n as u32,
        num_inserted_sites: n_other as u32,
        num_total_sites: n_full as u32,

        allele_freqs: allele_freqs,
        pbwt_positions: pbwt_positions,
        inserted_positions: inserted_positions,
        all_positions: data_positions,
        pbwt_col_flags: is_pbwt_col,

        bin_pbwt: binaries,
        count: count_vec,
        occ_list: occ_vec,
        fm_gap: fm_gap,
        
    };
}


// pub fn spaced_get_position(pbwt_data : &SpacedPbwt,
//     is_pbwt_col: u8, col_number: usize, location: isize, val: u8, side_distance: usize,
//     divergence_distance: usize, current_divergence: &Vec<Vec<u32>>) -> PositionData {
//     /**
//     Function to get the position a sequence moves to based on its current position
//     and upcoming value in a SpacedPBWT, also collects the n neighbouring upper and lower
//     values as well as a representation of how long the neighbours have been moving with
//     the PBWT.

//     col_number is the index of the upcoming column
//     */

//     let mut final_position: isize = -2;

//     if (current_divergence.len() != 2) || (current_divergence[0].len() != divergence_distance) {
//         println!("Size of {}, {}, {}",current_divergence.len(),current_divergence[0].len(),divergence_distance);
        
//         panic!("Current divergence array input does not match stated size");
//     }

//     if side_distance < divergence_distance {
//         panic!("Side length must be at least equal to divergence length")
//     }

//     if is_pbwt_col == 1 {

//         if (val != 0) && (val != 1) {
//             panic!("Trying to update PBWT position based on missing data");
//         }

//         if val == 0 {
//             //let occ_index = &pbwt_data.occ_list[i][0];
//             let new_occ_index = &pbwt_data.occ_list[col_number][0];

//             let fm = pbwt_data.fm_gap;

//             let mut tot_add = 0;
//             let mut cur_loc = location;

//             let mut rem_diff = cur_loc.rem_euclid(fm as isize);

//             while (rem_diff != 0) && (cur_loc != -1) {
//                 let cur_val = pbwt_data.bin_pbwt[col_number][cur_loc as usize];
//                 if cur_val == 0 {
//                     tot_add += 1;
//                 }
//                 cur_loc -= 1;
//                 rem_diff -= 1;
//             }
//             let extra_add = {
//                 if cur_loc == -1 {
//                     0
//                 } else if cur_loc.rem_euclid(fm as isize) == 0 {
//                     let index = (cur_loc as u32)/fm;
//                     new_occ_index[index as usize]
//                 } else {
//                     panic!("Something bad happened");
//                     0
//                 }
//             };

//             final_position = (extra_add as isize)+(tot_add as isize)-1;

//         }

//         if val == 1 {

//             //let occ_index = &pbwt_data.occ_list[i][1];
//             let new_occ_index = &pbwt_data.occ_list[col_number][1];

//             let fm = pbwt_data.fm_gap;

//             let mut tot_add = 0;
//             let mut cur_loc = location;

//             let mut rem_diff = cur_loc.rem_euclid(fm as isize);

//             while (rem_diff != 0) && (cur_loc != -1) {
//                 let cur_val = pbwt_data.bin_pbwt[col_number][cur_loc as usize];
//                 if cur_val == 1 {
//                     tot_add += 1;
//                 }
//                 cur_loc -= 1;
//                 rem_diff -= 1;
//             }
//             let extra_add = {
//                 if cur_loc == -1 {
//                     0
//                 } else if cur_loc.rem_euclid(fm as isize) == 0 {
//                     let index = (cur_loc as u32)/fm;
//                     new_occ_index[index as usize]
//                 } else {
//                     panic!("Something bad happened");
//                     0
//                 }
//             };

//             final_position = (extra_add as isize)+(tot_add as isize)+(pbwt_data.count[col_number as usize] as isize)-1;
//         }
//     } else {
//         final_position = location;
//     }

//     if (side_distance == 0) && (divergence_distance == 0) {
//         return PositionData{position:final_position,side_data:vec![Vec::new(),Vec::new()],
//             divergence_data: vec![Vec::new(),Vec::new()]}
//     }

//     let mut upper_data: Vec<u8>;
//     let mut lower_data: Vec<u8>;

//     if (final_position == -1) {
//         lower_data = Vec::new();
//         upper_data = pbwt_data.bin_pbwt[col_number][0..side_distance].to_vec();

//         if side_distance > pbwt_data.num_samples as usize {
//             panic!("Side distance is bigger than number of samples in PBWT");
//         }
//     } else {

//         let low_point: usize = cmp::max(0,location-(side_distance as isize)+1) as usize;
//         let high_point: usize = cmp::min(pbwt_data.num_samples as isize ,final_position+(side_distance as isize)+1) as usize;
        
//         if {col_number == 2} {
//             println!("COL {} {} {} {}",low_point,high_point,final_position,location);
//         }
//         lower_data = pbwt_data.bin_pbwt[col_number][low_point..(final_position as usize)+1].to_vec();
//         lower_data.reverse();
//         upper_data = pbwt_data.bin_pbwt[col_number][(final_position as usize) +1..high_point].to_vec();
        
//         if {col_number == 2} {
//             println!("COLEND {:?}",upper_data);
//             println!("Final {:?}",(final_position as usize) +1..high_point);
//         }
    
//     }

//     let mut new_lower_divergence: Vec<u32>;
//     let mut new_upper_divergence: Vec<u32>;

//     if is_pbwt_col == 1 {

//         let current_lower_divergence = &current_divergence[0];
//         let current_upper_divergence = &current_divergence[1];

//         new_lower_divergence = Vec::new();
//         new_upper_divergence = Vec::new();

//         for s in 0..cmp::min(lower_data.len(),divergence_distance) {
//             if lower_data[s] == val {
//                 new_lower_divergence.push(current_lower_divergence[s]+1);
//             }
//         }

//         for s in 0..cmp::min(upper_data.len(),divergence_distance) {
//             if upper_data[s] == val {
//                 new_upper_divergence.push(current_upper_divergence[s]+1);
//             }
//         }

//         new_lower_divergence.append(&mut vec![0;divergence_distance-new_lower_divergence.len()]);
//         new_upper_divergence.append(&mut vec![0;divergence_distance-new_upper_divergence.len()]);
//     } else {
//         new_lower_divergence = current_divergence[0].clone();
//         new_upper_divergence = current_divergence[1].clone();

//         if col_number != 0 {
//             for i in (0..new_lower_divergence.len()) {
//                 new_lower_divergence[i] += 0;
//                 new_upper_divergence[i] += 0;
//             }
//         }
//     }

//     if {col_number == 2} {
//         println!("SIDE RETURN {:?}",vec![lower_data.clone(),upper_data.clone()]);

//     }

//     return PositionData{position:final_position,side_data:vec![lower_data,upper_data],
//         divergence_data: vec![new_lower_divergence,new_upper_divergence]};


// }

pub fn spaced_get_position(pbwt_data : &SpacedPbwt,
    is_pbwt_col: u8, col_number: usize, location: isize, val: u8, side_distance: usize,
    divergence_distance: usize, current_divergence: &Vec<Vec<u32>>) -> isize {
    /**
    Function to get the position a sequence moves to based on its current position
    and upcoming value in a SpacedPBWT.

    col_number is the index of the upcoming column
    */

    let mut final_position: isize = -2;

    if (current_divergence.len() != 2) || (current_divergence[0].len() != divergence_distance) {
        println!("Size of {}, {}, {}",current_divergence.len(),current_divergence[0].len(),divergence_distance);
        
        panic!("Current divergence array input does not match stated size");
    }

    if side_distance < divergence_distance {
        panic!("Side length must be at least equal to divergence length")
    }

    if is_pbwt_col == 1 {

        if (val != 0) && (val != 1) {
            panic!("Trying to update PBWT position based on missing data");
        }

        if val == 0 {
            //let occ_index = &pbwt_data.occ_list[i][0];
            let new_occ_index = &pbwt_data.occ_list[col_number][0];

            let fm = pbwt_data.fm_gap;

            let mut tot_add = 0;
            let mut cur_loc = location;

            let mut rem_diff = cur_loc.rem_euclid(fm as isize);

            while (rem_diff != 0) && (cur_loc != -1) {
                let cur_val = pbwt_data.bin_pbwt[col_number][cur_loc as usize];
                if cur_val == 0 {
                    tot_add += 1;
                }
                cur_loc -= 1;
                rem_diff -= 1;
            }
            let extra_add = {
                if cur_loc == -1 {
                    0
                } else if cur_loc.rem_euclid(fm as isize) == 0 {
                    let index = (cur_loc as u32)/fm;
                    new_occ_index[index as usize]
                } else {
                    panic!("Something bad happened");
                    0
                }
            };

            final_position = (extra_add as isize)+(tot_add as isize)-1;

        }

        if val == 1 {

            //let occ_index = &pbwt_data.occ_list[i][1];
            let new_occ_index = &pbwt_data.occ_list[col_number][1];

            let fm = pbwt_data.fm_gap;

            let mut tot_add = 0;
            let mut cur_loc = location;

            let mut rem_diff = cur_loc.rem_euclid(fm as isize);

            while (rem_diff != 0) && (cur_loc != -1) {
                let cur_val = pbwt_data.bin_pbwt[col_number][cur_loc as usize];
                if cur_val == 1 {
                    tot_add += 1;
                }
                cur_loc -= 1;
                rem_diff -= 1;
            }
            let extra_add = {
                if cur_loc == -1 {
                    0
                } else if cur_loc.rem_euclid(fm as isize) == 0 {
                    let index = (cur_loc as u32)/fm;
                    new_occ_index[index as usize]
                } else {
                    panic!("Something bad happened");
                    0
                }
            };

            final_position = (extra_add as isize)+(tot_add as isize)+(pbwt_data.count[col_number as usize] as isize)-1;
        }
    } else {
        final_position = location;
    }


    return final_position as isize;

    

}

pub fn update_side_data(pbwt_data : &SpacedPbwt,
    is_pbwt_col: u8, col_number: usize, location: isize, val: u8, side_distance: usize,
    divergence_distance: usize, current_divergence: &Vec<Vec<u32>>) -> SideData {
    /**
    Function to get the position a sequence moves to based on its current position
    and upcoming value in a SpacedPBWT, also collects the n neighbouring upper and lower
    values as well as a representation of how long the neighbours have been moving with
    the PBWT.

    col_number is the index of the upcoming column
    */

    if (current_divergence.len() != 2) || (current_divergence[0].len() != divergence_distance) {
        println!("Size of {}, {}, {}",current_divergence.len(),current_divergence[0].len(),divergence_distance);
        
        panic!("Current divergence array input does not match stated size");
    }

    if side_distance < divergence_distance {
        panic!("Side length must be at least equal to divergence length")
    }

    if (side_distance == 0) && (divergence_distance == 0) {
        return SideData{side_data:vec![Vec::new(),Vec::new()],
            divergence_data: vec![Vec::new(),Vec::new()]}
    }

    let mut upper_data: Vec<u8>;
    let mut lower_data: Vec<u8>;

    if (location == -1) {
        lower_data = Vec::new();
        upper_data = pbwt_data.bin_pbwt[col_number][0..side_distance].to_vec();

        if side_distance > pbwt_data.num_samples as usize {
            panic!("Side distance is bigger than number of samples in PBWT");
        }
    } else {

        let low_point: usize = cmp::max(0,location-(side_distance as isize)+1) as usize;
        let high_point: usize = cmp::min(pbwt_data.num_samples as isize ,location+(side_distance as isize)+1) as usize;
        
        lower_data = pbwt_data.bin_pbwt[col_number][low_point..(location as usize)+1].to_vec();
        lower_data.reverse();
        upper_data = pbwt_data.bin_pbwt[col_number][(location as usize) +1..high_point].to_vec();
        
    
    }

    let mut new_lower_divergence: Vec<u32>;
    let mut new_upper_divergence: Vec<u32>;

    if is_pbwt_col == 1 {

        let current_lower_divergence = &current_divergence[0];
        let current_upper_divergence = &current_divergence[1];

        new_lower_divergence = Vec::new();
        new_upper_divergence = Vec::new();

        for s in 0..cmp::min(lower_data.len(),divergence_distance) {
            if lower_data[s] == val {
                new_lower_divergence.push(current_lower_divergence[s]+1);
            }
        }

        for s in 0..cmp::min(upper_data.len(),divergence_distance) {
            if upper_data[s] == val {
                new_upper_divergence.push(current_upper_divergence[s]+1);
            }
        }

        new_lower_divergence.append(&mut vec![0;divergence_distance-new_lower_divergence.len()]);
        new_upper_divergence.append(&mut vec![0;divergence_distance-new_upper_divergence.len()]);
    } else {
        new_lower_divergence = current_divergence[0].clone();
        new_upper_divergence = current_divergence[1].clone();

        if col_number != 0 {
            for i in (0..new_lower_divergence.len()) {
                new_lower_divergence[i] += 0;
                new_upper_divergence[i] += 0;
            }
        }
    }


    return SideData{side_data:vec![lower_data,upper_data],
        divergence_data: vec![new_lower_divergence,new_upper_divergence]};

    

}
pub fn spaced_insert_place(pbwt_data: &SpacedPbwt,
    test_sequence: &Vec<u8>, side_distance: usize, divergence_distance: usize) -> FullInsertData {
    /**
    Function to get the positions a sequence will insert into if it is placed into a constructed
    PBWT. It also computes the side_distance nearest neighbours both up and down and how long the 
    divergence_distance nearest neighbours have been moving with the sequence.

    col_number is the index of the upcoming column
    */    
    let mut insert_positions : Vec<isize> = Vec::with_capacity(test_sequence.len()+1);
    insert_positions.push((pbwt_data.bin_pbwt[0].len()-1) as isize);
       
    let mut side_values: Vec<Vec<Vec<u8>>> = Vec::with_capacity(test_sequence.len());
    let mut divergence_values: Vec<Vec<Vec<u32>>> = Vec::with_capacity(test_sequence.len());

    let mut current_divergence:Vec<Vec<u32>> = vec![vec![0;divergence_distance],vec![0;divergence_distance]];
    divergence_values.push(current_divergence.clone());

    for i in 0..test_sequence.len() {
        let col_type = pbwt_data.pbwt_col_flags[i];
        let cur_pos = insert_positions[insert_positions.len()-1];
        let cur_val = test_sequence[i];
            
        let upcoming_position = spaced_get_position(pbwt_data,col_type,
            i,cur_pos,cur_val,side_distance,divergence_distance,
            &current_divergence);

        let side_data = update_side_data(pbwt_data,col_type,
            i,cur_pos,cur_val,side_distance,divergence_distance,
            &current_divergence);
            
        insert_positions.push(upcoming_position);


        side_values.push(side_data.side_data);
        divergence_values.push(side_data.divergence_data.clone());

        current_divergence = side_data.divergence_data;

        }


    return FullInsertData{insert_positions: insert_positions, side_values: side_values,
        divergence_values:divergence_values};
}

pub fn spaced_recover_sequence(pbwt_data: &SpacedPbwt, index: isize) -> Vec<u8> {
    /**
    Function to recover the i^th sequence of the reference panel that makes up the PBWT given as
    input.
    */    
    let length = pbwt_data.bin_pbwt.len();
    let mut fin: Vec<u8> = Vec::with_capacity(length);
    let mut cur_loc: isize = index;

    if index < 0 {
        panic!("Index of sequence to recover must be non-negative");
    }
    
    let divergence_dummy: Vec<Vec<u32>> = vec![vec![],vec![]];

    for i in (0..length) {
        let val = pbwt_data.bin_pbwt[i][cur_loc as usize];
        fin.push(val);
        let positional_data = spaced_get_position(pbwt_data,pbwt_data.pbwt_col_flags[i],
            i,cur_loc,val,0,
            0,&divergence_dummy);

        cur_loc = positional_data;
    }
    return fin;
}

pub fn dual_pbwt(vcf_data: Arc<VCFData>, pbwt_cols: Arc<Vec<SiteRow>>, fm_gap: u32) -> DualPbwt{
    /**
    Function to compute both the forward and backward PBWTs for a reference panel.
    */    

    let pool = ThreadPool::new(1);

    let (tx,rx) = mpsc::channel();

    let vcf_copy = Arc::clone(&vcf_data);
    let pbwt_col_copy = Arc::clone(&pbwt_cols);


    let first_closure = closure!(move vcf_copy, move pbwt_col_copy, move tx, || {
        let now = std::time::Instant::now();
        let reverse_vcf = vcf_loader::reverse_vcf(&vcf_copy);

        let elapsed = now.elapsed();

        println!("Reverse time: {:.4?}", elapsed);

        let reverse_pbwt = spaced_pbwt(&reverse_vcf,&pbwt_col_copy,fm_gap);

        tx.send(reverse_pbwt).unwrap();

    });

    pool.execute(first_closure);

    let forward = spaced_pbwt(&vcf_data,&pbwt_cols,fm_gap);
    let mut reverse: SpacedPbwt = rx.recv().unwrap();


    return DualPbwt { forward_pbwt: forward, reverse_pbwt: reverse };
}
