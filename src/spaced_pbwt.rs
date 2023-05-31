use threadpool::ThreadPool;

use crate::vcf_loader;
use crate::vcf_structs::{VCFData, SiteRow};
use crate::pbwt_structs::{SpacedPbwt, DualPbwt};
use closure::closure;
use std::collections::HashSet;
use std::sync::{Arc, mpsc};

pub fn spaced_pbwt(vcf: &VCFData, pbwt_cols: &Vec<SiteRow>, fm_gap: u32) -> SpacedPbwt {

    let data_positions: Vec<u32> = vcf.positions.clone();
    let mut pbwt_positions: Vec<u32> = Vec::new();
    let mut insert_positions: Vec<u32> = Vec::new();
    let data: &Vec<Vec<u8>> = &vcf.vcf_data;

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

    println!("PBWT: {},{},{}",n_full,n_other,n);

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

    println!("Size: {},{}",binaries.len(),binaries[0].len());


    return SpacedPbwt {
        num_samples: m as u32,
        num_pbwt_sites: n as u32,
        num_inserted_sites: n_other as u32,
        num_total_sites: n_full as u32,

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


pub fn spaced_get_position(pbwt_data : &SpacedPbwt,
    is_pbwt_col: u8, i: usize, location: i32, val: u8) -> i32 {
    /**
    Function to get the position a sequence moves to based on its current position
    and upcoming value in a SpacedPBWT
    */

    if is_pbwt_col == 1 {

        if (val != 0) && (val != 1) {
            panic!("Trying to update PBWT position based on missing data");
        }

        let mut final_position : i32 = -2;

        if val == 0 {
            //let occ_index = &pbwt_data.occ_list[i][0];
            let new_occ_index = &pbwt_data.occ_list[i][0];

            let fm = pbwt_data.fm_gap;

            let mut tot_add = 0;
            let mut cur_loc = location;

            let mut rem_diff = cur_loc.rem_euclid(fm as i32);

            while (rem_diff != 0) && (cur_loc != -1) {
               let cur_val = pbwt_data.bin_pbwt[i][cur_loc as usize];
               if cur_val == 0 {
                    tot_add += 1;
               }
               cur_loc -= 1;
               rem_diff -= 1;
            }
            let extra_add = {
               if cur_loc == -1 {
                    0
               } else if cur_loc.rem_euclid(fm as i32) == 0 {
                    let index = (cur_loc as u32)/fm;
                    new_occ_index[index as usize]
               } else {
                    panic!("Something bad happened");
                    0
               }
            };

            final_position = (extra_add as i32)+(tot_add as i32)-1;

        }

        if val == 1 {

            //let occ_index = &pbwt_data.occ_list[i][1];
            let new_occ_index = &pbwt_data.occ_list[i][1];

            let fm = pbwt_data.fm_gap;

            let mut tot_add = 0;
            let mut cur_loc = location;

            let mut rem_diff = cur_loc.rem_euclid(fm as i32);

            while (rem_diff != 0) && (cur_loc != -1) {
                let cur_val = pbwt_data.bin_pbwt[i][cur_loc as usize];
                if cur_val == 1 {
                    tot_add += 1;
                }
                cur_loc -= 1;
                rem_diff -= 1;
            }
            let extra_add = {
                if cur_loc == -1 {
                    0
                } else if cur_loc.rem_euclid(fm as i32) == 0 {
                    let index = (cur_loc as u32)/fm;
                    new_occ_index[index as usize]
                } else {
                    panic!("Something bad happened");
                    0
                }
            };

            final_position = (extra_add as i32)+(tot_add as i32)+(pbwt_data.count[i as usize] as i32)-1;
        }
        return final_position;
    } else {
        return location;
    }
}

pub fn spaced_insert_place(pbwt_data: &SpacedPbwt,
    test_sequence: &Vec<u8>) -> Vec<i32> {
        let mut insert_positions : Vec<i32> = Vec::with_capacity(test_sequence.len()+1);
        insert_positions.push((pbwt_data.bin_pbwt[0].len()-1) as i32);
        for i in 0..test_sequence.len() {
            let col_type = pbwt_data.pbwt_col_flags[i];
            let cur_pos = insert_positions[insert_positions.len()-1];
            let cur_val = test_sequence[i];

            let next_pos = spaced_get_position(pbwt_data,col_type,i,cur_pos,cur_val);
            insert_positions.push(next_pos);
       }


       return insert_positions;
    }

pub fn spaced_recover_sequence(pbwt_data: &SpacedPbwt, index: u32) -> Vec<u8> {
    let length = pbwt_data.bin_pbwt.len();
    let mut fin: Vec<u8> = Vec::with_capacity(length);
    let mut cur_loc = index;
    for i in (0..length) {
        let val = pbwt_data.bin_pbwt[i][cur_loc as usize];
        fin.push(val);
        cur_loc = spaced_get_position(pbwt_data,pbwt_data.pbwt_col_flags[i],i,cur_loc as i32,val) as u32;
    }
    return fin;
}

pub fn dual_pbwt(vcf_data: Arc<VCFData>, pbwt_cols: Arc<Vec<SiteRow>>, fm_gap: u32) -> DualPbwt{


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
