use std::collections::HashMap;
use crate::vcf_structs::VCFData;
use serde::{Serialize,Deserialize};
use crate::pbwt_structs::PbwtInfo;

pub fn pbwt(haplotypes : &Vec<Vec<u8>>, fm_gap : u32) -> PbwtInfo{


    let m = haplotypes.len();

    let n = haplotypes[0].len();

    let mut prefixes : Vec<Vec<u32>> = Vec::with_capacity(n+1);
    let mut divergences : Vec<Vec<u32>> = Vec::with_capacity(n+1);
    let mut binaries: Vec<Vec<u8>> = Vec::with_capacity(n+1);
    let mut cur_prefix : Vec<u32> = Vec::from_iter(0..m as u32);
    let mut cur_divergence : Vec<u32> = vec![0; m];

    let mut count_vec: Vec<u32> = Vec::new();
    let mut occ_vec : Vec<Vec<Vec<u32>>> = Vec::new();

    prefixes.push(cur_prefix.clone());
    divergences.push(cur_divergence.clone());

    for i in 0..n as u32 {
        let mut a: Vec<u32> = Vec::new();
        let mut b: Vec<u32> = Vec::new();
        let mut d: Vec<u32> = Vec::new();
        let mut e: Vec<u32> = Vec::new();
        let mut p: u32 = i+1;
        let mut q: u32 = i+1;

        let mut zero_count_val: u32 = 0;
        let basemap : HashMap<i32,u32> = [(-1,0)].iter().cloned().collect();
        let mut occ_positions: Vec<Vec<u32>> = vec![Vec::new(),Vec::new()];

        let mut ct: i32 = 0;
        let mut zero_tot = 0;
        let mut one_tot = 0;

        for (idx,start_point) in
         cur_prefix.iter().zip(cur_divergence.iter()) {
            let cur_allele: u8 = haplotypes[*idx as usize][i as usize];
            
            if *start_point > p {
                p = *start_point;
            }

            if *start_point > q {
                q = *start_point;
            }

            if cur_allele == 0 {
                a.push(*idx);
                d.push(p);
                p = 0;

                zero_count_val += 1;
                zero_tot += 1;
            }

            if cur_allele == 1 {
                b.push(*idx);
                e.push(q);
                q = 0;

                one_tot += 1;
            }

            if ct.rem_euclid(fm_gap as i32) == 0 {
                occ_positions[0].push(zero_tot);
                occ_positions[1].push(one_tot);

            }
            ct += 1;

        }

        let mut new_prefix = a.clone();
        new_prefix.append(&mut b.clone());
        let mut new_divergence = d.clone();
        new_divergence.append(&mut e.clone());

        let mut bin_values: Vec<u8> = Vec::with_capacity(new_prefix.len());


        for v in (0..m) {
            bin_values.push(haplotypes[cur_prefix[v] as usize][i as usize]);
        }

        prefixes.push(new_prefix.clone());
        divergences.push(new_divergence.clone());
        binaries.push(bin_values);

        cur_prefix = new_prefix;
        cur_divergence = new_divergence;

        count_vec.push(zero_count_val);
        occ_vec.push(occ_positions);
        
    }

    return PbwtInfo {
        num_samples: binaries[0].len() as u32,
        num_sites: binaries.len() as u32,
        bin_pbwt : binaries,
        //pbwt_data : prefixes,
        //divergence_array : divergences,
        count : count_vec,
        //occ_list : occ_vec,
        occ_list : occ_vec,
        fm_gap : fm_gap,
    };
}

pub fn get_pbwt_vcf(vcf_data: &VCFData,fm_gap: u32) -> PbwtInfo {
    return pbwt(&vcf_data.vcf_data, fm_gap);
}

pub fn get_position(pbwt_data : &PbwtInfo,
     i: usize, location: i32, val: u8) -> i32 {

        let mut final_position : i32 = -2;

        if val == 0 {
            //let occ_index = &pbwt_data.occ_list[i][0];
            let new_occ_index = &pbwt_data.occ_list[i][0];

            let fm = pbwt_data.fm_gap;

            let mut tot_add = 0;
            let mut cur_loc = location;

            while (cur_loc.rem_euclid(fm as i32) != 0) && (cur_loc != -1) {
                let cur_val = pbwt_data.bin_pbwt[i][cur_loc as usize];
                if cur_val == 0 {
                    tot_add += 1;
                }
                cur_loc -= 1;
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

            while (cur_loc.rem_euclid(fm as i32) != 0) && (cur_loc != -1) {
                let cur_val = pbwt_data.bin_pbwt[i][cur_loc as usize];
                if cur_val == 1 {
                    tot_add += 1;
                }
                cur_loc -= 1;
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
     }

pub fn insert_place(pbwt_data: &PbwtInfo,
     test_sequence: &Vec<u8>) -> Vec<i32> {
        let mut insert_positions : Vec<i32> = Vec::with_capacity(test_sequence.len()+1);
        insert_positions.push((pbwt_data.bin_pbwt[0].len()-1) as i32);
        for i in 0..test_sequence.len() {
            let cur_pos = insert_positions[insert_positions.len()-1];
            let cur_val = test_sequence[i];

            let next_pos = get_position(pbwt_data,i,cur_pos,cur_val);
            insert_positions.push(next_pos);
        }


        return insert_positions;
     }

pub fn recover_sequence(pbwt_data: &PbwtInfo, index: u32) -> Vec<u8> {
    let length = pbwt_data.bin_pbwt.len();
    let mut fin: Vec<u8> = Vec::with_capacity(length);
    let mut cur_loc = index;
    for i in (0..length) {
        let val = pbwt_data.bin_pbwt[i][cur_loc as usize];
        fin.push(val);
        cur_loc = get_position(pbwt_data,i,cur_loc as i32,val) as u32;
    }
    return fin;
}


