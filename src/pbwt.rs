use std::collections::HashMap;
use crate::vcf_structs::VCFData;
use serde::{Serialize,Deserialize};


#[derive(Debug,Serialize,Deserialize)]
pub struct PbwtInfo {
    pub pbwt_data : Vec<Vec<u32>>,
    pub divergence_array : Vec<Vec<u32>>,
    pub count : Vec<u32>,
    pub occ_list : Vec<Vec<HashMap<i32,u32>>>,
}

pub fn pbwt(haplotypes : &Vec<Vec<u8>>, fm_gap : u32) -> PbwtInfo{


    let m = haplotypes.len();

    let n = haplotypes[0].len();

    let mut prefixes : Vec<Vec<u32>> = Vec::with_capacity(n+1);
    let mut divergences : Vec<Vec<u32>> = Vec::with_capacity(n+1);
    
    let mut cur_prefix : Vec<u32> = Vec::from_iter(0..m as u32);
    let mut cur_divergence : Vec<u32> = vec![0; m];

    let mut count_vec: Vec<u32> = Vec::new();
    let mut occ_vec : Vec<Vec<HashMap<i32,u32>>> = Vec::new();

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
        let mut occ_positions: Vec<HashMap<i32,u32>> = vec![basemap.clone(),basemap.clone()];

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
                occ_positions[0].insert(ct,zero_tot);
                occ_positions[1].insert(ct,one_tot);
            }
            ct += 1;

        }

        let mut new_prefix = a.clone();
        new_prefix.append(&mut b.clone());
        let mut new_divergence = d.clone();
        new_divergence.append(&mut e.clone());

        prefixes.push(new_prefix.clone());
        divergences.push(new_divergence.clone());

        cur_prefix = new_prefix;
        cur_divergence = new_divergence;

        count_vec.push(zero_count_val);
        occ_vec.push(occ_positions);
        
    }

    return PbwtInfo {
        pbwt_data : prefixes,
        divergence_array : divergences,
        count : count_vec,
        occ_list : occ_vec,
    };
}

pub fn get_pbwt_vcf(vcf_data: &VCFData,fm_gap: u32) -> PbwtInfo {
    return pbwt(&vcf_data.vcf_data, fm_gap);
}

pub fn get_position(haplotypes: &Vec<Vec<u8>>,pbwt_data : &PbwtInfo,
     i: usize, location: i32, val: u8) -> i32 {

        let mut final_position : i32 = -2;

        if val == 0 {
            let occ_index = &pbwt_data.occ_list[i][0];
            let mut tot_add = 0;
            let mut cur_loc = location;
            while !occ_index.contains_key(&cur_loc) {
                let row_pos = *(&pbwt_data.pbwt_data[i][cur_loc as usize]);
                let cur_val = haplotypes[row_pos as usize][i];
                if cur_val == 0 {
                    tot_add += 1;
                }
                cur_loc -= 1;
            }
            final_position = (tot_add as i32)+ (occ_index[&(cur_loc as i32)] as i32)-1;
        }

        if val == 1 {
            let occ_index = &pbwt_data.occ_list[i][1];
            let mut tot_add = 0;
            let mut cur_loc = location;
            while !occ_index.contains_key(&cur_loc) {
                let row_pos = *(&pbwt_data.pbwt_data[i][cur_loc as usize]);
                let cur_val = haplotypes[row_pos as usize][i];
                if cur_val == 1 {
                    tot_add += 1;
                }
                cur_loc -= 1;
            }
            final_position = (tot_add as i32) + (occ_index[&(cur_loc as i32)] as i32) +(pbwt_data.count[i as usize] as i32)-1;
        }
        return final_position;
     }

pub fn insert_place(haplotypes: &Vec<Vec<u8>>, pbwt_data: &PbwtInfo,
     test_sequence: &Vec<u8>) -> Vec<i32> {
        let mut insert_positions : Vec<i32> = Vec::with_capacity(test_sequence.len()+1);
        insert_positions.push((haplotypes.len()-1) as i32);
        for i in 0..test_sequence.len() {
            let cur_pos = insert_positions[insert_positions.len()-1];
            let cur_val = test_sequence[i];

            let next_pos = get_position(haplotypes,pbwt_data,i,cur_pos,cur_val);
            insert_positions.push(next_pos);
        }


        return insert_positions;
     }

