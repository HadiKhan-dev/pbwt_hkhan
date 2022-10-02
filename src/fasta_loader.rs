use std::cmp;
use bio::io::fasta;
use std::collections::HashMap;

pub fn get_variations(files_list : &Vec<&str>,zero_pad: bool) -> Vec<Vec<u8>> {
    let mut final_vec: Vec<Vec<u8>> = Vec::new();
    let mut genome_strings: Vec<String> = Vec::new();
    let mut max_len = 0;

    for filename in files_list {
        let reader = fasta::Reader::from_file(filename).unwrap();
        
        for result in reader.records() {
            let record = result.expect("FASTA parse error");
            let val_string = String::from_utf8_lossy(record.seq()).to_string();
            max_len = cmp::max(max_len,val_string.len());
            genome_strings.push(val_string);
        }
    
    };
    let mut hash_array: Vec<HashMap<char,u8>> = Vec::new();

    for _i in 0..max_len {
        hash_array.push(HashMap::new());
    }


    for strin in genome_strings {
        let mut cur_vec : Vec<u8> = Vec::new();

        let mut idx = 0;
        for ch in strin.chars() {
            let upch = ch.to_ascii_uppercase();
            if upch == '-' || upch == 'N' {
                cur_vec.push(255);
            } else if hash_array[idx].is_empty() {
                hash_array[idx].insert(upch,0);
                cur_vec.push(0);
            } else if hash_array[idx].contains_key(&upch) {
                cur_vec.push(*hash_array[idx].get(&upch).unwrap());
            } else {
                let mut fin_max: u8 = 0;
                for (_key,value) in &hash_array[idx] {
                    if value > &fin_max {
                        fin_max = *value;
                    }
                }
                fin_max += 1;
                hash_array[idx].insert(upch,fin_max);
                cur_vec.push(fin_max);
            }

            idx += 1;
        }
        if zero_pad {
            while idx < max_len {
                cur_vec.push(0);
                idx += 1;
            }
        }

        final_vec.push(cur_vec);
    }
    return final_vec;
}

pub fn bool_leveling(unleveled: &Vec<Vec<u8>>) -> Vec<Vec<u8>> {
    let mut new_vec: Vec<Vec<u8>> = Vec::new();
    for i in 0..unleveled.len() {
        new_vec.push(Vec::new());
        for j in 0..unleveled[i].len() {
            if unleveled[i][j] > 1 && unleveled[i][j] < 255 {
                new_vec[i].push(1);
            } else {
                new_vec[i].push(unleveled[i][j]);
            }
        }
    }
    return new_vec;
}

// Needs to have a square vector of vectors passed in, i.e. basic[i]
// must have the same length have the same length for all i
pub fn remove_all_zeros(basic: &Vec<Vec<u8>>) -> Vec<Vec<u8>> {
    let mut new_vec: Vec<Vec<u8>> = Vec::new();
    for _i in 0..basic.len() {
        new_vec.push(Vec::new());
    }
    for i in 0..basic[0].len() {
        let mut found : bool = false;
        for j in 0..basic.len() {
            if basic[j][i] != basic[0][i] {
                found = true;
                break;
            }
        }
        if found {
            for j in 0..basic.len() {
                new_vec[j].push(basic[j][i]);
            }
        }
    }
    return new_vec;
}