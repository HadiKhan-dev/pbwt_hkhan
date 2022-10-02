use ndarray::{array};

pub mod fasta_loader;


fn pbwt(haplotypes : &Vec<Vec<u8>>) -> ndarray::Array1::<i32>{
    let m = haplotypes.len();

    let n = haplotypes[0].len();

    let mut prefixes : Vec<Vec<u32>> = Vec::with_capacity(n+1);
    let mut divergences : Vec<Vec<u32>> = Vec::with_capacity(n+1);
    
    let mut cur_prefix : Vec<u32> = Vec::from_iter(0..m as u32);
    let mut cur_divergence : Vec<u32> = vec![0; m];

    prefixes.push(cur_prefix.clone());
    divergences.push(cur_divergence.clone());

    for i in 0..n as u32 {
        let mut a: Vec<u32> = Vec::new();
        let mut b: Vec<u32> = Vec::new();
        let mut d: Vec<u32> = Vec::new();
        let mut e: Vec<u32> = Vec::new();
        let mut p: u32 = i+1;
        let mut q: u32 = i+1;

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
            }

            if cur_allele == 1 {
                b.push(*idx);
                e.push(q);
                q = 0;
            }

        }

        let mut new_prefix = a.clone();
        new_prefix.append(&mut b.clone());
        let mut new_divergence = d.clone();
        new_divergence.append(&mut e.clone());

        prefixes.push(new_prefix.clone());
        divergences.push(new_divergence.clone());

        cur_prefix = new_prefix;
        cur_divergence = new_divergence;
        
    }

    println!("{:?}",prefixes);

    return array![1,2,3];
}
fn main() {
    let a = vec![
    vec![0, 1, 0, 1, 0, 1],
    vec![1, 1, 0, 0, 0, 1],
    vec![1, 1, 1, 1, 1, 1],
    vec![0, 1, 1, 1, 1, 0],
    vec![0, 0, 0, 0, 0, 0],
    vec![1, 0, 0, 0, 1, 0],
    vec![1, 1, 0, 0, 0, 1],
    vec![0, 1, 0, 1, 1, 0]];

    println!("{}",pbwt(&a));

    //let files_list = vec!["sequence.fasta", "sequence (1).fasta","sequence (2).fasta"];
    //let variation_vals = fasta_loader::get_variations(&files_list, true);
    //let filled = fasta_loader::bool_leveling(&variation_vals);
    //let fixed = fasta_loader::remove_all_zeros(&filled);

    //println!("");
    //for item in &fixed {
    //    println!("{:?}",item.get(..30).unwrap());
    //}
    //println!("{:?}",variation_vals);
}




