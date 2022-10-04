pub mod fasta_loader;

pub mod pbwt;


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

    println!("{:?}",pbwt::pbwt(&a, 1));

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




