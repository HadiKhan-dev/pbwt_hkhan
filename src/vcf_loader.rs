use vcf::{VCFReader, U8Vec, VCFHeaderFilterAlt, VCFError, VCFRecord};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::BufReader;


#[derive(Debug)]
pub struct VCFData {
    pub vcf_data: Vec<Vec<u8>>,
    pub positions: Vec<u64>,
    pub sample_names: Vec<String>,
    pub haplotype_names: Vec<String>,
}

pub fn read(filename: &str) -> Result<VCFData, VCFError> {
    let mut reader = VCFReader::new(BufReader::new(MultiGzDecoder::new(File::open(
        filename,
    )?)))?;

    let mut sample_names = Vec::new();
    let mut haplotype_names = Vec::new();

    //println!("{:?}",reader.header());


    for i in reader.header().samples() {
        //let j = i.line();
        let cur_sample = String::from_utf8(i.to_vec()).unwrap();

        let mut first = cur_sample.clone();
        let mut second = cur_sample.clone();

        first.push_str("[0]");
        second.push_str("[1]");

        sample_names.push(cur_sample);
        haplotype_names.push(first);
        haplotype_names.push(second);
    }

    let mut panel_matrix: Vec<Vec<u8>> = Vec::with_capacity(haplotype_names.len());

    for i in 0..haplotype_names.len() {
        panel_matrix.push(Vec::new());
    }

    let mut vcf_record = reader.empty_record();

    let mut last_position : u64 = 0;
    let mut cur_position : u64 = 0;

    let mut started = false;
    let mut finished = false;

    let mut position_list = Vec::new();

    while !started |!finished {

        last_position = cur_position;


        reader.next_record(&mut vcf_record);
        cur_position = vcf_record.position.clone();

        //println!("{:?}",vcf_record.genotype);
        
        if started & (cur_position == last_position) {
            finished = true;
            break;
        }
        started = true;

        position_list.push(cur_position);
        let mut current = 0;

        for sam in &sample_names {
            let byteval = sam.as_bytes();
            let vals = vcf_record.genotype(byteval,b"GT").unwrap();
            let value = String::from_utf8(vals[0].clone()).unwrap();
            let first_value = char::to_digit(value.as_bytes()[0] as char,10).unwrap() as u8;
            let second_value = char::to_digit(value.as_bytes()[2] as char,10).unwrap() as u8;
            let first_push = 2*current;
            let second_push = 2*current+1;

            panel_matrix[first_push].push(first_value);
            panel_matrix[second_push].push(second_value);

            current += 1;
        }

    }


    let loaded = VCFData { vcf_data: panel_matrix, positions: position_list,
        sample_names: sample_names, haplotype_names: haplotype_names};

    return Ok(loaded);
}