use serde::Deserialize;

#[derive(Debug,Deserialize,Clone)]
pub struct SiteRow {
    pub chromosome: u64,
    pub position: u64,
    pub reference: String,
    pub alternate: String,
}


#[derive(Debug,Clone)]
pub struct VCFData {
    pub vcf_data: Vec<Vec<u8>>,
    pub chromosomes: Vec<u64>,
    pub positions: Vec<u64>,
    pub sample_names: Vec<String>,
    pub haplotype_names: Vec<String>,
}