use serde::Deserialize;

#[derive(Debug,Deserialize,Clone)]
pub struct SiteRow {
    pub chromosome: String,
    pub position: u32,
    pub reference: String,
    pub alternate: String,
}


#[derive(Debug,Clone)]
pub struct VCFData {
    pub vcf_data: Vec<Vec<u8>>,
    pub chromosomes: Vec<String>,
    pub positions: Vec<u32>,
    pub sample_names: Vec<String>,
    pub haplotype_names: Vec<String>,
    pub alt_allele_freq: Vec<f64>
}