use serde::Deserialize;

#[derive(Deserialize)]
pub struct SiteRow {
    pub chromosome: u64,
    pub position: u64,
    pub reference: String,
    pub alternate: String,
}


#[derive(Debug)]
pub struct VCFData {
    pub vcf_data: Vec<Vec<u8>>,
    pub positions: Vec<u64>,
    pub sample_names: Vec<String>,
    pub haplotype_names: Vec<String>,
}