use serde::{Serialize,Deserialize};

#[derive(Debug,Serialize,Deserialize,Clone)]
pub struct PbwtInfo {
    pub num_samples: u32,
    pub num_sites: u32,
    pub bin_pbwt: Vec<Vec<u8>>,
    //pub pbwt_data : Vec<Vec<u32>>,
    //pub divergence_array : Vec<Vec<u32>>,
    pub count : Vec<u32>,
    pub occ_list: Vec<Vec<Vec<u32>>>,
    pub fm_gap: u32,
}

#[derive(Debug,Serialize,Deserialize,Clone)]
pub struct SpacedPbwt {
    pub num_samples: u32,
    pub num_pbwt_sites: u32,
    pub num_inserted_sites: u32,
    pub num_total_sites: u32,

    pub pbwt_positions: Vec<u32>,
    pub inserted_positions: Vec<u32>,
    pub all_positions: Vec<u32>,
    pub pbwt_col_flags: Vec<u8>,

    pub bin_pbwt: Vec<Vec<u8>>,
    pub count: Vec<u32>,
    pub occ_list: Vec<Vec<Vec<u32>>>,
    pub fm_gap: u32,
}

#[derive(Debug,Serialize,Deserialize,Clone)]
pub struct DualPbwt {
    pub forward_pbwt: SpacedPbwt,
    pub reverse_pbwt: SpacedPbwt,
}