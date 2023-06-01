pub struct PositionData {
    pub position: isize,
    pub side_data: Vec<Vec<u8>>,
    pub divergence_data: Vec<Vec<u32>>,
}

pub struct FullInsertData {
    pub insert_positions: Vec<isize>,
    pub side_values: Vec<Vec<Vec<u8>>>,
    pub divergence_values: Vec<Vec<Vec<u32>>>
}