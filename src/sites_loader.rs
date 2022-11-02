use std::fs::File;
use csv;
use std::error::Error;
use serde::Deserialize;

use crate::vcf_structs::SiteRow;


pub fn read_sites(filename: &str) -> Result<Vec<SiteRow>,Box<dyn Error>> {

    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(filename)?;

    let mut full_rows : Vec<SiteRow> = Vec::new();

    let header = csv::StringRecord::from(vec!["chromosome","position","reference","alternate"]);

    for result in rdr.records() {
        let record = result?;
        let dsr: SiteRow = record.deserialize(Some(&header))?;
        full_rows.push(dsr);
    }
    return Ok(full_rows);
}