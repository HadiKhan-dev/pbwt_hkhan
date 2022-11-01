use std::fs::File;
use csv;
use std::error::Error;
use serde::Deserialize;

#[derive(Deserialize)]
pub struct SiteRow {
    chromosome: u64,
    position: u64,
    reference: String,
    alternate: String,
}

pub fn read_sites(filename: &str) -> Result<Vec<SiteRow>,Box<dyn Error>> {

    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(false).from_path(filename)?;

    let mut full_rows : Vec<SiteRow> = Vec::new();

    //let headers = reader.headers()?;
    let header = csv::StringRecord::from(vec!["chromosome","position","reference","alternate"]);

    for result in rdr.records() {
        let record = result?;
        let dsr: SiteRow = record.deserialize(Some(&header))?;
        full_rows.push(dsr);
        println!("{:?}",&full_rows[full_rows.len()-1].alternate);
    }
    return Ok(full_rows);
}