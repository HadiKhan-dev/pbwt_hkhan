use serde;
use serde_json;
use crate::pbwt::PbwtInfo;
use std;
use flate2::Compression;
use flate2::write::ZlibEncoder;
use std::io::prelude::*;
use flate2::read::ZlibDecoder;

pub fn write_pbwt(data: &PbwtInfo, filename: &str) -> () {
    let basic_write = serde_json::to_string_pretty(data).unwrap();
    let mut e = ZlibEncoder::new(Vec::new(), Compression::default());
    e.write_all(basic_write.as_bytes());
    let compressed = e.finish().unwrap();
    
    std::fs::write(filename, compressed);
}

pub fn read_pbwt(filename: &str) -> PbwtInfo {
    let mut read_data: Vec<u8> = Vec::new();

    let mut read_file = std::fs::File::open(filename).unwrap();
    read_file.read_to_end(&mut read_data);
    let mut zlib = ZlibDecoder::new(&read_data[..]);
    let mut raw_pbwt_info = String::new();
    zlib.read_to_string(&mut raw_pbwt_info).unwrap();

    let mut pbwt_info = serde_json::from_str::<PbwtInfo>(&raw_pbwt_info).unwrap();

    return pbwt_info;

    
}