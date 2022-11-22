use serde;
use serde_json;
use crate::pbwt_structs::PbwtInfo;
use std;
use flate2::Compression;
use flate2::write::ZlibEncoder;
use std::io::prelude::*;
use flate2::read::ZlibDecoder;
use bincode;
use serde::{Serialize,Deserialize};
use xz2::write::XzEncoder;
use xz2::read::XzDecoder;
use std::fs::File;

#[derive(Debug,Serialize,Deserialize)]
pub struct Basic {
    pub bin_pbwt_data : Vec<Vec<u8>>,
}

pub fn write_basic(data: &Basic, filename: &str) -> () {

    let basic_write = bincode::serialize(&data).unwrap();

    //let mut e = ZlibEncoder::new(Vec::new(), Compression::new(9));
    //e.write_all(&basic_write[..]);
    //let compressed = e.finish().unwrap();
    
    let mut e = XzEncoder::new(Vec::new(),6);
    e.write_all(&basic_write[..]);
    let compressed = e.finish().unwrap();
    std::fs::write(filename, compressed);

}

pub fn write_pbwt(data: &PbwtInfo, filename: &str) -> () {
    let basic_write = bincode::serialize(&data).unwrap();

    //let mut e = ZlibEncoder::new(Vec::new(), Compression::new(9));
    //e.write_all(&basic_write[..]);
    //let compressed = e.finish().unwrap();

    let mut e = XzEncoder::new(Vec::new(), 6);
    e.write_all(&basic_write[..]);
    let compressed = e.finish().unwrap();
    
    std::fs::write(filename, compressed);
}

pub fn read_pbwt(filename: &str) -> PbwtInfo {

    let mut read_data: Vec<u8> = Vec::new();


    let mut read_file = std::fs::File::open(filename).unwrap();
    read_file.read_to_end(&mut read_data);

    // let mut zlib = ZlibDecoder::new(&read_data[..]);
    // let mut raw_pbwt_info: Vec<u8> = Vec::new();
    // zlib.read_to_end(&mut raw_pbwt_info);

    let mut xz = XzDecoder::new(&read_data[..]);
    let mut raw_pbwt_info: Vec<u8>  = Vec::new();
    xz.read_to_end(&mut raw_pbwt_info);



    //let mut pbwt_info = serde_json::from_str::<PbwtInfo>(&raw_pbwt_info).unwrap();
    let mut pbwt_info = bincode::deserialize::<PbwtInfo>(&raw_pbwt_info).unwrap();

    return pbwt_info;

    
}