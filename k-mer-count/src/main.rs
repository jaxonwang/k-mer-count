extern crate bio;

use std::env;
use std::fs::File;
use std::io::prelude::*;
use bio::io::fasta;
use bio::io::fasta::{Record, FastaRead, Reader};


fn main() {
    let file = File::open("sample.fasta").expect("Error during opening the file");
    let mut reader = Box::new(fasta::Reader::new(file));
    let mut record = Record::new();
    loop {
        reader.read(&mut record).unwrap();
        if record.is_empty(){
            break;
        }
        println!("{}\t{:?}\t{:?}", record.id(), record.desc().unwrap(), std::str::from_utf8(record.seq()).unwrap());
    }
}