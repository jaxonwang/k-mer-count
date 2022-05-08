extern crate bio;
extern crate rdxsort;

use std::fs::File;
use bio::io::fasta::{Record, FastaRead, Reader};
use std::str::from_utf8;

fn main() {
    let file = File::open("sample.fasta").expect("Error during opening the file");
    let mut reader = Reader::new(file);
    let mut record = Record::new();

    let l_len = 27;
    let r_len = 27;
    let mut lr_chunk:Vec<String> = Vec::new();
    let mut window_start: usize;
    let mut l_start: usize;
    let mut l_end:   usize;
    let mut r_start: usize;
    let mut r_end:   usize;
    let mut m_len:   usize;

    loop {
        reader.read(&mut record).unwrap();
        if record.is_empty(){
            break;
        }
        for dna_chunk_size in 80..141 {
            window_start = 0;
            loop{
                m_len = dna_chunk_size - l_len - r_len;
                l_start = window_start;
                l_end   = l_start + l_len;
                r_start = l_end + m_len;
                r_end   = r_start + r_len;
                window_start += 1;

                if r_end > record.seq().len(){
                    break;
                }
                /*
                print!("record id: {}\twindow start: {}\tDNA chunk size: {}", record.id(), window_start, dna_chunk_size);
                print!("\tL start: {}\t L end: {}\tR start: {}\tR end: {}", l_start, l_end, r_start, r_end);
                println!("\trecord length: {}\t {}", record.seq().len(), r_end > record.seq().len());
                */
                let l = &record.seq()[l_start..l_end];
                let r = &record.seq()[r_start..r_end];
                //print!("{:?}", from_utf8(l).unwrap());
                //println!("{:?}", from_utf8(r).unwrap());
                let tmp_lr_chunk = from_utf8(l).unwrap().clone().to_owned() + from_utf8(r).unwrap();
                lr_chunk.push(tmp_lr_chunk);
            }
            //println!("{} loop {} DONE", record.id(), dna_chunk_size);
        }
    }
    lr_chunk.sort();
    for each_chunk in lr_chunk.iter() {
        println!("{}", each_chunk);
    }
}