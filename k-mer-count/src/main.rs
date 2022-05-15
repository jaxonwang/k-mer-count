extern crate bio;
extern crate rdxsort;

use std::io;
use std::env;
use std::fs;
use std::fs::File;
use std::str::from_utf8;
use std::io::{BufRead, BufReader};
use std::path::Path;

use bio::io::fastq::Reader as fqReader;
use bio::io::fastq::Record as fqRecord;
use bio::io::fasta::Reader as faReader;
use bio::io::fasta::Record as faRecord;
use crate::bio::io::fastq::FastqRead;
use crate::bio::io::fasta::FastaRead;

use rdxsort::*;
use bit_reverse::ParallelReverse;
use anyhow::Result;
use flate2::read::MultiGzDecoder;

use voracious_radix_sort::{RadixSort};

/*
fn bucket_sort(source: Vec<&str>, place: usize) -> Vec<&str>{
    assert!(source[0].len() >= place);
    let mut result_a: Vec<&str> = Vec::new();
    let mut result_c: Vec<&str> = Vec::new();
    let mut result_g: Vec<&str> = Vec::new();
    let mut result_t: Vec<&str> = Vec::new();
    let index:  usize = source[0].len() - place;
    for each_chunk in source.iter(){
        let current_char = each_chunk.chars().nth(index).unwrap();
        match current_char{
            'A' => {result_a.push(each_chunk);}
            'C' => {result_c.push(each_chunk);}
            'G' => {result_g.push(each_chunk);}
            'T' => {result_t.push(each_chunk);}
            _   => {panic!("Unexpected charactor {} appears in {}", current_char, each_chunk);}
        }
    }
    let mut result: Vec<&str> = Vec::new();
    result.append(&mut result_a);
    result.append(&mut result_c);
    result.append(&mut result_g);
    result.append(&mut result_t);
    return result;
}

fn radix_sort(mut source: Vec<&str>) -> Vec<&str>{
    let length: usize = source[0].len();
    for i in 1..length{
        source = bucket_sort(source, i);
    }
    return source;
}
*/

fn encode_DNA_seq_2_u64(sequence: &[u8]) -> u64{
    let mut result: u64 = 0;
    for each_base in sequence.iter(){
        match each_base{
            b'A' => {result |= 0;}
            b'C' => {result |= 1;}
            b'G' => {result |= 2;}
            b'T' => {result |= 3;}
            _   => {panic!("Unexpected character: {}", each_base);}
        }
        result = result << 2;
    }
    result = result >> 2;
    return result;
}


fn decode_u64_2_DNA_seq(source:u64, index: usize, length: usize) ->u8{
    let mut result: u8 = 0;
    let mut tmp: u64 = source;
    //println!("{:#066b}", tmp);
    tmp = source << (64 - length * 2);
    //println!("{:#066b}", tmp);
    tmp = tmp << index  * 2;
    //println!("{:#066b}", tmp);
    tmp = tmp >> 62;
    //println!("{:#066b}", tmp);

    match tmp{
        0 => {result = b'A';}
        1 => {result = b'C';}
        2 => {result = b'G';}
        3 => {result = b'T';}
        _ => {panic!("Never reached!!!tmp: {}", tmp);}
    }
    return result;

}

pub fn open_with_gz<P: AsRef<Path>>(p: P) -> Result<Box<dyn BufRead>> {
    let r = std::fs::File::open(p.as_ref())?;
    let ext = p.as_ref().extension();

    if ext == Some(std::ffi::OsStr::new("gz")) {
        let gz = MultiGzDecoder::new(r)?;
        let buf_reader = BufReader::new(gz);
        Ok(Box::new(buf_reader))
    } else {
        let buf_reader = BufReader::new(r);
        Ok(Box::new(buf_reader))
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let path = &args[1];
    eprintln!("input {:?}", path);
    //let mut reader = fastq::Reader::new(open_with_gz(path).unwrap());
    let file = File::open(path).expect("Error during opening the file");
    let mut reader = faReader::new(file);
    let mut record = faRecord::new();

    eprintln!("loading {:?} done", path);

    let l_len = 27;
    let r_len = 27;
    let mut lr_chunk:Vec<[u64;2]> = Vec::new();
    let mut window_start: usize;
    let mut l_start: usize;
    let mut l_end:   usize;
    let mut r_start: usize;
    let mut r_end:   usize;
    let mut m_len:   usize;
    let mut loop_cnt:usize = 0;
    loop {
        reader.read(&mut record).unwrap();
        if record.is_empty(){
            break;
        }
        eprintln!("loop: {:09?}, current record id:{:?}\tlength: {:?}", loop_cnt, record.id(), record.seq().len());
        loop_cnt += 1;
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
                let l = &record.seq()[l_start..l_end];
                let r = &record.seq()[r_start..r_end];
                let l_u64: u64 = encode_DNA_seq_2_u64(l);
                let r_u64: u64 = encode_DNA_seq_2_u64(r);
                //println!("{} => {:#066b}", from_utf8(l).unwrap(), l_u64);
/*
                for i in 0..27{
                    let tmp = decode_u64_2_DNA_seq(l_u64, i, l_len);
                    //println!("{}th base: {}", i, tmp);
                }
*/
                //let tmp_lr_chunk = from_utf8(l).unwrap().clone().to_owned() + from_utf8(r).unwrap();
                lr_chunk.push([l_u64, r_u64]);
            }
        }
    }
/*
1/10 or 1/100のリードを調べて当たりをつける
hyper log log counter......?
bloom filterで出現回数の少ないものをカットする
*/

    lr_chunk.rdxsort();
    //lr_chunk.voracious_stable_sort();
    let mut buf: [u8;54] = [64; 54];
    let mut cnt = 0;
    for each_chunk in lr_chunk.iter() {
        for i in 0..l_len{
            buf[cnt] = decode_u64_2_DNA_seq(each_chunk[0], i, l_len);
            cnt+=1;
            //print!("{:?}", char::from(decode_u64_2_DNA_seq(each_chunk[0], i, l_len)));
        }
        for i in 0..r_len{
            buf[cnt] = decode_u64_2_DNA_seq(each_chunk[1], i, l_len);
            cnt+=1;
            //print!("{:?}", char::from(decode_u64_2_DNA_seq(each_chunk[1], i, r_len)));
        }
        cnt = 0;
        println!("{:?}", std::str::from_utf8(&buf).unwrap());
    }
}