extern crate bio;
extern crate rdxsort;

//use std::io;
use std::env;
//use std::fs;
use std::fs::File;
//use std::str::from_utf8;
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

const L_LEN: usize = 27;
const R_LEN: usize = 27;
const TOW_SQ20: u128 = 2_u128.pow(20);

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

fn count_occurence(input_sequence: Vec<u128>) -> Vec<u128>{
    let mut current_sequence: u128 = input_sequence[0];
    let mut ret_vec: Vec<u128> = Vec::new();
    let mut index: usize = 0;
    let mut counter: u128 = 1;
    let mut buf: u128;
    let goal: usize = input_sequence.len();
    loop{
        index += 1;
        if index >= goal{
            break ret_vec;
        }
        //in case that index reaches the border between different sequence.
        if input_sequence[index] != current_sequence{
        let dna_string = String::from_utf8(decode_u128_2_dna_seq(&current_sequence).to_vec()).unwrap();
        if counter > TOW_SQ20{
            eprintln!("count_occurence unexpected situation: {} appears more than {}", dna_string, TOW_SQ20);
            eprintln!("{}\t{:?}", counter, dna_string);
            break ret_vec;
        }

/*
        println!("counter: {}", counter);
        println!("current_sequence: {}", dna_string);
        println!("counter:          {:#0130b}", counter);//バイナリ列で表示する
*/
        counter = counter << (L_LEN + R_LEN) * 2;
/*
        println!("counter(shifted): {:#0130b}", counter);//バイナリ列で表示する
        println!("current_sequence: {:#0130b}", current_sequence);//バイナリ列で表示する
*/
        buf = counter + current_sequence;
/*
        println!("sum of them:      {:#0130b}",buf);//バイナリ列で表示する
        println!("{}", decode_u128_2_occurence(&buf));
        println!("{}", String::from_utf8(decode_u128_2_dna_seq(&buf).to_vec()).unwrap());
        println!();
*/
        ret_vec.push(buf);
        counter = 1;
        current_sequence = input_sequence[index];
        }else{
            counter += 1;
        }
    }
}

fn encode_dna_seq_2_u64(sequence: &[u8]) -> u64{
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

fn encode_dna_seq_2_u128(sequence: &[u8]) -> u128{
    let mut result: u128 = 0;
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



fn decode_u64_2_dna_seq(source:u64, index: usize, length: usize) ->u8{
    let mut result: u8 = 0;
    let mut tmp: u64 = source << (64 - length * 2);
    tmp = tmp << index * 2;
    tmp = tmp >> 62;
    match tmp{
        0 => {result = b'A';}
        1 => {result = b'C';}
        2 => {result = b'G';}
        3 => {result = b'T';}
        _ => {panic!("Never reached!!!tmp: {}", tmp);}
    }
    return result;
}

fn decode_u128_2_dna_seq(source:&u128) -> [u8; L_LEN + R_LEN]{
    let mut result: [u8; L_LEN + R_LEN] = [b'X'; L_LEN + R_LEN];
    let mut tmp: u128 = source.clone() << 20;
    let mut base: u128;
    for i in 0..L_LEN + R_LEN{
        base = (tmp & 0xC000_0000_0000_0000_0000_0000_0000_0000) >> 126;
        //println!("{}, {:#130b}", base, tmp);
        tmp = tmp << 2;
        match base{
            0 => {result[i] = b'A';}
            1 => {result[i] = b'C';}
            2 => {result[i] = b'G';}
            3 => {result[i] = b'T';}
            _ => {panic!("Never reached!!!tmp: {}", tmp);}
        }
    }
    return result;
}
fn decode_u128_2_occurence(source: &u128) -> u32{
    return TryFrom::try_from(source >> (L_LEN + R_LEN) * 2).unwrap();
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

    //let mut lr_chunk:Vec<[u64;2]> = Vec::new();
    let mut lr_chunk:Vec<u128> = Vec::new();
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
                m_len = dna_chunk_size - L_LEN - R_LEN;
                l_start = window_start;
                l_end   = l_start + L_LEN;
                r_start = l_end + m_len;
                r_end   = r_start + R_LEN;
                window_start += 1;

                if r_end > record.seq().len(){
                    break;
                }
                let l = &record.seq()[l_start..l_end];
                let r = &record.seq()[r_start..r_end];
                let mut lr_string: [u8;L_LEN + R_LEN] = [64; L_LEN + R_LEN];
                for i in 0..L_LEN{
                    lr_string[i] = l[i];
                }
                for i in 0..R_LEN{
                    lr_string[i + L_LEN] = r[i];
                }
                let lr_u128 = encode_dna_seq_2_u128(&lr_string);
                //let l_u64: u64 = encode_dna_seq_2_u64(l);
                //let r_u64: u64 = encode_dna_seq_2_u64(r);
                //lr_chunk.push([l_u64, r_u64]);
                lr_chunk.push(lr_u128);
            }
        }
    }
/*
1/10 or 1/100のリードを調べて当たりをつける
hyper log log counter......?
bloom filterで出現回数の少ないものをカットする
*/

    //lr_chunk.rdxsort();
    lr_chunk.voracious_mt_sort(8);
    let sorted_lr_chunk: Vec<u128> = count_occurence(lr_chunk);
    for each_chunk in sorted_lr_chunk.iter() {
        let dna_string = String::from_utf8(decode_u128_2_dna_seq(each_chunk).to_vec()).unwrap();
        let occurrence = decode_u128_2_occurence(&each_chunk);
        println!("{}\t{}", occurrence, dna_string);

    }
/*
    for each_chunk in lr_chunk.iter() {
        let dna_string = String::from_utf8(decode_u128_2_dna_seq(each_chunk).to_vec()).unwrap();
        println!("{:?}", dna_string);
    }
*/
/*    let mut cnt = 0;
    for each_chunk in lr_chunk.iter() {
        for i in 0..L_LEN{
            buf[cnt] = decode_u64_2_dna_seq(each_chunk[0], i, L_LEN);
            cnt+=1;
            //print!("{:?}", char::from(decode_u64_2_dna_seq(each_chunk[0], i, L_LEN)));
        }
        for i in 0..R_LEN{
            buf[cnt] = decode_u64_2_dna_seq(each_chunk[1], i, L_LEN);
            cnt+=1;
            //print!("{:?}", char::from(decode_u64_2_dna_seq(each_chunk[1], i, R_LEN)));
        }
        cnt = 0;
        println!("{:?}", std::str::from_utf8(&buf).unwrap());
    }*/

}