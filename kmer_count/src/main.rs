extern crate bio;
extern crate rdxsort;

use std::env;
//use std::fs;
//use std::fs::File;
//use std::str::from_utf8;
//use std::io::{BufRead, BufReader};
//use std::path::Path;
//use std::collections::hash_map::DefaultHasher;
//use sha2::{Sha256, Sha512, Digest};


//use bio::io::fastq::Reader as fqReader;
//use bio::io::fastq::Record as fqRecord;
//use bio::io::fasta::Reader as faReader;
//use bio::io::fasta::Record as faRecord;
//use crate::bio::io::fastq::FastqRead;
//use crate::bio::io::fasta::FastaRead;

//use rdxsort::*;
//use bit_reverse::ParallelReverse;
//use anyhow::Result;
//use flate2::read::MultiGzDecoder;

use voracious_radix_sort::{RadixSort};

//use rand::prelude::*;

//use kmer_count::encoder_util::decode_u128_2_dna_seq;
//use kmer_count::encoder_util::encode_dna_seq_2_u128;
//use kmer_count::encoder_util::decode_u128_2_occurence;
use kmer_count::counting_bloomfilter_util::L_LEN;
use kmer_count::counting_bloomfilter_util::R_LEN;
use kmer_count::counting_bloomfilter_util::BLOOMFILTER_TABLE_SIZE;
use kmer_count::counting_bloomfilter_util::{build_counting_bloom_filter, number_of_high_occurence_kmer, pick_up_high_occurence_kmer};


//const TOW_SQ20: u128 = 2_u128.pow(20);


fn decode_u128_2_dna_seq(source:&u128, char_size: usize) -> Vec<u8>{
    let mut result: Vec<u8> = Vec::new();
    let mut base;
    for i in 0..char_size{
        base = source >> 2 * (char_size - 1 - i) & 3;
        match base{
            0 => {result.push(b'A');}
            1 => {result.push(b'C');}
            2 => {result.push(b'G');}
            3 => {result.push(b'T');}
            _ => {panic!("Never reached!!!base: {}", base);}
        }
    }
    return result;
}

fn decode_u128_l(source: &u128) -> [u8; L_LEN]{
    let mut result: [u8; L_LEN] = [b'X'; L_LEN];
    let mut base;
    for i in 0..L_LEN{
        base = source >> ((R_LEN + i) * 2) & 3;
        match base{
            0 => {result[i] = b'A';}
            1 => {result[i] = b'C';}
            2 => {result[i] = b'G';}
            3 => {result[i] = b'T';}
            _ => {panic!("Never reached!!!base: {}", base);}
        }
    }
    return result;
}

fn decode_u128_r(source: &u128) -> [u8; R_LEN]{
    let mut result: [u8; R_LEN] = [b'X'; R_LEN];
    let mut base;
    for i in 0..R_LEN{
        base = source >> (i * 2) & 3;
        match base{
            0 => {result[i] = b'A';}
            1 => {result[i] = b'C';}
            2 => {result[i] = b'G';}
            3 => {result[i] = b'T';}
            _ => {panic!("Never reached!!!base: {}", base);}
        }
    }
    return result;
}



fn main() {
    let args: Vec<String> = env::args().collect();
    let  input_path = &args[1];
    //let output_path = &args[2];

    eprintln!("input  file: {:?}",  input_path);
    //eprintln!("output file: {:?}", output_path);

//    let file = File::open(input_path).expect("Error during opening the file");
//    let mut reader = faReader::new(file);
//    let mut record = faRecord::new();

    eprintln!("loading {:?} done", input_path);

    //1段目
    eprintln!("start calling build_counting_bloom_filter");
    let counting_bloom_filter_table: Box<[u64; BLOOMFILTER_TABLE_SIZE]> = build_counting_bloom_filter(input_path);
    eprintln!("finish calling build_counting_bloom_filter");

    //2段目
    eprintln!("start calling number_of_high_occurence_kmer");
    let (high_occr_bloomfilter_table, occurence) = number_of_high_occurence_kmer(&counting_bloom_filter_table, input_path, 8);
    eprintln!("finish calling number_of_high_occurence_kmer");
    //3段目

    eprintln!("Vec size is {}", occurence);
    eprintln!("start calling pick_up_high_occurence_kmer");
    let occr_with_mergin = ((occurence as f64) * 1.2).ceil() as usize;
    let mut high_occurence_kmer: Vec<u128> = pick_up_high_occurence_kmer(&high_occr_bloomfilter_table, input_path, occr_with_mergin);
    eprintln!("finish calling pick_up_high_occurence_kmer");


    high_occurence_kmer.voracious_mt_sort(8);
    let mut previous_kmer: u128 = 0;
    let mut previous_l_kmer: [u8; L_LEN] = [b'X'; L_LEN];
    let mut current_l_kmer:  [u8; L_LEN] = [b'X'; L_LEN];

    for each_kmer in high_occurence_kmer{
        current_l_kmer = decode_u128_l(&each_kmer);
        if current_l_kmer != previous_l_kmer{
            println!("{:?}", current_l_kmer);
        }
        previous_l_kmer = current_l_kmer;

    }


    let mut cnt = 0;
/*
    for each_kmer in high_occurence_kmer{
        if previous_kmer != each_kmer{
            println!("{:?}", String::from_utf8(decode_u128_2_dna_seq(&each_kmer, 54)).unwrap());
            cnt += 1;
            if cnt >= 1000{
                break;
            }
        }
        previous_kmer = each_kmer;
    }
*/


}