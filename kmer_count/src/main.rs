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

//use voracious_radix_sort::{RadixSort};

//use rand::prelude::*;

//use kmer_count::encoder_util::decode_u128_2_dna_seq;
//use kmer_count::encoder_util::encode_dna_seq_2_u128;
//use kmer_count::encoder_util::decode_u128_2_occurence;
//use kmer_count::encoder_util::L_LEN;
//use kmer_count::encoder_util::R_LEN;
use kmer_count::counting_bloomfilter_util::BLOOMFILTER_TABLE_SIZE;
use kmer_count::counting_bloomfilter_util::{build_counting_bloom_filter};


//const TOW_SQ20: u128 = 2_u128.pow(20);




fn main() {
    let args: Vec<String> = env::args().collect();
    let path = &args[1];
    eprintln!("input {:?}", path);
//    let file = File::open(path).expect("Error during opening the file");
//    let mut reader = faReader::new(file);
//    let mut record = faRecord::new();

    eprintln!("loading {:?} done", path);

    //1段目
    eprintln!("start calling build_counting_bloom_filter");
    let counting_bloom_filter_table: Box<[u64; BLOOMFILTER_TABLE_SIZE]> = build_counting_bloom_filter(path);
    eprintln!("finish calling build_counting_bloom_filter");

    //2段目
    //let high_occr_cnt: u64 = number_of_high_occurence_kmer(&counting_bloom_filter_table, path);
    //3段目
    //let mut high_occurence_kmer: Vec<u128> = pick_up_high_occurence_kmer(&counting_bloom_filter_table, path, high_occr_cnt);

    //どんなふうに出力しようか？

/*
    high_occurence_kmer.voracious_mt_sort(8);
    for each_kmer in high_occurence_kmer{
        println!("{:?}", each_kmer);
    }
*/
}