extern crate bio;
extern crate rdxsort;
extern crate getopts;
use getopts::Options;
use std::{env, process};
use std::fs;
use std::io::{BufWriter, Write};
use kmer_count::counting_bloomfilter_util::BLOOMFILTER_TABLE_SIZE;
use kmer_count::counting_bloomfilter_util::{L_LEN, M_LEN, R_LEN};
use kmer_count::sequence_encoder_util::{decode_u128_2_dna_seq, DnaSequence};

use std::fs::File;
use bio::io::fasta::Reader as faReader;
use bio::io::fasta::Record as faRecord;
use bio::io::fasta::FastaRead;
use std::collections::HashSet;

fn print_usage(program: &str, opts: &Options) {
    let brief = format!("Usage: {} FILE", program);
    print!("{}", opts.usage(&brief));
    process::exit(0);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => { m }
        Err(f) => { panic!("{}", f.to_string()) }
    };
    if matches.opt_present("h") {
        print_usage(&program, &opts);
        return;
    }

    let input_file = if !matches.free.is_empty() {
        matches.free[0].clone()
    } else {
        print_usage(&program, &opts);
        return;
    };

    let file = File::open(input_file).expect("Error during opening the file");
    let mut reader = faReader::new(file);
    let mut record = faRecord::new();
    let mut hash_count: usize = 0;
    let mut hash_from_u128: HashSet<u128> = HashSet::new();
    'each_read: loop {
        reader.read(&mut record).unwrap();
        if record.is_empty(){
            break 'each_read;
        }
        let sequence_as_vec: Vec<u8> = record.seq().to_vec();
        let current_sequence = DnaSequence::new(&sequence_as_vec);
        for i in 0..current_sequence.len(){
            let subseq: u128 = current_sequence.subsequence_as_u128(vec!([i, i + 64]));
            hash_from_u128.insert(subseq);
            hash_count += 1;
            eprintln!("hash_count: {}", hash_count);
        }
    }
}