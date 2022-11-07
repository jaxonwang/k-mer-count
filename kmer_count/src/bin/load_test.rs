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








fn print_usage(program: &str, opts: &Options) {
    let brief = format!("Usage: {} FILE [options]", program);
    print!("{}", opts.usage(&brief));
    process::exit(0);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();
    opts.optopt("o", "output", "set output file name", "NAME");
    opts.optopt("t", "thread", "number of threads to use for radix sort. default value is 8.", "THREAD");
    opts.optopt("a", "threshold", "threshold for hyper log counter. default value is 8.", "THRESHOLD");
    opts.optflag("b", "binary", "outputs binary file");
    opts.optflag("r", "only-num", "outputs only total number of k-mer");
    opts.optflag("h", "help", "print this help menu");

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

    let threads: usize = if matches.opt_present("t") {
        matches.opt_str("t").unwrap().parse::<usize>().unwrap()
    }else{
        8
    };

    let threshold:u32 = if matches.opt_present("a") {
        matches.opt_str("a").unwrap().parse::<u32>().unwrap()
    }else{
        8
    };

    let output_file = if matches.opt_present("o") {
        matches.opt_str("o").unwrap()
    }else{
        format!("{:?}_threshold{}_threads{}.out", input_file, threshold, threads)
    };


    let file = File::open(input_file).expect("Error during opening the file");
    let mut reader = faReader::new(file);
    let mut record = faRecord::new();
    let mut sequences: Vec<DnaSequence> = Vec::new();
    'each_read: loop {
        reader.read(&mut record).unwrap();
        if record.is_empty(){
            break 'each_read;
        }
        let sequence_as_vec: Vec<u8> = record.seq().to_vec();
        let current_sequence = DnaSequence::new(&sequence_as_vec);
        sequences.push(current_sequence);
    }
}