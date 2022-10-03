extern crate kmer_count;
extern crate getopts;
use std::{env, process};
use std::fs::File;
use std::io::{Read, BufReader};
use std::io::prelude::*;
use std::process::{Command, Stdio};
use std::thread;
use std::sync::Arc;
use std::sync::Mutex;
use getopts::Options;
use kmer_count::sequence_encoder_util::{decode_u128_l, decode_u128_m, decode_u128_r};



fn sequence_dump(sequences: &Vec<u128>) -> Vec<String>{
    let mut str_vec: Vec<String> = Vec::new();
    for each_seq in sequences {
        let l_u8_array  = decode_u128_l(&each_seq);
        let m_u8_array  = decode_u128_m(&each_seq);
        let r_u8_array  = decode_u128_r(&each_seq);
        let l_str: &str = std::str::from_utf8(&l_u8_array).unwrap();
        let m_str: &str = std::str::from_utf8(&m_u8_array).unwrap();
        let r_str: &str = std::str::from_utf8(&r_u8_array).unwrap();
        let sequence    = format!("{}{}{}", l_str, m_str, r_str);
        str_vec.push(sequence);
    }
    return str_vec;
}




fn print_usage(program: &str, opts: &Options) {
    let brief = format!("Usage: {} FILE ", program);
    print!("{}", opts.usage(&brief));
    process::exit(0);
}

fn main(){
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();
    opts.optflag("h", "help", "print this help menu");
    opts.optopt("t", "thread", "number of thread to use for radix sort. default value is 8.", "THREAD");


    let matches = match opts.parse(&args[1..]) {
        Ok(m) => { m }
        Err(f) => { panic!("{}", f.to_string()) }
    };
    if matches.opt_present("h") {
        print_usage(&program, &opts);
        return;
    }

    let thread_number: usize = if matches.opt_present("t") {
        matches.opt_str("t").unwrap().parse::<usize>().unwrap()
    }else{
        4
    };

    let input_file = if !matches.free.is_empty() {
        matches.free[0].clone()
    } else {
        print_usage(&program, &opts);
        return;
    };
    let f: File = File::open(&input_file).unwrap();
    let mut reader = BufReader::new(f);
    let mut buf: [u8; 16] = [0; 16];
    let mut tmp_seq_as_u128: u128 = 0;
    let mut candidates: Vec<u128> = Vec::new();

    loop {
        match reader.read(&mut buf).unwrap() {
            0 => break,
            n => {
                let buf = &buf[..n];
                for i in 0..16 {
                    tmp_seq_as_u128 <<= 8;
                    tmp_seq_as_u128 += u128::from(buf[i]);
                }
                candidates.push(tmp_seq_as_u128);
            }
        }
    }

    let sequence_dump_string: Vec<String> = sequence_dump(&candidates);
    for candidate in sequence_dump_string{
        println!("{}", candidate);
    }
}