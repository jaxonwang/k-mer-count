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



fn primer3_core_input_sequence(sequences: &Vec<u128>) -> Vec<String>{
    let mut str_vec: Vec<String> = Vec::new();
    let many_n = "N".to_string().repeat(38);

    for each_seq in sequences {
        let l_u8_array = decode_u128_l(&each_seq);
        let m_u8_array = decode_u128_m(&each_seq);
        let r_u8_array = decode_u128_r(&each_seq);
        let l_str: &str = std::str::from_utf8(&l_u8_array).unwrap();
        let m_str: &str = std::str::from_utf8(&m_u8_array).unwrap();
        let r_str: &str = std::str::from_utf8(&r_u8_array).unwrap();
        let sequence_with_internal_n = format!("{}{}{}{}{}", l_str, many_n, m_str, many_n, r_str);
        let primer3_fmt_str = format!("SEQUENCE_ID={:0x}
SEQUENCE_TEMPLATE={}
PRIMER_TASK=pick_pcr_primers_and_hyb_probe
PRIMER_OPT_SIZE=18
PRIMER_MIN_SIZE=15
PRIMER_MAX_SIZE=21
P3_FILE_FLAG=0
PRIMER_EXPLAIN_FLAG=1
=", each_seq, sequence_with_internal_n);
        str_vec.push(primer3_fmt_str);
    }
    return str_vec;
}
//SEQUENCE_TARGET=28,86
//PRIMER_PRODUCT_SIZE_RANGE=100-300


fn execute_primer3(formatted_string: String) -> String{
    let process = match Command::new("primer3_core")
    .stdin(Stdio::piped())
    .stdout(Stdio::piped())
    .spawn() {
        Err(why) => panic!("couldn't spawn primer3: {}", why),
        Ok(process) => process,
    };

    match process.stdin.as_ref().unwrap().write_all(formatted_string.as_bytes()) {
        Err(why) => panic!("couldn't write to primer3_core stdin: {}", why),
        Ok(_) => eprintln!("sent pangram to primer3_core"),
    }

    let output = process.wait_with_output().expect("Failed to wait on child");
    let result = String::from_utf8(output.stdout).unwrap();
    //println!("{}", result);
    return result;
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

    let primer3_fmt_string: Vec<String> = primer3_core_input_sequence(&candidates);

    let mut chunks_of_input: Vec<String> = Vec::new();
    for _i in 0..thread_number{
        chunks_of_input.push(String::new());
    }
    for (index, string) in primer3_fmt_string.iter().enumerate(){
        chunks_of_input[index % thread_number] += string;
        chunks_of_input[index % thread_number] += "\n";
    }


    let arc_chunks_of_input: Arc<Vec<String>> = Arc::new(chunks_of_input);
    let final_result: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
    let mut children = Vec::new();
    for i in 0..thread_number{
        let chunks_of_input  = Arc::clone(&arc_chunks_of_input);
        let arc_final_result = Arc::clone(&final_result);
        //eprintln!("about to spawn thread {}", i);
        children.push(
            thread::spawn(move|| {
                //eprintln!("thread {}: ready to run primer3", i);
                let primer3_results: String = execute_primer3((*chunks_of_input[i]).to_string());
                //eprintln!("thread {}: finish primer3", i);
                arc_final_result.lock().unwrap().push(primer3_results);
        })
        );
    }
    for child in children{
        let _ = child.join();
    }
    for i in final_result.lock().unwrap().iter(){
        println!("{}", i);
    }
/*
    let primer3_fmt_string: Vec<String> = primer3_core_input_sequence(&candidates);
    let stdin_txt: String = primer3_fmt_string.join("\n");
    //println!("{}", stdin_txt);
    let primer3_results: String = execute_primer3(stdin_txt);
    print!("{}", primer3_results);
*/
}