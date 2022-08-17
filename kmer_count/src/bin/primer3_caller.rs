extern crate kmer_count;
extern crate getopts;
use std::{env, process};
use std::fs::File;
use std::io::{Read, BufReader};
use std::io::prelude::*;
use std::error::Error;
use std::process::{Command, Stdio};
use std::borrow::Cow;
use std::collections::HashMap;
use std::thread;
use std::sync::Arc;
use std::sync::Mutex;
use getopts::Options;

use kmer_count::sequence_encoder_util::{decode_u128_l, decode_u128_r, decode_u128_2_dna_seq};



/*
1. バイナリファイルを読み込んで、128bitごとに切る。
2. primer3_coreに投げられるように、文字列を生成する。
3. primer3_coreに投げる。

投げる際、上手に並列にやりたい。


*/

/*
SEQUENCE_ID=example
SEQUENCE_TEMPLATE=GTAGTCAGTAGACNATGACNACTGACGATGCAGACNACACACACACACACAGCACACAGGTATTAGTGGGCCATTCGATCCCGACCCAAATCGATAGCTACGATGACG
SEQUENCE_TARGET=37,21
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=18
PRIMER_MIN_SIZE=15
PRIMER_MAX_SIZE=21
PRIMER_MAX_NS_ACCEPTED=1
PRIMER_PRODUCT_SIZE_RANGE=100-150
P3_FILE_FLAG=1
SEQUENCE_INTERNAL_EXCLUDED_REGION=37,21
PRIMER_EXPLAIN_FLAG=1
=

*/

/*OUTPUT example:
SEQUENCE_ID=example1
SEQUENCE_TEMPLATE=tattggtgaagcctcaggtagtgcagaatatgaaacttcaggatccagtgggcatgctactggtagtgctgccggccttacaggcattatggtggcaaagtcgacagagttta
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=22
PRIMER_PRODUCT_SIZE_RANGE=75-150
PRIMER_EXPLAIN_FLAG=1
P3_FILE_FLAG=1
PRIMER_EXPLAIN_FLAG=1
PRIMER_LEFT_EXPLAIN=considered 195, low tm 122, high tm 6, ok 67
PRIMER_RIGHT_EXPLAIN=considered 195, low tm 47, high tm 51, high hairpin stability 26, ok 71
PRIMER_PAIR_EXPLAIN=considered 7, unacceptable product size 2, ok 5
PRIMER_LEFT_NUM_RETURNED=5
PRIMER_RIGHT_NUM_RETURNED=5
PRIMER_INTERNAL_NUM_RETURNED=0
PRIMER_PAIR_NUM_RETURNED=5
PRIMER_PAIR_0_PENALTY=0.571539
PRIMER_LEFT_0_PENALTY=0.321671
PRIMER_RIGHT_0_PENALTY=0.249868
PRIMER_LEFT_0_SEQUENCE=gtgaagcctcaggtagtgca
PRIMER_RIGHT_0_SEQUENCE=ctctgtcgactttgccacca
PRIMER_LEFT_0=5,20
PRIMER_RIGHT_0=108,20
PRIMER_LEFT_0_TM=59.678
PRIMER_RIGHT_0_TM=60.250
PRIMER_LEFT_0_GC_PERCENT=55.000
PRIMER_RIGHT_0_GC_PERCENT=55.000
PRIMER_LEFT_0_SELF_ANY_TH=0.00
PRIMER_RIGHT_0_SELF_ANY_TH=17.92
PRIMER_LEFT_0_SELF_END_TH=0.00
PRIMER_RIGHT_0_SELF_END_TH=0.00
PRIMER_LEFT_0_HAIRPIN_TH=30.14
PRIMER_RIGHT_0_HAIRPIN_TH=0.00
PRIMER_LEFT_0_END_STABILITY=4.5700
PRIMER_RIGHT_0_END_STABILITY=4.1700
PRIMER_PAIR_0_COMPL_ANY_TH=0.00
PRIMER_PAIR_0_COMPL_END_TH=2.16
PRIMER_PAIR_0_PRODUCT_SIZE=104
PRIMER_PAIR_0_PRODUCT_TM=84.0
=

*/


fn primer3_core_input_sequence(sequences: &Vec<u128>) -> Vec<String>{
    let mut str_vec: Vec<String> = Vec::new();
    let many_n = "N".to_string().repeat(86);

    for each_seq in sequences {
        let l_u8_array = decode_u128_l(&each_seq);
        let r_u8_array = decode_u128_r(&each_seq);
        let l_str: &str = std::str::from_utf8(&l_u8_array).unwrap();
        let r_str: &str = std::str::from_utf8(&r_u8_array).unwrap();
        let sequence_with_internal_n = format!("{}{}{}", l_str, many_n, r_str);
        let primer3_fmt_str = format!("SEQUENCE_ID={:0x}
SEQUENCE_TEMPLATE={}
SEQUENCE_TARGET=28,86
PRIMER_TASK=pick_pcr_primers
PRIMER_OPT_SIZE=18
PRIMER_MIN_SIZE=15
PRIMER_MAX_SIZE=21
PRIMER_PRODUCT_SIZE_RANGE=201-300 101-200
P3_FILE_FLAG=0
PRIMER_EXPLAIN_FLAG=1
=", each_seq, sequence_with_internal_n);
        str_vec.push(primer3_fmt_str);
    }
    return str_vec;
}


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
/*
struct Candidate{
    SEQUENCE_ID                  : String,
    SEQUENCE_TEMPLATE            : String,
    PRIMER_TASK                  : String,
    PRIMER_PICK_LEFT_PRIMER      : u8,
    PRIMER_PICK_INTERNAL_OLIGO   : u8,
    PRIMER_PICK_RIGHT_PRIMER     : u8,
    PRIMER_OPT_SIZE              : u8,
    PRIMER_MIN_SIZE              : u8,
    PRIMER_MAX_SIZE              : u8,
    PRIMER_PRODUCT_SIZE_RANGE    : Vec<[u8; 2]>,
    P3_FILE_FLAG                 : bool,
    PRIMER_EXPLAIN_FLAG          : bool,
    PRIMER_LEFT_EXPLAIN          : String,
    PRIMER_RIGHT_EXPLAIN         : String,
    PRIMER_PAIR_EXPLAIN          : String,
    PRIMER_LEFT_NUM_RETURNED     : u8,
    PRIMER_RIGHT_NUM_RETURNED    : u8,
    PRIMER_INTERNAL_NUM_RETURNED : u8,
    PRIMER_PAIR_NUM_RETURNED     : u8
}

impl Candidate{
    fn new(source: String) -> Candidate{

        /*
        タグごとに処理が変わる
        _4_みたいなのは同じグループにしたい。
        Pythonで書きたい......。
        */


        //println!("Candidate::new() recieved \n{}\n\n\n1", source);
        let mut dict = HashMap::new();
        for each_line in source.split("\n"){
            if each_line == ""{continue}
            let key = each_line.split("=").collect::<Vec<&str>>()[0];
            let val = each_line.split("=").collect::<Vec<&str>>()[1];
            dict.insert(key, val);
        }
        return Candidate{
SEQUENCE_ID                  : "example1".to_string(),
SEQUENCE_TEMPLATE            : "agtcgacagagttta".to_string(),
PRIMER_TASK                  : "generic".to_string(),
PRIMER_PICK_LEFT_PRIMER      : 1,
PRIMER_PICK_INTERNAL_OLIGO   : 0,
PRIMER_PICK_RIGHT_PRIMER     : 1,
PRIMER_OPT_SIZE              : 20,
PRIMER_MIN_SIZE              : 18,
PRIMER_MAX_SIZE              : 22,
PRIMER_PRODUCT_SIZE_RANGE    : vec![[75, 150]],
P3_FILE_FLAG                 : false,
PRIMER_EXPLAIN_FLAG          : false,
PRIMER_LEFT_EXPLAIN          : "considered 195, low tm 122, high tm 6, ok 67".to_string(),
PRIMER_RIGHT_EXPLAIN         : "considered 195, low tm 47, high tm 51, high hairpin stability 26, ok 71".to_string(),
PRIMER_PAIR_EXPLAIN          : "considered 7, unacceptable product size 2, ok 5".to_string(),
PRIMER_LEFT_NUM_RETURNED     : 5,
PRIMER_RIGHT_NUM_RETURNED    : 5,
PRIMER_INTERNAL_NUM_RETURNED : 0,
PRIMER_PAIR_NUM_RETURNED     : 5,
            }
        /*
            Candidate１個分のStringを入力とする
            末尾に"="は含まれない
        */
    }
}



struct PrimerPair{
    primer_pair_penalty        : String,
    primer_left_penalty        : String,
    primer_right_penalty       : String,
    primer_left_sequence       : String,
    primer_right_sequence      : String,
    primer_left                : String,
    primer_right               : String,
    primer_left_tm             : String,
    primer_right_tm            : String,
    primer_left_gc_percent     : String,
    primer_right_gc_percent    : String,
    primer_left_self_any_th    : String,
    primer_right_self_any_th   : String,
    primer_left_self_end_th    : String,
    primer_right_self_end_th   : String,
    primer_left_hairpin_th     : String,
    primer_right_hairpin_th    : String,
    primer_left_end_stability  : String,
    primer_right_end_stability : String,
    primer_pair_compl_any_th   : String,
    primer_pair_compl_end_th   : String,
    primer_pair_product_size   : String,
    primer_pair_product_tm     : String
}



fn primer3_result_parser(source: String) -> Vec<Candidate>{
    let mut ret_val: Vec<Candidate> = Vec::new();
    let mut buf_str: String = String::new();
    for each_line in source.split("\n") {
        if each_line != "="{
            buf_str += each_line;
            buf_str += "\n";
        }else{
            //println!("{}\n\n", buf_str);
            ret_val.push(Candidate::new(buf_str));
            buf_str = String::new();
        }
    }
    return ret_val
}
*/




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
    for i in 0..thread_number{
        chunks_of_input.push(String::new());
    }
    for (index, string) in primer3_fmt_string.iter().enumerate(){
        chunks_of_input[index % thread_number] += string;
        chunks_of_input[index % thread_number] += "\n";
    }

    let arc_chunks_of_input: Arc<Vec<String>> = Arc::new(chunks_of_input);
    let mut final_result: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
    let mut children = Vec::new();
    for i in 0..thread_number{
        let chunks_of_input  = Arc::clone(&arc_chunks_of_input);
        let arc_final_result = Arc::clone(&final_result);
        eprintln!("about to spawn thread {}", i);
        children.push(
            thread::spawn(move|| {
                eprintln!("thread {}: ready to run primer3", i);
                let primer3_results: String = execute_primer3((*chunks_of_input[i]).to_string());
                eprintln!("thread {}: finish primer3", i);
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