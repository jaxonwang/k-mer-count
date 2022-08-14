extern crate bio;
extern crate rdxsort;
extern crate getopts;
use getopts::Options;
use std::{env, process};
use std::fs;
use std::fs::File;
use std::io::{BufWriter, Write};
use voracious_radix_sort::{RadixSort};
use kmer_count::counting_bloomfilter_util::L_LEN;
use kmer_count::counting_bloomfilter_util::R_LEN;
use kmer_count::counting_bloomfilter_util::BLOOMFILTER_TABLE_SIZE;
use kmer_count::counting_bloomfilter_util::{build_counting_bloom_filter, number_of_high_occurence_kmer, pick_up_high_occurence_kmer};
use kmer_count::sequence_encoder_util::{decode_u128_l, decode_u128_r, decode_u128_2_dna_seq};


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

    let threads = if matches.opt_present("t") {
        matches.opt_str("t").unwrap().parse::<usize>().unwrap()
    }else{
        8
    };

    let threshold = if matches.opt_present("a") {
        matches.opt_str("a").unwrap().parse::<u64>().unwrap()
    }else{
        8
    };

    let output_file = if matches.opt_present("o") {
        matches.opt_str("o").unwrap()
    }else{
        format!("{:?}_threshold{}_threads{}.out", input_file, threshold, threads)
    };


    eprintln!("input  file: {:?}",  input_file);
    eprintln!("loading {:?} done", input_file);

    //1段目
    eprintln!("start calling build_counting_bloom_filter");
    let counting_bloom_filter_table: Box<[u64; BLOOMFILTER_TABLE_SIZE]> = build_counting_bloom_filter(&input_file);
    eprintln!("finish calling build_counting_bloom_filter");

    //2段目
    eprintln!("start calling number_of_high_occurence_kmer");
    let (high_occr_bloomfilter_table, occurence) = number_of_high_occurence_kmer(&counting_bloom_filter_table, &input_file, threshold);
    eprintln!("finish calling number_of_high_occurence_kmer");
    //3段目

    eprintln!("Vec size is {}", occurence);
    eprintln!("start calling pick_up_high_occurence_kmer");
    let occr_with_mergin = ((occurence as f64) * 1.2).ceil() as usize;
    let mut high_occurence_kmer: Vec<u128> = pick_up_high_occurence_kmer(&high_occr_bloomfilter_table, &input_file, occr_with_mergin);
    eprintln!("finish calling pick_up_high_occurence_kmer");

    //sortする
    eprintln!("start voracious_mt_sort({})", threads);
    high_occurence_kmer.voracious_mt_sort(threads);
    eprintln!("finish voracious_mt_sort({})", threads);

/*
    let mut previous_l_kmer: [u8; L_LEN] = [b'A'; L_LEN];
    let mut current_l_kmer:  [u8; L_LEN] = [b'A'; L_LEN];
    for each_kmer in high_occurence_kmer{
        current_l_kmer = decode_u128_l(&each_kmer);
        if current_l_kmer != previous_l_kmer{
            println!("{:?}", String::from_utf8(current_l_kmer.to_vec()).unwrap());
        }
        previous_l_kmer = current_l_kmer;
    }
*/
    eprintln!("start  writing to output file: {:?}", &output_file);

    //let mut w = File::create(&output_file).unwrap();
    let mut w = BufWriter::new(fs::File::create(&output_file).unwrap());
    let mut w_kensho = BufWriter::new(fs::File::create("./kensho_out").unwrap());

    let mut previous_kmer: u128 = 0;
    let mut cnt = 0;
    let mut buf_array: [u8; 16] = [0; 16];
    let mut buf_num: u128;

    if matches.opt_present("r") {
        eprintln!("matches.opt_present('r'): {}\tmatches.opt_present('b'): {}", matches.opt_present("r"), matches.opt_present("b"));
        for each_kmer in &high_occurence_kmer{
            if previous_kmer != *each_kmer{
                cnt += 1;
            }
            previous_kmer = *each_kmer;
        }
        writeln!(&mut w, "k-mer count: {}\tthreshold: {}\tinput file {:?}", cnt, threshold, &input_file).unwrap();
    }
    if !matches.opt_present("r") && matches.opt_present("b"){
        eprintln!("matches.opt_present('r'): {}\tmatches.opt_present('b'): {}", matches.opt_present("r"), matches.opt_present("b"));
        for each_kmer in &high_occurence_kmer{
            if previous_kmer != *each_kmer{
                cnt += 1;
                buf_num = *each_kmer;
                for i in 0..16{
                    buf_array[15 - i] = u8::try_from(buf_num & 0xFF).unwrap();
                    buf_num >>= 8;
                }
                w.write(&buf_array).unwrap();
            }
            previous_kmer = *each_kmer;
        }
    }
    if !matches.opt_present("r") && !matches.opt_present("b"){
        eprintln!("matches.opt_present('r'): {}\tmatches.opt_present('b'): {}", matches.opt_present("r"), matches.opt_present("b"));
        for each_kmer in &high_occurence_kmer{
            if previous_kmer != *each_kmer{
                cnt += 1;
                writeln!(&mut w, "{:?}", String::from_utf8(decode_u128_2_dna_seq(&each_kmer, 54)).unwrap()).unwrap();
            }
            previous_kmer = *each_kmer;
        }
    }



    eprintln!("finish writing to output file: {:?}", &output_file);
    eprintln!("total cardinarity of 54-mer: {}", cnt);
    eprintln!("threads: {}\tthreshold: {}\tinput file {:?}", threads, threshold, &input_file);

}