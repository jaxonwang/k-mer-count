extern crate bio;
extern crate rdxsort;
extern crate getopts;
use getopts::Options;
use std::{env, process};
use std::fs;
use std::io::{BufWriter, Write};
use std::collections::HashSet;
use std::thread;
use std::sync::Arc;
use voracious_radix_sort::{RadixSort};
use kmer_count::counting_bloomfilter_util::BLOOMFILTER_TABLE_SIZE;
use kmer_count::counting_bloomfilter_util::{L_LEN, M_LEN, R_LEN};
use kmer_count::counting_bloomfilter_util::{build_counting_bloom_filter, number_of_high_occurence_kmer};
use kmer_count::sequence_encoder_util::{decode_u128_2_dna_seq};
use kmer_count::sequence_encoder_util::DnaSequence;
use bio::io::fasta::Reader as faReader;
use bio::io::fasta::Record as faRecord;
use std::fs::File;
use crate::bio::io::fasta::FastaRead;



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
    eprintln!("input  file: {:?}",  input_file);
    eprintln!("loading {:?} done", input_file);


    let file = File::open(&input_file).expect("Error during opening the file");
    let mut reader = faReader::new(file);
    let mut record = faRecord::new();
    let mut sequences: Vec<DnaSequence> = Vec::new();
    eprintln!("loading {:?} done", input_file);
    'each_read: loop {
        reader.read(&mut record).unwrap();
        if record.is_empty(){
            break 'each_read;
        }
        let sequence_as_vec: Vec<u8> = record.seq().to_vec();
        let current_sequence = DnaSequence::new(&sequence_as_vec);
        sequences.push(current_sequence);
    }





/*

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


*/










/*
ここにマルチスレッド処理を書く
*/

/* 
    let sequences_len: usize = sequences.len();
    let chunk_size: usize = sequences_len / threads;
    let mut cbf_oyadama: Box<[u32; BLOOMFILTER_TABLE_SIZE]> = Box::new([0; BLOOMFILTER_TABLE_SIZE]);
    let iter_cnt = Arc::new([1..threads]);
    for i in iter_cnt {
        let handle = thread::spawn(|| {
            let start: usize = i * chunk_size;
            let end: usize;
            if i != threads - 1{
                end = (i + 1) * chunk_size;
            }else{
                end = sequences_len - 1;
            }
            eprintln!("start calling build_counting_bloom_filter[{}]", i);
            let cbf: Box<[u32; BLOOMFILTER_TABLE_SIZE]> = build_counting_bloom_filter(&sequences, start, end);
            eprintln!("finish calling build_counting_bloom_filter");
        });
    }
 */
    eprintln!("start calling build_counting_bloom_filter[{}]", i);
    let cbf: Box<[u32; BLOOMFILTER_TABLE_SIZE]> = build_counting_bloom_filter(&sequences, start, end);
    eprintln!("finish calling build_counting_bloom_filter");

    eprintln!("start calling number_of_high_occurence_kmer");
    let h_cbf_h: HashSet<u128> = number_of_high_occurence_kmer(&cbf_oyadama, &sequences, threshold);
    eprintln!("finish calling number_of_high_occurence_kmer");
    let mut high_occurence_kmer: Vec<u128> = h_cbf_h.into_iter().collect();

/*
ここまで
*/




    //sortする
    eprintln!("start voracious_mt_sort({})", threads);
    high_occurence_kmer.voracious_mt_sort(threads);
    eprintln!("finish voracious_mt_sort({})", threads);

    eprintln!("start  writing to output file: {:?}", &output_file);

    //let mut w = File::create(&output_file).unwrap();
    let mut w = BufWriter::new(fs::File::create(&output_file).unwrap());
    //let mut w_kensho = BufWriter::new(fs::File::create("./kensho_out").unwrap());

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
    eprintln!("L:{}, M:{}, R{}, threshold: {}({}x63)\tcardinarity: {}", L_LEN, M_LEN, R_LEN, threshold, threshold / 63, cnt);
    eprintln!("total cardinarity: {}", cnt);
    eprintln!("threads: {}\tthreshold: {}\tinput file {:?}", threads, threshold, &input_file);

}