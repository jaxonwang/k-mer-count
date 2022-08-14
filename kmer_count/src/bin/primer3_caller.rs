extern crate kmer_count;
extern crate getopts;
use std::{env, process};
use std::fs::File;
use std::io::{Read, BufReader};
use getopts::Options;
use kmer_count::sequence_encoder_util::{decode_u128_l, decode_u128_r, decode_u128_2_dna_seq};

/*
1. バイナリファイルを読み込んで、128bitごとに切る。
2. primer3_coreに投げられるように、文字列を生成する。
3. primer3_coreに投げる。

投げる際、上手に並列にやりたい。


*/

fn print_usage(program: &str, opts: &Options) {
	let brief = format!("Usage: {} FILE [options]", program);
	print!("{}", opts.usage(&brief));
	process::exit(0);
}


fn main(){
	let args: Vec<String> = env::args().collect();
	let program = args[0].clone();

	let mut opts = Options::new();
	opts.optopt("o", "output", "set output file name", "NAME");
	opts.optopt("t", "thread", "number of threads to use for radix sort. default value is 8.", "THREAD");
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
	let f: File = File::open(&input_file).unwrap();
	let mut reader = BufReader::new(f);
	let mut buf: [u8; 16] = [0; 16];
	let mut tmp_seq_as_u128: u128 = 0;
	loop {
		match reader.read(&mut buf).unwrap() {
			0 => break,
			n => {
				let buf = &buf[..n];
				for i in 0..16 {
					tmp_seq_as_u128 <<= 8;
					tmp_seq_as_u128 += u128::from(buf[i]);
				}
				//println!("{:#0130b}", tmp_seq_as_u128);
				let l_u8_array = decode_u128_l(&tmp_seq_as_u128);
				let r_u8_array = decode_u128_r(&tmp_seq_as_u128);
				let l_str: &str = std::str::from_utf8(&l_u8_array).unwrap();
				let r_str: &str = std::str::from_utf8(&r_u8_array).unwrap();
				print!("{:?}", l_str);
				println!("{:?}", r_str);
			}
		}
	}
}