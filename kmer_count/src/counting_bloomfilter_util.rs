//use crate::encoder_util::decode_u128_2_dna_seq;
//use crate::encoder_util::encode_dna_seq_2_u128;
//use crate::encoder_util::decode_u128_2_occurence;
pub const L_LEN: usize = 27;
pub const R_LEN: usize = 27;

use crate::sequence_encoder_util::DnaSequence;
use rand::Rng;
use std::fs::File;
use sha2::Sha256;
use sha2::Digest;

use std::time::{/*Duration, */Instant};

//use bio::io::fastq::Reader as fqReader;
//use bio::io::fastq::Record as fqRecord;
use bio::io::fasta::Reader as faReader;
use bio::io::fasta::Record as faRecord;
//use bio::io::fastq::FastqRead;
use bio::io::fasta::FastaRead;


use std::io::{self, Read, Write, BufReader, BufWriter};


pub const BLOOMFILTER_TABLE_SIZE: usize = u32::MAX as usize + 1;
pub const THRESHOLD_OCCURENCE: u64 = 100;



//全てのL, Rと、hash値を出力する
//部分配列のdecoderを書き、テストする
pub fn build_counting_bloom_filter(path: &str) -> Box<[u64; BLOOMFILTER_TABLE_SIZE]>{
    let mut l_window_start: usize;
    let mut l_window_end:   usize;
    let mut r_window_start: usize;
    let mut r_window_end:   usize;

    let mut loop_cnt:usize = 0;
    eprintln!("Allocating Box<[u64; BLOOMFILTER_TABLE_SIZE]> where BLOOMFILTER_TABLE_SIZE = {}", BLOOMFILTER_TABLE_SIZE);
    let mut ret_array: Box<[u64; BLOOMFILTER_TABLE_SIZE]> = Box::new([0; BLOOMFILTER_TABLE_SIZE]);
    eprintln!("finish allocating");
    let file = File::open(path).expect("Error during opening the file");

    let mut reader = faReader::new(file);
    let mut record = faRecord::new();

    let mut rng = rand::thread_rng();
    let start = Instant::now();
    let mut previous_time = start.elapsed();
    'each_read: loop {
        reader.read(&mut record).unwrap();
        if record.is_empty(){
            break 'each_read;
        }
        l_window_start = 0;
        let mut add_bloom_filter_cnt: usize = 0;
        let mut l_window_cnt: usize = 0;

        eprint!("1st loop: {:09?}, current record id:{:?}\tlength: {:?}\t", loop_cnt, record.id(), record.seq().len());
        loop_cnt += 1;
        let sequence_as_vec: Vec<u8> = record.seq().to_vec();
        let current_sequence = DnaSequence::new(&sequence_as_vec);
        'each_l_window: loop{
            l_window_end = l_window_start + L_LEN;
            if l_window_end >= current_sequence.len(){
                break 'each_l_window;
            }
            l_window_cnt += 1;
            let l_has_poly_base_or_simple_repeat: bool = current_sequence.has_poly_base_or_simple_repeat(l_window_start, l_window_end);
            if  l_has_poly_base_or_simple_repeat == true{
                l_window_start += L_LEN - 4;
                continue 'each_l_window;
            }
            'each_r_window: for dna_chunk_size in 80..141 {
                r_window_start = l_window_start + dna_chunk_size - R_LEN;
                r_window_end   = r_window_start + R_LEN;
                if r_window_end >= current_sequence.len(){
                    break 'each_r_window;
                }
                let r_has_poly_base_or_simple_repeat: bool = current_sequence.has_poly_base_or_simple_repeat(r_window_start, r_window_end);
                if r_has_poly_base_or_simple_repeat == true{
                    r_window_start += R_LEN - 4;
                    continue 'each_r_window;
                }else{
                    add_bloom_filter_cnt += 1;
                    let lr_string = current_sequence.subsequence_as_u128(vec![[l_window_start, l_window_end], [r_window_start, r_window_end]]);
                    let table_indice:[u32;8] = hash_from_u128(lr_string);//u128を受けてhashを返す関数
                    let occurence = count_occurence_from_counting_bloomfilter_table(&ret_array, table_indice);
                    if rng.gen::<u64>() < (u64::MAX >> (64 - occurence.leading_zeros())){
                        for i in 0..8{
                            let idx: usize = table_indice[i] as usize;
                            if ret_array[idx] == u64::MAX{
                                //incrementしない
                            }else{
                                ret_array[idx] += 1;
                            }
                        }
                    }
                }
            }
            l_window_start += 1;
        }
        let end = start.elapsed();
        eprintln!("sec: {}.{:03}\t subject to add bloom filter: {}\tl_window_cnt: {}", end.as_secs() - previous_time.as_secs(),end.subsec_nanos() - previous_time.subsec_nanos(),  add_bloom_filter_cnt, l_window_cnt);
        previous_time = end;
    }
    return ret_array;
}

fn hash_from_u128(source: u128) -> [u32; 8]{
    let mut ret_val: [u32;8] = [0;8];
    let mut hasher = Sha256::new();
    let mut u8_array: [u8; 16] = [0; 16];
    let mut src_copy: u128 = source;
    for i in 0..16{
        u8_array[i] = (src_copy & 255).try_into().unwrap();
        src_copy >>= 8;
    }
    hasher.update(u8_array);
    let result = hasher.finalize();
    let sha256_bit_array = result.as_slice();//&[u8;32]
    for i in 0..8{
        for j in 0..4{
            ret_val[i] += sha256_bit_array[i * 4 + j] as u32;
            ret_val[i] <<= 8;
        }
    }
    return ret_val;

}

fn count_occurence_from_counting_bloomfilter_table(counting_bloomfilter_table: &Box<[u64; BLOOMFILTER_TABLE_SIZE]>, indice: [u32; 8]) -> u64{
    let mut retval: u64 = u64::MAX;
    for index in indice{
        if counting_bloomfilter_table[index as usize] < retval{
            retval = counting_bloomfilter_table[index as usize];
        }
    }
    return retval;
}


pub fn number_of_high_occurence_kmer(source_table: &Box<[u64; BLOOMFILTER_TABLE_SIZE]>, path: &str, threshold: u64) -> (Box<[bool; BLOOMFILTER_TABLE_SIZE]>, usize){
    let mut ret_table: Box<[bool; BLOOMFILTER_TABLE_SIZE]> = Box::new([false; BLOOMFILTER_TABLE_SIZE]);
    let mut ret_val: usize = 0;
    let mut l_window_start: usize;
    let mut l_window_end:   usize;
    let mut r_window_start: usize;
    let mut r_window_end:   usize;
    let mut loop_cnt:usize = 0;
    let file = File::open(path).expect("Error during opening the file");

    let mut reader = faReader::new(file);
    let mut record = faRecord::new();

    let start = Instant::now();
    let mut previous_time = start.elapsed();
    'each_read: loop {
        reader.read(&mut record).unwrap();
        if record.is_empty(){
            break 'each_read;
        }
        l_window_start = 0;

        eprint!("2nd loop: {:09?}, current record id:{:?}\tlength: {:?}\t", loop_cnt, record.id(), record.seq().len());
        loop_cnt += 1;
        let sequence_as_vec: Vec<u8> = record.seq().to_vec();
        let current_sequence = DnaSequence::new(&sequence_as_vec);
        'each_l_window: loop{
            l_window_end = l_window_start + L_LEN;
            if l_window_end >= current_sequence.len(){
                break 'each_l_window;
            }
            let l_has_poly_base_or_simple_repeat: bool = current_sequence.has_poly_base_or_simple_repeat(l_window_start, l_window_end);
            if  l_has_poly_base_or_simple_repeat == true{
                l_window_start += 1;
                continue 'each_l_window;
            }
            'each_r_window: for dna_chunk_size in 80..141 {
                r_window_start = l_window_start + dna_chunk_size - R_LEN;
                r_window_end   = r_window_start + R_LEN;
                if r_window_end >= current_sequence.len(){
                    break 'each_r_window;
                }
                let r_has_poly_base_or_simple_repeat: bool = current_sequence.has_poly_base_or_simple_repeat(r_window_start, r_window_end);
                if r_has_poly_base_or_simple_repeat != true{
                    let lr_string: u128 = current_sequence.subsequence_as_u128(vec![[l_window_start, l_window_end], [r_window_start, r_window_end]]);
                    let table_indice:[u32;8] = hash_from_u128(lr_string);
                    let occurence: u64 = count_occurence_from_counting_bloomfilter_table(source_table, table_indice);
                    if occurence >= threshold{
                        ret_val += 1;
                        for i in 0..8{
                            let idx: usize = table_indice[i] as usize;
                            ret_table[idx] |= true;
                        }
                    }
                }else{//ポリ塩基を持ってるとき
                    continue 'each_r_window;
                }
            }
            l_window_start += 1;
        }
        let end = start.elapsed();
        eprintln!("sec: {}.{:03}", end.as_secs() - previous_time.as_secs(),end.subsec_nanos() - previous_time.subsec_nanos());
        previous_time = end;
    }
    return (ret_table, ret_val);
}



fn refer_bloom_filter_table(bloomfilter_table: &Box<[bool; BLOOMFILTER_TABLE_SIZE]>, query: u128) -> bool {
    let indice: [u32;8] = hash_from_u128(query);
    let mut retval: bool = true;
    for index in indice{
        retval &= bloomfilter_table[index as usize]
    }
    return retval;
}


//3週目。
/*
ファイルを舐めて、個数を確定させてからVecを確保する。
*/

pub fn pick_up_high_occurence_kmer(source_table: &Box<[bool; BLOOMFILTER_TABLE_SIZE]>, path: &str, max_size_of_list: usize) -> Vec<u128>{
    //let mut ret_table: Box<[bool; BLOOMFILTER_TABLE_SIZE]> = Box::new([false; BLOOMFILTER_TABLE_SIZE]);
    let mut l_window_start: usize;
    let mut l_window_end:   usize;
    let mut r_window_start: usize;
    let mut r_window_end:   usize;
    let mut loop_cnt:usize = 0;
    let file = File::open(path).expect("Error during opening the file");

    let mut ret_vec = vec![0 as u128; max_size_of_list + 1];
    let mut ret_vec_cnt: usize = 0;

    let mut reader = faReader::new(file);
    let mut record = faRecord::new();

    let start = Instant::now();
    let mut previous_time = start.elapsed();
    'each_read: loop {
        reader.read(&mut record).unwrap();
        if record.is_empty(){
            break 'each_read;
        }
        l_window_start = 0;

        eprint!("3rd loop: {:09?}, current record id:{:?}\tlength: {:?}\t", loop_cnt, record.id(), record.seq().len());
        loop_cnt += 1;
        let sequence_as_vec: Vec<u8> = record.seq().to_vec();
        let current_sequence = DnaSequence::new(&sequence_as_vec);
        'each_l_window: loop{
            l_window_end = l_window_start + L_LEN;
            if l_window_end >= current_sequence.len(){
                break 'each_l_window;
            }
            let l_has_poly_base_or_simple_repeat: bool = current_sequence.has_poly_base_or_simple_repeat(l_window_start, l_window_end);
            if  l_has_poly_base_or_simple_repeat == true{
                l_window_start += 1;
                continue 'each_l_window;
            }
            'each_r_window: for dna_chunk_size in 80..141 {
                r_window_start = l_window_start + dna_chunk_size - R_LEN;
                r_window_end   = r_window_start + R_LEN;
                if r_window_end >= current_sequence.len(){
                    break 'each_r_window;
                }
                let r_has_poly_base_or_simple_repeat: bool = current_sequence.has_poly_base_or_simple_repeat(r_window_start, r_window_end);
                if r_has_poly_base_or_simple_repeat != true{
                    let lr_string: u128         = current_sequence.subsequence_as_u128(vec![[l_window_start, l_window_end], [r_window_start, r_window_end]]);
                    //let table_indice:[u32;8]    = hash_from_u128(lr_string);
                    let is_high_occr_kmer: bool = refer_bloom_filter_table(source_table, lr_string);
                    if is_high_occr_kmer == true{
                        ret_vec[ret_vec_cnt] = lr_string;
                        ret_vec_cnt += 1;
                    }
                }else{//ポリ塩基を持ってるとき
                    continue 'each_r_window;
                }
            }
            l_window_start += 1;
        }
        let end = start.elapsed();
        eprintln!("sec: {}.{:03}", end.as_secs() - previous_time.as_secs(),end.subsec_nanos() - previous_time.subsec_nanos());
        previous_time = end;
    }
    return ret_vec;
}

/*
pub fn write_to_file_as_binary(sorted_source: Vec<u128>, path: &str) -> (){
    let mut writer = BufWriter::new(File::create(path));
    let mut previous_kmer: u128 = 0;
    for each_kmer in sorted_source{
        if previous_kmer != each_kmer{
            writer.write_all(each_kmer);
            //println!("{:?}\t{:0128b}", String::from_utf8(decode_u128_2_dna_seq(&each_kmer, 54)).unwrap(), each_kmer);
        }
        previous_kmer = each_kmer;
    }
}
*/










