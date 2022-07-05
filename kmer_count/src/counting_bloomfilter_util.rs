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

pub const BLOOMFILTER_TABLE_SIZE: usize = u32::MAX as usize + 1;
pub const THRESHOLD_OCCURENCE: u64 = 100;


pub fn build_counting_bloom_filter(path: &str) -> Box<[u64; BLOOMFILTER_TABLE_SIZE]>{
    let mut l_window_start: usize = 0;
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
    loop {
        reader.read(&mut record).unwrap();
        if record.is_empty(){
            continue;
        }
        let mut add_bloom_filter_cnt: usize = 0;
        eprint!("1st loop: {:09?}, current record id:{:?}\tlength: {:?}\t", loop_cnt, record.id(), record.seq().len());
        loop_cnt += 1;
        //recordをVec<u8>に変更して、DNA_sequence.new()に渡す
        //&[u8] -> Vec<u8>
        let sequence_as_vec: Vec<u8> = record.seq().to_vec();
        let current_sequence = DnaSequence::new(&sequence_as_vec);
        //for dna_chunk_size in 80..141 {
        loop{
            l_window_end = l_window_start + L_LEN;
            if l_window_end >= current_sequence.len(){
                break;
            }
            let l_has_poly_base: bool = current_sequence.has_poly_base(l_window_start, l_window_end);
            if  l_has_poly_base == true{
                l_window_start += 1;
                continue;
            }
            for dna_chunk_size in 80..141 {
                r_window_start = l_window_start + dna_chunk_size - R_LEN;
                r_window_end   = r_window_start + R_LEN;
                if r_window_end >= current_sequence.len(){
                    continue;
                }
                let r_has_poly_base: bool = current_sequence.has_poly_base(r_window_start, r_window_end);
                if r_has_poly_base != true{
                    add_bloom_filter_cnt += 1;
                    //counting bloom_filterに追加する
                    let lr_string = current_sequence.subsequence_as_u128(vec![[l_window_start, l_window_end], [r_window_start, r_window_end]]);
                    let table_indice:[u32;8] = hash_from_u128(lr_string);//u128を受けてhashを返す関数
                    for i in 0..8{
                        if ret_array[table_indice[i] as usize] == u64::MAX{
                            //incrementしない
                        }else{
                            if rng.gen::<u64>() < (u64::MAX >> (64 - ret_array[table_indice[i] as usize].leading_zeros())){
                                ret_array[table_indice[i] as usize] += 1;
                            }
                        }
                    }
                }
            }
            l_window_start += 1;
        }
        let end = start.elapsed();
        eprintln!("sec: {}, subject to add bloom filter: {}", end.as_secs() - previous_time.as_secs(), add_bloom_filter_cnt);
        previous_time = end;
    }
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
/*
pub fn count_occurence_from_counting_bloomfilter_table(counting_bloomfilter_table: &Box<[u64; BLOOMFILTER_TABLE_SIZE]>, query: &[u8;L_LEN + R_LEN]) -> u64{
    let indice: [u32;8] = hasher(query);
    let mut retval: u64 = u64::MAX;
    for index in indice{
        if counting_bloomfilter_table[index as usize] < retval{
            retval = counting_bloomfilter_table[index as usize];
        }
    }
    return retval;
}
*/

/*
pub fn hasher(source: &[u8;L_LEN + R_LEN]) -> [u32;8]{
    let mut ret_val: [u32;8] = [0;8];
    let mut hasher = Sha256::new();
    hasher.update(source);
    let result = hasher.finalize();
    let sha256_bit_array = result.as_slice();//&[u8;32]
    for i in 0..8{
        for j in 0..4{
            ret_val[i] += sha256_bit_array[i * 8 + j] as u32;
            ret_val[i] <<= 8;
        }
    }
    return ret_val;
}
*/
//2週目。出現頻度がある閾値を越えるk-merの個数を返す。
//3週目ではこの個数を受けて、vecのメモリを確保して出力用vecを用意して、再びファイルを舐める。
/*
pub fn number_of_high_occurence_kmer(source_table: &Box<[u64; BLOOMFILTER_TABLE_SIZE]>, path: &str) -> u64{
    let mut retval: u64 = 0;
    let mut window_start: usize;
    let mut l_start: usize;
    let mut l_end:   usize;
    let mut r_start: usize;
    let mut r_end:   usize;
    let mut m_len:   usize;
    let mut loop_cnt:usize = 0;

    let file = File::open(path).expect("Error during opening the file");
    let mut reader = faReader::new(file);
    let mut record = faRecord::new();
    let mut buf: u64 = 0;
    let mut lr_string: [u8;L_LEN + R_LEN] = [64; L_LEN + R_LEN];

    loop {
        reader.read(&mut record).unwrap();
        if record.is_empty(){
            break;
        }
        eprintln!("2nd loop: {:09?}, current record id:{:?}\tlength: {:?}", loop_cnt, record.id(), record.seq().len());
        loop_cnt += 1;
        for dna_chunk_size in 80..141 {
            window_start = 0;
            loop{
                m_len = dna_chunk_size - L_LEN - R_LEN;
                l_start = window_start;
                l_end   = l_start + L_LEN;
                r_start = l_end + m_len;
                r_end   = r_start + R_LEN;
                window_start += 1;

                if r_end > record.seq().len(){
                    break;
                }
                let l = &record.seq()[l_start..l_end];
                let r = &record.seq()[r_start..r_end];
                for i in 0..L_LEN{
                    lr_string[i] = l[i];
                }
                for i in 0..R_LEN{
                    lr_string[i + L_LEN] = r[i];
                }
                let table_indice:[u32;8] = hasher(&lr_string);
                let tmp: u64 = count_occurence_from_counting_bloomfilter_table(&source_table, &lr_string);
                if tmp >= THRESHOLD_OCCURENCE{
                    retval = retval + 1;
                }
            }
        }
    }
    return retval;
}

*/
//3週目。

/*
pub fn pick_up_high_occurence_kmer(source_table: &Box<[u64; BLOOMFILTER_TABLE_SIZE]>, path: &str, max_size_of_vec: u64) -> Vec<u128>{
    let mut retval: Vec<u128> = vec![0; max_size_of_vec.try_into().unwrap()];
    let mut retval_index: usize = 0;
    let mut window_start: usize;
    let mut l_start: usize;
    let mut l_end:   usize;
    let mut r_start: usize;
    let mut r_end:   usize;
    let mut m_len:   usize;
    let mut loop_cnt:usize = 0;

    let file = File::open(path).expect("Error during opening the file");
    let mut reader = faReader::new(file);
    let mut record = faRecord::new();
    let mut buf: u64 = 0;
    let mut lr_string: [u8;L_LEN + R_LEN] = [64; L_LEN + R_LEN];

    loop {
        reader.read(&mut record).unwrap();
        if record.is_empty(){
            break;
        }
        eprintln!("3rd loop: {:09?}, current record id:{:?}\tlength: {:?}", loop_cnt, record.id(), record.seq().len());
        loop_cnt += 1;
        for dna_chunk_size in 80..141 {
            window_start = 0;
            loop{
                m_len = dna_chunk_size - L_LEN - R_LEN;
                l_start = window_start;
                l_end   = l_start + L_LEN;
                r_start = l_end + m_len;
                r_end   = r_start + R_LEN;
                window_start += 1;

                if r_end > record.seq().len(){
                    break;
                }
                let l = &record.seq()[l_start..l_end];
                let r = &record.seq()[r_start..r_end];
                for i in 0..L_LEN{
                    lr_string[i] = l[i];
                }
                for i in 0..R_LEN{
                    lr_string[i + L_LEN] = r[i];
                }
                //ここら辺に、閾値回数以上出現するk-merを処理するコードを書く
                let table_indice:[u32;8] = hasher(&lr_string);
                let tmp: u64 = count_occurence_from_counting_bloomfilter_table(&source_table, &lr_string);
                if tmp >= THRESHOLD_OCCURENCE{//2^10相当。本当は引数で基準を変えられるようにしたい。
                    retval[retval_index] = encode_dna_seq_2_u128(&lr_string);
                    retval_index += 1;
                }
            }
        }
    }
    return retval;
}
*/