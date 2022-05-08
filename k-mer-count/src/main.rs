extern crate bio;
extern crate rdxsort;

use std::fs::File;
use bio::io::fasta::{Record, FastaRead, Reader};
use std::str::from_utf8;


fn bucket_sort(source: Vec<&str>, place: usize) -> Vec<&str>{
    assert!(source[0].len() >= place);
    let mut result_a: Vec<&str> = Vec::new();
    let mut result_c: Vec<&str> = Vec::new();
    let mut result_g: Vec<&str> = Vec::new();
    let mut result_t: Vec<&str> = Vec::new();
    let index:  usize = source[0].len() - place;
    for each_chunk in source.iter(){
        let current_char = each_chunk.chars().nth(index).unwrap();
        match current_char{
            'A' => {result_a.push(each_chunk);}
            'C' => {result_c.push(each_chunk);}
            'G' => {result_g.push(each_chunk);}
            'T' => {result_t.push(each_chunk);}
            _   => {panic!("Unexpected charactor {} appears in {}", current_char, each_chunk);}
        }
    }
    let mut result: Vec<&str> = Vec::new();
    result.append(&mut result_a);
    result.append(&mut result_c);
    result.append(&mut result_g);
    result.append(&mut result_t);
    return result;
}

fn radix_sort(mut source: Vec<&str>) -> Vec<&str>{
    let length: usize = source[0].len();
    for i in 1..length{
        source = bucket_sort(source, i);
    }
    return source;
}


fn main() {
    let file = File::open("sample.fasta").expect("Error during opening the file");
    let mut reader = Reader::new(file);
    let mut record = Record::new();

    let l_len = 27;
    let r_len = 27;
    let mut lr_chunk:Vec<String> = Vec::new();
    let mut window_start: usize;
    let mut l_start: usize;
    let mut l_end:   usize;
    let mut r_start: usize;
    let mut r_end:   usize;
    let mut m_len:   usize;

    loop {
        reader.read(&mut record).unwrap();
        if record.is_empty(){
            break;
        }
        for dna_chunk_size in 80..141 {
            window_start = 0;
            loop{
                m_len = dna_chunk_size - l_len - r_len;
                l_start = window_start;
                l_end   = l_start + l_len;
                r_start = l_end + m_len;
                r_end   = r_start + r_len;
                window_start += 1;

                if r_end > record.seq().len(){
                    break;
                }
                let l = &record.seq()[l_start..l_end];
                let r = &record.seq()[r_start..r_end];
                let tmp_lr_chunk = from_utf8(l).unwrap().clone().to_owned() + from_utf8(r).unwrap();
                lr_chunk.push(tmp_lr_chunk);
            }
        }
    }
    
    let sorted_lr_chunk = radix_sort(lr_chunk.iter().map(|s| &**s).collect());


    lr_chunk.sort();
    for each_chunk in lr_chunk.iter() {
        println!("{}", each_chunk);
    }
    
}