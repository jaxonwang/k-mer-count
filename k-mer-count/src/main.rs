extern crate bio;

use std::fs::File;
use bio::io::fasta;
use bio::io::fasta::{Record, FastaRead};
use std::str::from_utf8;

fn main() {
    let file = File::open("sample.fasta").expect("Error during opening the file");
    let mut reader = Box::new(fasta::Reader::new(file));
    let mut record = Record::new();

    let l_len = 27;
    let r_len = 27;
    let mut lr_chunk = Vec::new();
    for dna_chunk_size in 80..141 {
        let m_len = dna_chunk_size - l_len - r_len;
        let mut window_start = 0;
        loop{
            reader.read(&mut record).unwrap();
            if record.is_empty(){
                break;
            }
            let l_start: usize = window_start;
            let l_end:   usize = l_start + l_len;
            let r_start: usize = l_end   + m_len;
            let r_end:   usize = r_start + r_len;
            window_start += 1;
            if r_end > record.seq().len().try_into().unwrap(){
                break;
            }
            let l = &record.seq()[l_start..l_end];
            let r = &record.seq()[r_start..r_end];
/*
            let spacer1: String = " ".repeat(l_start);
            let spacer2: String = " ".repeat(m_len);
            println!("  {:?}\n{:?}{:?}{:?}{:?}",
                from_utf8(&record.seq()).unwrap(),
                spacer1,
                from_utf8(l).unwrap(),
                spacer2,
                from_utf8(r).unwrap()
                );
*/
            let tmp_lr_chunk = from_utf8(l).unwrap().clone().to_owned() + from_utf8(r).unwrap();
            //lr_chunk.push(format!("{:?}{:?}", l, r));
            lr_chunk.push(tmp_lr_chunk);
        }
    }
    lr_chunk.sort();
    for each_chunk in lr_chunk.iter() {
        //println!("{:?}", from_utf8(each_chunk).unwrap());
        println!("{}", each_chunk);
    }
    //println!("{:?}", lr_chunk);
}

/*
    for DNA_chunk_size in range(80, 140 + 1):
        M_len = DNA_chunk_size - L_len - R_len
        for each_read in reads:
            read_length = len(each_read)
            i = 0
            while True:
                L_start = i
                L_end   = L_start + L_len
                R_start = L_end + M_len
                R_end   = R_start + R_len
                L = each_read[L_start:L_end]
                R = each_read[R_start:R_end]
                if len(R) != R_len:
                    break
                #print(f"L_start: {L_start}\tL_end: {L_end}\tlen(L): {len(L)}\tR_start: {R_start}\tR_end: {R_end}\tlen(R): {len(R)}")
                reads_LR.append(L + R)
                i = i + 1
*/