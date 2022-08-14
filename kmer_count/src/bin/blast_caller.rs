extern crate kmer_count;
use kmer_count::sequence_encoder_util::{decode_u128_l, decode_u128_r, decode_u128_2_dna_seq};

/*
blast検索...
1. L-27-merの完全マッチを探す-> L-hit
2. R-27-merの完全マッチを探す-> R-hit
L-hit., R-hitが同じ染色体にある場合、両者の距離を比較する
近いとアニーリングしちゃうかも
（Pythonでよくないか？）
*/


fn main(){}
