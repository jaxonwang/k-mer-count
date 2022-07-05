pub struct DnaSequence{
    sequence: Vec<u64>,
    length:   usize
}

impl DnaSequence{
    pub fn new(source: &Vec<u8>) -> DnaSequence{
        let mut retval:Vec<u64> = Vec::new();
        let mut buf: u64 = 0;
        let mut cnt: u64 = 0;
        for each_base in source.iter(){
            match each_base{
                b'A' => {buf |= 0;}
                b'C' => {buf |= 1;}
                b'G' => {buf |= 2;}
                b'T' => {buf |= 3;}
                _   => {panic!("Unexpected character: {}", each_base);}
            }
            //buf <<= 2;//境界を跨ぐ場合、シフトしてはいけない。
            if cnt == 31{//cnt == 31のときはbufが全部埋まってる。
                retval.push(buf);
                buf = 0;
                cnt = 0;
            }else{
                buf <<= 2;
                cnt += 1;
            }
        }

        if cnt != 32{
            buf <<= 2 * (31 - cnt);//塩基を表すbitを上位に寄せる。
            retval.push(buf);//32塩基の節目で切れなかった時に備えてpushする。
        }
/*
        for each_buf in retval.iter(){
            println!("{:064b}", each_buf);
        }
*/
        return DnaSequence{sequence: retval, length: source.len()};
    }

    pub fn encode(source: &Vec<u8>) -> DnaSequence{
        return DnaSequence::new(source);
    }

    pub fn len(&self) -> usize{
        return self.length;
    }

    pub fn decode(&self) -> Vec<u8>{
        let mut retval = Vec::new();
        let mut buf: u8 = 0;
        for i in 0..self.length{
            buf = ((self.sequence[i / 32] >> (2 * (i % 32))) & 3).try_into().unwrap();
            match buf{
                0 => {retval.push('A' as u8);}
                1 => {retval.push('C' as u8);}
                2 => {retval.push('G' as u8);}
                3 => {retval.push('T' as u8);}
                _ => {panic!("Never reached!!!buf: {}", buf);}
            }
        }
        return retval;
    }

    pub fn subsequence(&self, ranges: Vec<[usize; 2]>) -> DnaSequence{
        let mut buf: u64 = 0;
        let mut cnt: u64 = 0;
        let mut retval: Vec<u64> = Vec::new();
        let mut length: usize = 0;
        for each_range in ranges.iter(){
            let start: usize = each_range[0];
            let end:   usize = each_range[1];
            assert!(start < end, "DnaSequence::subsequence assertion failed: {} !< {}", start, end);
            assert!(start >= 0, "DnaSequence::subsequence assertion failed: {} >= 0", start);
            assert!(end < self.length, "DnaSequence::subsequence assertion failed: {} < {}", end, self.length);
            length += end - start;
            for i in start..end{
                buf += (self.sequence[i / 32] >> (62 - 2 * (i % 32))) & 3;
                buf <<= 2;
                cnt += 1;
                if cnt == 31{
                    retval.push(buf);
                    buf = 0;
                    cnt = 0;
                }
            }
        }
        if cnt != 31{
            buf <<= 2 * (31 - cnt);//塩基を表すbitを上位に寄せる。
            retval.push(buf);//32塩基の節目で切れなかった時に備えてpushする。
        }
        return DnaSequence{sequence: retval, length: length}
    }

    pub fn subsequence_as_u128(&self, ranges: Vec<[usize; 2]>) -> u128{
        let mut buf: u128 = 0;
        let mut cnt: usize = 0;
        for each_range in ranges.iter(){
            let start: usize = each_range[0];
            let end:   usize = each_range[1];
            assert!(start < end, "DnaSequence::subsequence_as_u128 assertion failed: {} !< {}", start, end);
            assert!(start >= 0, "DnaSequence::subsequence_as_u128 assertion failed: {} >= 0", start);
            assert!(end < self.length, "DnaSequence::subsequence_as_u128 assertion failed: {} < {}", end, self.length);
            for i in start..end{
                buf += ((self.sequence[i / 32] >> (62 - 2 * (i % 32))) & 3) as u128;
                buf <<= 2;
                cnt += 1;
            }
            assert!(cnt < 64, "DnaSequence::subsequence_as_u128 assertion failed: too many DNA size. cnt reaches {}", cnt);
        }
        return buf
    }


    pub fn has_poly_base(&self, start: usize, end: usize) -> bool {
        assert!(start < end, "DnaSequence::has_poly_base assertion failed: {} !< {}", start, end);
        assert!(end - start > 3, "DnaSequence::has_poly_base assertion failed: {} - {} < 4", end, start);
        assert!(end - start < 32, "DnaSequence::has_poly_base assertion failed: length of the evaluation subject must be shorter than 32");
        assert!(end < self.length, "DnaSequence::has_poly_base assertion failed: end coordinate must be smaller than length of the sequence. end: {}, self.lngth: {}", end, self.length);
        let mut original:  u64 = 0;
        let mut zero_ichi: u64 = 1;
        for i in start..end{
            original += (self.sequence[i / 32] >> (62 - 2 * (i % 32))) & 3;
            if i != end -1{
                original <<= 2;
                zero_ichi <<= 2;
                zero_ichi += 1;
            }
        }//ここまでで、originalに右詰で対象の領域がコピーされる。
        let val1 = original;
        let val2 = original >> 2;
        let val3 = val1 ^ val2;
        let val4 = val3 >> 1;
        let val5 = val3 ^ val4;
        let val6 = !val5;
        let val7 = val6 & zero_ichi;
        let val8 = val7 >> 2;
        let val9 = val7 >> 4;
/*
        println!("start:    {}", start);
        println!("end:      {}", end);
        println!("0101: {:064b}", zero_ichi);
        println!("val1: {:064b}", val1);
        println!("val2: {:064b}", val2);
        println!("val3: {:064b}", val3);
        println!("val4: {:064b}", val4);
        println!("val5: {:064b}", val5);
        println!("val6: {:064b}", val6);
        println!("val7: {:064b}", val7);
        println!("val8: {:064b}", val8);
        println!("val9: {:064b}", val9);
*/
        return val7 & val8 & val9 != 0;
        //shift演算でポリ塩基の情報がおっこちてる
    }
}


#[cfg(test)]
mod tests{
    use crate::sequence_encoder_util::DnaSequence;
    #[test]
    fn encode_test_4A(){
        let source: Vec<u8> = vec![b'A', b'A', b'A', b'A'];
        let obj = DnaSequence::new(&source);
        assert!(obj.sequence[0] == 0x0000000000000000, "encode_test_1 failed");
    }

    #[test]
    fn encode_test_4C(){
        let source: Vec<u8> = vec![b'C', b'C', b'C', b'C'];
        let obj = DnaSequence::new(&source);
        assert!(obj.sequence[0] == 0x5500000000000000, "encode_test_1 failed");
    }

    #[test]
    fn encode_test_4G(){
        let source: Vec<u8> = vec![b'G', b'G', b'G', b'G'];
        let obj = DnaSequence::new(&source);
        assert!(obj.sequence[0] == 0xAA00000000000000, "encode_test_1 failed");
    }


    #[test]
    fn encode_test_4T(){
        let source: Vec<u8> = vec![b'T', b'T', b'T', b'T'];
        let obj = DnaSequence::new(&source);
        assert!(obj.sequence[0] == 0xFF00000000000000, "encode_test_1 failed");
    }

    #[test]
    fn encode_test_ACGT(){
        let source: Vec<u8> = vec![b'A', b'C', b'G', b'T'];
        let obj = DnaSequence::new(&source);
        assert!(obj.sequence[0] == 0x1b00000000000000, "encode_test_2 failed. {:08b}", obj.sequence[0]);//上位bitに寄せる
    }
    #[test]
    fn encode_test_16C(){
        let source: String = "CCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555500000000, "encode_test_3 failed. index {} {:08b}", 0, obj.sequence[0]);
    }

    #[test]
    fn encode_test_31C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555554, "encode_test_3 failed. index {} {:08b}", 0, obj.sequence[0]);
    }


    #[test]
    fn encode_test_32C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555555, "encode_test_3 failed. index {} {:08x}", 0, obj.sequence[0]);
    }
    #[test]
    fn encode_test_33C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555555, "encode_test_3 failed. index {} {:08x}", 0, obj.sequence[0]);
        assert!(obj.sequence[1] == 0x4000000000000000, "encode_test_3 failed. index {} {:08x}", 1, obj.sequence[1]);
    }

    #[test]
    fn encode_test_34C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555555, "encode_test_3 failed. index {} {:08x}", 0, obj.sequence[0]);
        assert!(obj.sequence[1] == 0x5000000000000000, "encode_test_3 failed. index {} {:08x}", 1, obj.sequence[1]);
    }


    #[test]
    fn encode_test_63C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555555, "encode_test_3 failed. index {} {:08x}", 0, obj.sequence[0]);
        assert!(obj.sequence[1] == 0x5555555555555554, "encode_test_3 failed. index {} {:08x}", 1, obj.sequence[1]);
    }

    #[test]
    fn encode_test_64C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555555, "encode_test_3 failed. index {} {:08x}", 0, obj.sequence[0]);
        assert!(obj.sequence[1] == 0x5555555555555555, "encode_test_3 failed. index {} {:08x}", 1, obj.sequence[1]);
    }

    #[test]
    fn encode_test_65C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555555, "encode_test_3 failed. index {} {:08x}", 0, obj.sequence[0]);
        assert!(obj.sequence[1] == 0x5555555555555555, "encode_test_3 failed. index {} {:08x}", 1, obj.sequence[1]);
        assert!(obj.sequence[2] == 0x4000000000000000, "encode_test_3 failed. index {} {:08x}", 2, obj.sequence[2]);
    }


    #[test]
    fn encode_test_32N(){
        let source: String = "GAACGACTGTTTCTACTATAAATCCTTCCTTC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0b1000000110000111101111110111000111001100000011010111110101111101, "encode_test_3 failed. index {} {:08b}", 0, obj.sequence[0]);
    }//                               G A A C G A C T G T T T C T A C T A T A A A T C C T T C C T T C

    #[test]
    fn encode_test_120N(){
        let source: String = "GAACGACTGTTTCTACTATAAATCCTTCCTTCCTAGCCTATCATTTCTGGAGTCCTTGGTGAACTGTAGGAAGCTCTGAACACACACGTTCCCTTGGATTCGTACCTATGAATACTCCGT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x8187bf71cc0d7d7d, "encode_test_3 failed. index {} {:08b}", 0, obj.sequence[0]);
        assert!(obj.sequence[1] == 0x725cd3f7a2d7eb81, "encode_test_3 failed. index {} {:08b}", 1, obj.sequence[1]);
        assert!(obj.sequence[2] == 0xeca09de04446f57e, "encode_test_3 failed. index {} {:08b}", 2, obj.sequence[2]);
        assert!(obj.sequence[3] == 0x8f6c5ce0c75b0000, "encode_test_3 failed. index {} {:08b}", 3, obj.sequence[3]);
    }

    #[test]
    fn has_poly_base_test_8C(){
        let source: Vec<u8> = vec![b'C', b'C', b'C', b'C', b'C', b'C', b'C', b'C'];
        let obj = DnaSequence::new(&source);
        assert!(obj.has_poly_base(0, 4) == true, "has_poly_base_test_8C failed: obj.has_poly_base(0, 4) returns {}", obj.has_poly_base(0, 4));
        assert!(obj.has_poly_base(0, 5) == true, "has_poly_base_test_8C failed: obj.has_poly_base(0, 5) returns {}", obj.has_poly_base(0, 5));
        assert!(obj.has_poly_base(0, 6) == true, "has_poly_base_test_8C failed: obj.has_poly_base(0, 6) returns {}", obj.has_poly_base(0, 6));
        assert!(obj.has_poly_base(0, 7) == true, "has_poly_base_test_8C failed: obj.has_poly_base(0, 7) returns {}", obj.has_poly_base(0, 7));
    }

    #[test]
    fn has_poly_base_test_8G(){
        let source: Vec<u8> = vec![b'G', b'G', b'G', b'G', b'G', b'G', b'G', b'G'];
        let obj = DnaSequence::new(&source);
        assert!(obj.has_poly_base(0, 4) == true, "has_poly_base_test_8G failed: obj.has_poly_base(0, 4) returns {}", obj.has_poly_base(0, 4));
        assert!(obj.has_poly_base(0, 5) == true, "has_poly_base_test_8G failed: obj.has_poly_base(0, 5) returns {}", obj.has_poly_base(0, 5));
        assert!(obj.has_poly_base(0, 6) == true, "has_poly_base_test_8G failed: obj.has_poly_base(0, 6) returns {}", obj.has_poly_base(0, 6));
        assert!(obj.has_poly_base(0, 7) == true, "has_poly_base_test_8G failed: obj.has_poly_base(0, 7) returns {}", obj.has_poly_base(0, 7));
    }


    #[test]
    fn has_poly_base_test_8T(){
        let source: Vec<u8> = vec![b'T', b'T', b'T', b'T', b'T', b'T', b'T', b'T'];
        let obj = DnaSequence::new(&source);
        assert!(obj.has_poly_base(0, 4) == true, "has_poly_base_test_8T failed: obj.has_poly_base(0, 4) returns {}", obj.has_poly_base(0, 4));
        assert!(obj.has_poly_base(0, 5) == true, "has_poly_base_test_8T failed: obj.has_poly_base(0, 5) returns {}", obj.has_poly_base(0, 5));
        assert!(obj.has_poly_base(0, 6) == true, "has_poly_base_test_8T failed: obj.has_poly_base(0, 6) returns {}", obj.has_poly_base(0, 6));
        assert!(obj.has_poly_base(0, 7) == true, "has_poly_base_test_8T failed: obj.has_poly_base(0, 7) returns {}", obj.has_poly_base(0, 7));
    }


    #[test]
    fn has_poly_base_test_8N(){
        let source: Vec<u8> = vec![b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T'];
        let obj = DnaSequence::new(&source);
        assert!(obj.has_poly_base(0, 4) == false, "has_poly_base_test_8N failed: obj.has_poly_base(0, 4) returns {}", obj.has_poly_base(0, 4));
        assert!(obj.has_poly_base(0, 5) == false, "has_poly_base_test_8N failed: obj.has_poly_base(0, 5) returns {}", obj.has_poly_base(0, 5));
        assert!(obj.has_poly_base(0, 6) == false, "has_poly_base_test_8N failed: obj.has_poly_base(0, 6) returns {}", obj.has_poly_base(0, 6));
        assert!(obj.has_poly_base(0, 7) == false, "has_poly_base_test_8N failed: obj.has_poly_base(0, 7) returns {}", obj.has_poly_base(0, 7));
    }

    #[test]
    fn has_poly_base_test_120N(){
        let source: String = "GAACGACTGTTTTTACTATAAATCCTTCCTTCCTAGCCTATCATTTCTGGAGTCCTTGGTGAACTGTAGGAAGCTCTGAACACACACGTTCCCTTGGATTCGTACCTATGAATACTCCGT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 9)   == false, "has_poly_base_test_120N failed: obj.has_poly_base (0, 9)   returns {}", obj.has_poly_base(0, 4));
        assert!(obj.has_poly_base(0, 19)  == true,  "has_poly_base_test_120N failed: obj.has_poly_base (0, 19)  returns {}", obj.has_poly_base(0, 5));
        assert!(obj.has_poly_base(9, 20)  == true,  "has_poly_base_test_120N failed: obj.has_poly_base (9, 20)  returns {}", obj.has_poly_base(0, 6));
        assert!(obj.has_poly_base(28, 39) == false, "has_poly_base_test_120N failed: obj.has_poly_base (28, 39) returns {}", obj.has_poly_base(0, 7));
    }




}

