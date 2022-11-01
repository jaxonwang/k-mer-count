use crate::counting_bloomfilter_util::L_LEN;
use crate::counting_bloomfilter_util::M_LEN;
use crate::counting_bloomfilter_util::R_LEN;

pub fn decode_u128_2_dna_seq(source:&u128, char_size: usize) -> Vec<u8>{
    let mut result: Vec<u8> = Vec::new();
    let mut base;
    for i in 0..char_size{
        base = source >> 2 * (char_size - 1 - i) & 3;
        match base{
            0 => {result.push(b'A');}
            1 => {result.push(b'C');}
            2 => {result.push(b'G');}
            3 => {result.push(b'T');}
            _ => {panic!("Never reached!!!base: {}", base);}
        }
    }
    return result;
}

pub fn decode_u128_l(source: &u128) -> [u8; L_LEN]{
    let mut result: [u8; L_LEN] = [b'X'; L_LEN];
    let mut base;
    for i in 0..L_LEN{
        base = source >> (((L_LEN + M_LEN + R_LEN) - i - 1) * 2) & 3;
        match base{
            0 => {result[i] = b'A';}
            1 => {result[i] = b'C';}
            2 => {result[i] = b'G';}
            3 => {result[i] = b'T';}
            _ => {panic!("Never reached!!!base: {}", base);}
        }
    }
    return result;
}

pub fn decode_u128_m(source: &u128) -> [u8; M_LEN]{
    let mut result: [u8; M_LEN] = [b'X'; M_LEN];
    let mut base;
    for i in 0..M_LEN{
        base = source >> (((M_LEN + R_LEN) - i - 1) * 2) & 3;
        match base{
            0 => {result[i] = b'A';}
            1 => {result[i] = b'C';}
            2 => {result[i] = b'G';}
            3 => {result[i] = b'T';}
            _ => {panic!("Never reached!!!base: {}", base);}
        }
    }
    return result;
}

pub fn decode_u128_r(source: &u128) -> [u8; R_LEN]{
    let mut result: [u8; R_LEN] = [b'X'; R_LEN];
    let mut base;
    for i in 0..R_LEN{
        base = source >> ((R_LEN - i - 1) * 2) & 3;
        match base{
            0 => {result[i] = b'A';}
            1 => {result[i] = b'C';}
            2 => {result[i] = b'G';}
            3 => {result[i] = b'T';}
            _ => {panic!("Never reached!!!base: {}", base);}
        }
    }
    return result;
}




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

    pub fn decode(&self, start: usize, end: usize) -> Vec<u8>{
        assert!(start < end, "DnaSequence::decode assertion failed: {} !< {}", start, end);
        assert!(end <= self.length, "DnaSequence::decode assertion failed: {} !< {}", end, self.length);
        let mut retval = Vec::new();
        let mut buf: u8;
        for i in start..end{
            buf = ((self.sequence[i / 32] >> (2 * (31 - i % 32))) & 3).try_into().unwrap();
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
            //assert!(start >= 0, "DnaSequence::subsequence assertion failed: {} >= 0", start);
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

//subsequence_as_u128は右詰め
    pub fn subsequence_as_u128(&self, ranges: Vec<[usize; 2]>) -> u128{
        let mut buf: u128 = 0;
        let mut cnt: usize = 0;
        for each_range in ranges.iter(){
            let start: usize = each_range[0];
            let end:   usize = each_range[1];
            assert!(start < end, "DnaSequence::subsequence_as_u128 assertion failed: {} !< {}", start, end);
            //assert!(start >= 0, "DnaSequence::subsequence_as_u128 assertion failed: {} >= 0", start);
            assert!(end <= self.length, "DnaSequence::subsequence_as_u128 assertion failed: {} < {}", end, self.length);
            for i in start..end{
                buf <<= 2;
                cnt += 1;
                buf += ((self.sequence[i / 32] >> (62 - 2 * (i % 32))) & 3) as u128;
            }
            assert!(cnt <= 64, "DnaSequence::subsequence_as_u128 assertion failed: too many DNA size. cnt reaches {}", cnt);
        }
        return buf
    }

    pub fn has_poly_base_or_simple_repeat(&self, start: usize, end: usize) -> bool {
        return self.has_poly_base(start, end, false).0 | self.has_simple_repeat(start, end).0 | self.has_2base_repeat(start, end).0;
    }

    pub fn has_poly_base(&self, start: usize, end: usize, flag: bool) -> (bool, usize) {
        if flag {
            eprintln!("has_poly_base recieved start: {} end: {} nd self.length is {}", start, end, self.length);
            assert!(end <= self.length, "DnaSequence::has_poly_base assertion failed: end coordinate must be smaller than length of the sequence. start: {}, end: {}, self.lngth: {}", start, end, self.length);
            eprintln!("has_poly_base recieved start: {} end: {} nd self.length is {}", start, end, self.length);
        }
        assert!(start < end, "DnaSequence::has_poly_base assertion failed: {} !< {}", start, end);
        assert!(end - start > 3, "DnaSequence::has_poly_base assertion failed: {} - {} < 4", end, start);
        assert!(end - start < 32, "DnaSequence::has_poly_base assertion failed: length of the evaluation subject must be shorter than 32");
        assert!(end <= self.length, "DnaSequence::has_poly_base assertion failed: end coordinate must be smaller than length of the sequence. start: {}, end: {}, self.lngth: {}", start, end, self.length);
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
        let val1  = original;
        let val2  = original << 2;
        let val3  = val1 ^ val2;
        let val4  = val3 >> 1;
        let val5  = val3 | val4;
        let val6  = !val5;
        let val7  = val6 & zero_ichi;
        let val8  = val7 << 2;
        let val9  = val7 & val8;
        let val10 = (val9 << (2 * (32 + start - end))).leading_zeros() / 2;

        #[cfg(test)]{
            println!("{}", std::str::from_utf8(&self.decode(start, end)).unwrap());
            println!("start:    {}", start);
            println!("end:      {}", end);
            println!(" 0101: {:064b}", zero_ichi);
            println!(" val1: {:064b}", val1);
            println!(" val2: {:064b}", val2);
            println!(" val3: {:064b}", val3);
            println!(" val4: {:064b}", val4);
            println!(" val5: {:064b}", val5);
            println!(" val6: {:064b}", val6);
            println!(" val7: {:064b}", val7);
            println!(" val8: {:064b}", val8);
            println!(" val9: {:064b}", val9);
            println!("val10: {}", val10);
        }
        if val10 == 32{
            return (false, 0)
        }else{
            return (true, val10.try_into().unwrap())
        }
        //shift演算でポリ塩基の情報がおっこちてる
    }

    pub fn has_simple_repeat(&self, start: usize, end: usize) -> (bool, usize) {
        assert!(start < end, "DnaSequence::has_simple_repeat assertion failed: {} !< {}", start, end);
        assert!(end - start > 3, "DnaSequence::has_simple_repeat assertion failed: {} - {} < 4", end, start);
        assert!(end - start < 32, "DnaSequence::has_simple_repeat assertion failed: length of the evaluation subject must be shorter than 32");
        assert!(end <= self.length, "DnaSequence::has_simple_repeat assertion failed: end coordinate must be smaller than length of the sequence. end: {}, self.lngth: {}", end, self.length);
        let mut original:  u64 = 0;
        let mut zero_ichi: u64 = 1;
        for i in start..end{
            original += (self.sequence[i / 32] >> (62 - 2 * (i % 32))) & 3;
            if i != end - 1{
                original <<= 2;
                zero_ichi <<= 2;
                zero_ichi += 1;
            }
        }//ここまでで、originalに右詰で対象の領域がコピーされる。
        let val1 = original;
        let val2 = original << 6;
        let val3 = val1 ^ val2;
        let val4 = val3 >> 1;
        let val5 = val3 | val4;
        let val6 = !val5;
        let val7 = val6 & zero_ichi;
        let val8 = val7 << 2;
        let val9 = val7 << 4;
        let val10 = val7 << 6;
        let val11 = val7 << 8;
        let val12 = val7 << 10;
        let last  = val7 & val8 & val9 & val10 & val11 & val12;
        let leading0 = (last << (2 * (32 + start - end))).leading_zeros() / 2;

        #[cfg(test)]{
            println!("start: {}", start);
            println!("end:   {}", end);
            println!("0101:  {:064b}", zero_ichi);
            println!("val1:  {:064b}", val1);
            println!("val2:  {:064b}", val2);
            println!("val3:  {:064b}", val3);
            println!("val4:  {:064b}", val4);
            println!("val5:  {:064b}", val5);
            println!("val6:  {:064b}", val6);
            println!("val7:  {:064b}", val7);
            println!("val8:  {:064b}", val8);
            println!("val9:  {:064b}", val9);
            println!("val10: {:064b}", val10);
            println!("val11: {:064b}", val11);
            println!("val12: {:064b}", val12);
            println!("last:  {:064b}", last);
        }
        if leading0 == 32{
            return (false, 0)
        }else{
            return (true, leading0.try_into().unwrap())
        }
    }

    pub fn has_2base_repeat(&self, start: usize, end: usize) -> (bool, usize) {
        assert!(start < end, "DnaSequence::has_2base_repeat assertion failed: {} !< {}", start, end);
        assert!(end - start > 3, "DnaSequence::has_2base_repeat assertion failed: {} - {} < 4", end, start);
        assert!(end - start < 32, "DnaSequence::has_2base_repeat assertion failed: length of the evaluation subject must be shorter than 32");
        assert!(end <= self.length, "DnaSequence::has_2base_repeat assertion failed: end coordinate must be smaller than length of the sequence. end: {}, self.lngth: {}", end, self.length);
        let mut original:  u64 = 0;
        //let mut zero_ichi: u64 = 1;
        for i in start..end{
            original += (self.sequence[i / 32] >> (62 - 2 * (i % 32))) & 3;
            if i != end - 1{
                original <<= 2;
                //zero_ichi <<= 2;
                //zero_ichi += 1;
            }
        }//ここまでで、originalに右詰で対象の領域がコピーされる。
        let val1 = original;
        let val2 = original >> 4;
        let mut val3 = val1 ^ val2;//val3で0が20個並んでるのを検出したい。これで2x6の単調反復を検出できる。
        //２個隣の塩基が自身と異なればnon-zero, 同じなら00が立つ
        //意味のあるbitは下位(end - start) * 2bit
        let mut ret_flag = false;
        val3 <<= 2 * (32 - (end - start));
        #[cfg(test)]{
            println!("{}", std::str::from_utf8(&self.decode(start, end)).unwrap());
            println!("start: {}", start);
            println!("end:   {}", end);
            //println!("0101:  {:064b}", zero_ichi);
            println!("val1:  {:064b}", val1);
            println!("val2:  {:064b}", val2);
            println!("val3:  {:064b}", val3);
        }
        let mut leading0: i32 = 0;
        for _i in 0..(end - start){
            #[cfg(test)]{
                println!("val3:  {:064b}", val3);
                println!("val3.leading_zeros(): {}", val3.leading_zeros());//上位のbitから見る。
            }
            if val3.leading_zeros() >= 20{
                ret_flag = true;
                break;
            }
            leading0 = _i.try_into().unwrap();
            val3 <<= 2;
        }
        #[cfg(test)]{
            println!("leading0: {}", leading0);
        }

        if leading0 == (end - start - 1).try_into().unwrap(){
            return (false, 0)
        }else{
            return (true, leading0.try_into().unwrap())
        }
    }
}


fn get_function_name<F>(_: F) -> &'static str
where
    F: Fn(),
{
    std::any::type_name::<F>()
}


#[cfg(test)]
mod tests{
    use crate::sequence_encoder_util::DnaSequence;
    use ::function_name::named;
//Encode test
    #[test]
    #[named]
    fn encode_test_4A(){
        let source: Vec<u8> = vec![b'A', b'A', b'A', b'A'];
        let obj = DnaSequence::new(&source);
        assert!(obj.sequence[0] == 0x0000000000000000, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_4C(){
        let source: Vec<u8> = vec![b'C', b'C', b'C', b'C'];
        let obj = DnaSequence::new(&source);
        assert!(obj.sequence[0] == 0x5500000000000000, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_4G(){
        let source: Vec<u8> = vec![b'G', b'G', b'G', b'G'];
        let obj = DnaSequence::new(&source);
        assert!(obj.sequence[0] == 0xAA00000000000000, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_4T(){
        let source: Vec<u8> = vec![b'T', b'T', b'T', b'T'];
        let obj = DnaSequence::new(&source);
        assert!(obj.sequence[0] == 0xFF00000000000000, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_ACGT(){
        let source: Vec<u8> = vec![b'A', b'C', b'G', b'T'];
        let obj = DnaSequence::new(&source);
        assert!(obj.sequence[0] == 0x1b00000000000000, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_16C(){
        let source: String = "CCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555500000000, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_31C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555554, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_32C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555555, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_33C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555555, "{} failed", function_name!());
        assert!(obj.sequence[1] == 0x4000000000000000, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_34C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555555, "{} failed", function_name!());
        assert!(obj.sequence[1] == 0x5000000000000000, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_63C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555555, "{} failed", function_name!());
        assert!(obj.sequence[1] == 0x5555555555555554, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_64C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555555, "{} failed", function_name!());
        assert!(obj.sequence[1] == 0x5555555555555555, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_65C(){
        let source: String = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x5555555555555555, "{} failed", function_name!());
        assert!(obj.sequence[1] == 0x5555555555555555, "{} failed", function_name!());
        assert!(obj.sequence[2] == 0x4000000000000000, "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn encode_test_32N(){
        let source: String = "GAACGACTGTTTCTACTATAAATCCTTCCTTC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0b1000000110000111101111110111000111001100000011010111110101111101, "{} failed", function_name!());
    }//                               G A A C G A C T G T T T C T A C T A T A A A T C C T T C C T T C
    #[test]
    #[named]
    fn encode_test_120N(){
        let source: String = "GAACGACTGTTTCTACTATAAATCCTTCCTTCCTAGCCTATCATTTCTGGAGTCCTTGGTGAACTGTAGGAAGCTCTGAACACACACGTTCCCTTGGATTCGTACCTATGAATACTCCGT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.sequence[0] == 0x8187bf71cc0d7d7d, "{} failed", function_name!());
        assert!(obj.sequence[1] == 0x725cd3f7a2d7eb81, "{} failed", function_name!());
        assert!(obj.sequence[2] == 0xeca09de04446f57e, "{} failed", function_name!());
        assert!(obj.sequence[3] == 0x8f6c5ce0c75b0000, "{} failed", function_name!());
    }

// Detection of repetitive subsequence
    #[test]
    #[named]
    fn has_poly_base_test_8C(){
        let source: Vec<u8> = vec![b'C', b'C', b'C', b'C', b'C', b'C', b'C', b'C'];
        let obj = DnaSequence::new(&source);
        assert!(obj.has_poly_base(0, 4) == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 5) == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 6) == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 7) == (true, 0),  "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn has_poly_base_test_8G(){
        let source: Vec<u8> = vec![b'G', b'G', b'G', b'G', b'G', b'G', b'G', b'G'];
        let obj = DnaSequence::new(&source);
        assert!(obj.has_poly_base(0, 4) == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 5) == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 6) == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 7) == (true, 0),  "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn has_poly_base_test_8T(){
        let source: Vec<u8> = vec![b'T', b'T', b'T', b'T', b'T', b'T', b'T', b'T'];
        let obj = DnaSequence::new(&source);
        assert!(obj.has_poly_base(0, 4) == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 5) == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 6) == (true, 0),  "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 7) == (true, 0),  "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn has_poly_base_test_8N(){
        let source: Vec<u8> = vec![b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T'];
        let obj = DnaSequence::new(&source);
        assert!(obj.has_poly_base(0, 4) == (false, 0), "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 5) == (false, 0), "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 6) == (false, 0), "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 7) == (false, 0), "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn has_poly_base_test_120N(){
        let source: String = "GAACGACTGTTTTTACTATAAATCCTTCCTTCCTAGCCTATCATTTCTGGAGTCCTTGGTGAACTGTAGGAAGCTCTGAACACACACGTTCCCTTGGATTCGTACCTATGAATACTCCGT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 9)   == (false, 0), "{} failed", function_name!());
        assert!(obj.has_poly_base(0, 19)  == (true, 9),   "{} failed", function_name!());
        assert!(obj.has_poly_base(9, 20)  == (true, 0),   "{} failed", function_name!());
        assert!(obj.has_poly_base(28, 39) == (false, 0), "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn has_poly_base_test_27N_1(){
        let source: String = "ATTCATACTTAATACTGTATCAGTTGA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 27) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_27N_2(){
        let source: String = "TTCATACTTAATACTGTATCAGTTGAG".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 27) == (false, 0), "{} failed", function_name!());
    }
    #[test]
    #[named]
    fn has_poly_base_test_27N_3(){
        let source: String = "TCATACTTAATACTGTATCAGTTGAGT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 27) == (false, 0), "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn has_poly_base_test_27N_4(){
        let source: String = "TTCGAAAATCATCATCATCATCATCAC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 27) == (true, 4), "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn has_poly_base_test_27N_5(){
        let source: String = "TTCGAAATTCATCATCATCATCATCAC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 27) == (true, 4), "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn has_poly_base_test_27N_6(){
        let source: String = "TTCGTAATTCATCATCATCATCATCAC".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_poly_base(0, 27) == (false, 0), "{} failed", function_name!());
    }



    #[test]
    #[named]
    fn has_simple_repeat_27N_1(){
        let source: String = "TCATATGCAACAACAACTCATACTTAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_simple_repeat(0, 27) == (true, 7), "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn has_simple_repeat_27N_2(){
        let source: String = "TCATATGCTACAACAACTCATACTTAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_simple_repeat(0, 27) == (false, 0), "{} failed", function_name!());
    }


    #[test]
    #[named]
    fn has_2base_repeat_27N_1(){
        let source: String = "TCATATGCTACAACAACTCATACTTAA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_2base_repeat(0, 27) == (false, 0), "{} failed", function_name!());
    }


    #[test]
    #[named]
    fn has_2base_repeat_27N_2(){
        let source: String = "CGTACGCTTATATATATATATACCGCA".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_2base_repeat(0, 27) == (true, 9), "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn has_2base_repeat_27N_2x6(){
        let source: String = "TATATATATATAGCCCGCACGTACGCT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_2base_repeat(0, 27) == (true, 0), "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn has_2base_repeat_27N_2x5(){
        let source: String = "CGTATATATATAGCCCGCACGTACGCT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_2base_repeat(0, 27) == (false, 0), "{} failed", function_name!());
    }


    #[test]
    #[named]
    fn has_2base_repeat_27N_5(){
        let source: String = "ACGTTATACCTGTACCGCACGTACGCT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_2base_repeat(0, 27) == (false, 0), "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn has_2base_repeat_27N_6(){
        let source: String = "GCGGTATATACGTACCGCACGTACGCT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.has_2base_repeat(0, 27) == (false, 0), "{} failed", function_name!());
    }



//Decode test
    #[test]
    #[named]
    fn decode_test_8C(){
        let source: Vec<u8> = vec![b'C', b'C', b'C', b'C', b'C', b'C', b'C', b'C'];
        let obj = DnaSequence::new(&source);
        assert!(obj.decode(0, 4) == vec![67, 67, 67, 67            ], "{} failed", function_name!());
        assert!(obj.decode(0, 5) == vec![67, 67, 67, 67, 67        ], "{} failed", function_name!());
        assert!(obj.decode(0, 6) == vec![67, 67, 67, 67, 67, 67    ], "{} failed", function_name!());
        assert!(obj.decode(0, 7) == vec![67, 67, 67, 67, 67, 67, 67], "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn decode_test_120N(){
        let source: String = "GAACGACTGTTTTTACTATAAATCCTTCCTTCCTAGCCTATCATTTCTGGAGTCCTTGGTGAACTGTAGGAAGCTCTGAACACACACGTTCCCTTGGATTCGTACCTATGAATACTCCGT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        assert!(obj.decode(0, 120) == v, "{} failed", function_name!());
    }

//Subsequence test
    #[test]
    #[named]
    fn subsequence_as_u128_test_64G(){
        let source: String = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        let retval = obj.subsequence_as_u128(vec![[0 as usize, 64 as usize]]);
        let expectedvalue: u128 = 0xAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA;
        assert!(retval == expectedvalue, "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn subsequence_as_u128_test_10G(){
        let source: String = "GGGGGGGGGG".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        let retval = obj.subsequence_as_u128(vec![[0 as usize, 10 as usize]]);
        let expectedvalue: u128 = 0b10101010101010101010;
        assert!(retval == expectedvalue, "{} failed", function_name!());
    }

    #[test]
    #[named]
    fn subsequence_as_u128_test_12N(){
        let source: String = "ACGTAACCGGTT".to_string();
        let v: Vec<u8> = source.into_bytes();
        let obj = DnaSequence::new(&v);
        let retval = obj.subsequence_as_u128(vec![[0 as usize, 12 as usize]]);
        let expectedvalue: u128 = 0b000110110000010110101111;
        assert!(retval == expectedvalue, "{} failed", function_name!());
    }

}

