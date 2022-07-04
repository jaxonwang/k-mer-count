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
            buf <<= 2;
            cnt += 1;
            if cnt == 31{//cnt == 31のときはbufが全部埋まってる。
                retval.push(buf);
                buf = 0;
                cnt = 0;
            }
        }
        if cnt != 31{
            buf <<= 2 * (31 - cnt);//塩基を表すbitを上位に寄せる。
            retval.push(buf);//32塩基の節目で切れなかった時に備えてpushする。
        }
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
        assert!(end - start < 32, "DnaSequence::has_poly_base assertion failed: length of the evaluation subject must be shorter than 32");
        assert!(end < self.length, "DnaSequence::has_poly_base assertion failed: end coordinate must be smaller than length of the sequence. wnd: {}, self.lngth: {}", end, self.length);
        let mut original:  u64 = 0;
        let mut zero_ichi: u64 = 0;
        for i in start..end{
            original += (self.sequence[i / 32] >> (62 - 2 * (i % 32))) & 3;
            original <<= 2;
            zero_ichi += 1;
            zero_ichi << 2;
        }//ここまでで、originalに右詰で対象の領域がコピーされる。
        let val4 = original ^ (original >> 2);
        let val5 = val4 | (val4 >> 1);
        let val7 = !val5 & zero_ichi;
        return (val7 & (val7 >> 2)) & (val7 >> 4) != 0;
    }
}
