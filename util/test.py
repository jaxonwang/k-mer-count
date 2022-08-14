#! /usr/bin/env python3
from Bio import SeqIO
import argparse
import pprint
pp = pprint.PrettyPrinter(indent = 2)

def read_sequence(filename):
	retArray = []
	for seq_record in SeqIO.parse(filename, "fasta"):
		retArray.append(str(seq_record.seq))
	return retArray


def main():
	parser = argparse.ArgumentParser(description = "k-mer counter")
	parser.add_argument("fastafile", type=str, metavar="fastafile", help="input fasta file")
	args = parser.parse_args()
	filename = args.fastafile
	reads = read_sequence(filename)
	L_len = R_len = 27
	reads_LR = []
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
	reads_LR.sort()
	print("\n".join(reads_LR))




if __name__ == "__main__":
	main()



