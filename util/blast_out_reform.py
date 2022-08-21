#! /usr/bin/env python3
import argparse
import pprint
import re
from abc import ABC
from dataclasses import dataclass

pp = pprint.PrettyPrinter(indent = 2)

@dataclass
class Record(ABC):
	seqid    : int
	role     : str
	qseqid   : str
	sseqid   : str
	pident   : float
	length   : int
	mismatch : int
	gapopen  : int
	qstart   : int
	qend     : int
	sstart   : int
	send     : int
	evalue   : str
	bitscore : str

def main():
	parser = argparse.ArgumentParser(description = "blastn-short result reformer")
	parser.add_argument("blast_result", type=str, metavar="Blast result", help="Blast results file (outfmt 6)")
	parser.add_argument("primer_candidates", type=str, metavar="Primer candidates", help="Primer candidate name list file. One candidate in u128 in a line.")
	args = parser.parse_args()
	filename = args.blast_result
	candidates_pack_raw = []
	records = []
	primer_pairs = set()
	with open(filename) as f:
		for line in f:
			elm = line.strip().split("\t")
			seqid = int(elm[0].split("-")[0], 16)
			role = elm[0].split("-")[1]
			qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore = elm
			pident   = float(pident)
			length   = int(length)
			mismatch = int(mismatch)
			gapopen  = int(gapopen)
			qstart   = int(qstart)
			qend     = int(qend)
			sstart   = int(sstart)
			send     = int(send)
			primer_pairs.add(seqid)
			records.append(Record(seqid, role, qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore))

	discarded_primer_pairs = set()
	for each_primer_pair in primer_pairs:
		subject_to_be_checked = list(filter(lambda x: x.seqid == each_primer_pair, records))
		roles_in_each_primer_pair = set()
		for i in subject_to_be_checked:
			roles_in_each_primer_pair.add(i.role)
		if len(roles_in_each_primer_pair) == 3:
			sseqids = set([x.sseqid for x in subject_to_be_checked])
			for each_sseqid in sseqids:
				role_L_coverage = [(x.sstart, x.send) for x in subject_to_be_checked if x.role == "L" and x.sseqid == each_sseqid]
				role_M_coverage = [(x.sstart, x.send) for x in subject_to_be_checked if x.role == "M" and x.sseqid == each_sseqid]
				role_R_coverage = [(x.sstart, x.send) for x in subject_to_be_checked if x.role == "R" and x.sseqid == each_sseqid]
				if role_L_coverage != [] and role_M_coverage != [] and role_R_coverage != []:
					# test all combination
					for l in role_L_coverage:
						for m in role_M_coverage:
							for r in role_R_coverage:
								if 0 < m[0] - l[1] < 76 and 0 < r[0] - m[1] < 76:
									discarded_primer_pairs.add(each_primer_pair)

	primer_candidates_set = set()
	with open(args.primer_candidates) as f:
		for line in f:
			primer_candidates_set.add(int(line, 16))

	for i in primer_candidates_set - discarded_primer_pairs:
		print(f'{hex(i).replace("0x", "")}')
	print(f"{primer_candidates_set} - {discarded_primer_pairs} = {primer_candidates_set - discarded_primer_pairs} ({discarded_primer_pairs / (primer_candidates_set - discarded_primer_pairs)})", file = sys.stderr)



if __name__ == "__main__":
	main()



