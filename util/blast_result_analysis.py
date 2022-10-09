#! /usr/bin/env python3

import argparse
import re
import pprint
import json
import sys
import csv
from Bio import SeqIO
pp = pprint.PrettyPrinter(indent = 2)
from primer3_result_parser import primer3_result_parser, extract_primer_pairs


class BlastResult():
	def __init__(self, qseqid, sseqid, sacc, slen, qstart, qend, sstart, send, qseq, sseq, evalue, length, staxid, staxids, ssciname, scomname):
		self.qseqid   = qseqid
		self.sseqid   = sseqid
		self.sacc     = sacc
		self.slen     = slen
		self.qstart   = qstart
		self.qend     = qend
		self.sstart   = sstart
		self.send     = send
		self.qseq     = qseq
		self.sseq     = sseq
		self.evalue   = evalue
		self.length   = length
		self.staxid   = staxid
		self.staxids  = staxids
		self.ssciname = ssciname
		self.scomname = scomname


def grep_input_record_name(primer_pairs_dict):
	retset = set()
	for k, v in primer_pairs_dict.items():
		for idx, each_pair in enumerate(v["Primer3_output"]):
			retset.add(f">{k}_{idx}_L")
			retset.add(f">{k}_{idx}_M")
			retset.add(f">{k}_{idx}_R")
	return retset


def main():
	parser = argparse.ArgumentParser(description = "compare input fasta file and blast result")
	parser.add_argument("fasta", metavar = "fasta", type = str, help = "fasta file name")
	parser.add_argument("blast", metavar = "blast", type = str, help = "blast output file name. outfmt must be '6 qseqid sseqid sacc slen qstart qend sstart send qseq sseq evalue length staxid staxids ssciname scomname'")
	parser.add_argument("-o",    metavar = "output_file",    type = str, default = "sys.stdout", help = "output file name (default = sys.stdout)")
	args = parser.parse_args()
	filename = args.fasta
	fasta_ids = set()
	with open(filename) as handle:
		for record in SeqIO.parse(handle, "fasta"):
			fasta_ids.add(record.id)

	blast_results = []
	with open(args.blast) as f:
		reader = csv.reader(f, delimiter='\t')
		for each_record in reader:
			blast_results.append(BlastResult(*each_record))

	blast_trapped_seq_ids = set()
	for each_hit in blast_results:
		blast_trapped_seq_ids.add(each_hit.qseqid)

	print(f"total count of input sequence                  : {len(fasta_ids)}")
	print(f"total count of blast hits                      : {len(blast_results)}")
	print(f"cardinarity of input which was trapped by blast: {len(blast_trapped_seq_ids)}")
	print(f"survivor                                       : {len(fasta_ids) - len(blast_trapped_seq_ids)}")
	pp.pprint(fasta_ids - blast_trapped_seq_ids)


if __name__ == '__main__':
	main()