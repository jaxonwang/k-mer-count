#! /usr/bin/env python3
import argparse
import pprint
import sys
pp = pprint.PrettyPrinter(indent = 2)

def main():
	parser = argparse.ArgumentParser(description = "return intersection of two file")
	parser.add_argument("-file1", "--file1", required=True, type=str, metavar="file 1", help="file1")
	parser.add_argument("-file2", "--file2", required=True, type=str, metavar="file 2", help="file2")
	args = parser.parse_args()
	file1_set = set()
	file2_set = set()
	with open(args.file1, 'r') as f1:
		for line in f1:
			file1_set.add(line.strip())

	with open(args.file2, 'r') as f2:
		for line in f2:
			file2_set.add(line.strip())
	for i in file1_set & file2_set:
		print(i)

	

if __name__ == "__main__":
	main()