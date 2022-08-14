#! /usr/bin/env python3
import argparse
import pprint
import re
pp = pprint.PrettyPrinter(indent = 2)


def main():
	parser = argparse.ArgumentParser(description = "primer3 result parser")
	parser.add_argument("input", type=str, metavar="file", help="input file")
	args = parser.parse_args()
	filename = args.input
	candidates_pack_raw = []
	with open(filename) as f:
		current_pack = []
		for i in f:
			each_line = i.strip()
			if each_line != "=":
				current_pack.append(each_line)
			else:
				candidates_pack_raw.append(current_pack)
				current_pack = []

	p = re.compile(r'_\d')

	candidates_pack = {}
	for each_pack in candidates_pack_raw:
		candidate_dict = {}
		primer_pair_dict = {}
		for each_tag in each_pack:
			tag, value = each_tag.split("=")
			if re.search(p, tag) == None:#_1_とかがないとき
				candidate_dict[tag] = value
			else:
				group = re.search(p, tag).group(0).replace("_", "")
				if group in primer_pair_dict.keys():
					primer_pair_dict[group][tag] = value
				else:
					primer_pair_dict[group] = {}
					primer_pair_dict[group][tag] = value
			candidate_dict["Primers"] = primer_pair_dict
		candidates_pack[candidate_dict["SEQUENCE_ID"]] = candidate_dict
	pp.pprint(candidates_pack)


if __name__ == "__main__":
	main()

