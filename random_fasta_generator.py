#! /usr/bin/env python3
import random

for i in range(1, 201):
	num = str(i).zfill(3)
	print(f">dummy_sequence_{num} {i}th record")
	for j in range(5):
		random_str = "".join([random.choice(["A", "C", "G", "T"]) for x in range(80)])
		print(random_str)