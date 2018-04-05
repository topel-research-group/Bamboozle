#!/usr/bin/env python3

import argparse
import sys

parser = argparse.ArgumentParser(description='Check number of mutations appearing in exons')
parser.add_argument("-i", "--input", help="List of mutation coverages")
args = parser.parse_args()

coverage = []

with open(args.input,'r') as in_file:
	for line in in_file:
		coverage.append(line)

#print(coverage[0].strip("\n"))
#del coverage[0]
#print(coverage[0].strip("\n"))

while len(coverage) >= 2:
	before = coverage[0].split("\t")
	after = coverage[1].split("\t")
	if len(coverage) >= 3:
		next = coverage[2].split("\t")
	if int(after[1]) - int(before[1]) == 1:
		if (int(before[2]) == 0 and int(after[2]) == 0) or (int(before[2]) == 0):
			print(after[0] + "\t" + after[1] + "\tN/A")
		else:
			print(after[0] + "\t" + after[1] + "\t" + str(round(((100.0/int(before[2]))*int(after[2])), 2)))
		del coverage[0]
		if int(next[1]) - int(after[1]) != 1:
			del coverage[0]
	else:
		print("Whoops, something went wrong!")
		sys.exit()
