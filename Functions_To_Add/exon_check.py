#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description='Check number of mutations appearing in exons')
parser.add_argument("-b", "--bed", help="Bed file")
parser.add_argument("-i", "--input", help="List of mutations")
args = parser.parse_args()

frameshifts = 0
exons = []
mutations = []

with open(args.bed,'r') as bed_file:
	for line in bed_file:
		exons.append(line)

with open(args.input,'r') as in_file:
	for line in in_file:
		mutations.append(line)

for mutation in mutations:
	mutation = mutation.split("\t")
	for exon in exons:
		exon = exon.split("\t")
		if (mutation[1] == exon[0]) and (int(exon[1]) <= int(mutation[2]) <= int(exon[2])):
			frameshifts += 1
			print(mutation[3].strip("\n") + "bp mutation at " + mutation[1] + " " + mutation[2] + \
			      " hits exon " + exon[3] + ".")
			break

print("Total number of frameshifts: " + str(frameshifts))
