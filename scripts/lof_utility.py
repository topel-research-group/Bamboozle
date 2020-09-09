#!/usr/bin/env python3
# coding=utf-8

#André Soares - 03/09/2020

#       This script will generate a samples/genes LOF matrix for a multi-sample VCF
#
#       Copyright (C) 2020 André Soares, andre.soares@bioenv.gu.se
#
#       This program is free software: you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation, either version 3 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program.  If not, see <https://www.gnu.org/licenses/>.



#1. INPUT: multi-sample filtered VCF
#2. Extract information from multi-VCF:
#	probably pandas
#	using pysam's VariantFile
#		for each sample in the VCF
#			for each gene in each sample
#				populate the sample/gene index with SO Annotation
#
# what to do for >1 events per gene? separate by , in the same cell?
# would be ready for UpsetR like viz for example



import argparse
import os
import pandas as pd
import subprocess
from pysam import VariantFile

#read arguments
parser = argparse.ArgumentParser(description='Provide an occurrence matrix of LOF-affected genes given a multi-sample VCF')

#args in
parser.add_argument("-v", "--input_vcf", \
	nargs= '*', type=str, required=True, \
	help="Input full path to VCF file or files.")
parser.add_argument("-g", "--gff", \
	nargs= 1, type=str, \
	help="Reference GFF.")
parser.add_argument("-o", "--out_prefix", \
	nargs = 1, type=str, \
	help="For multiple samples, add an output prefix")

args = parser.parse_args()

def gen_matrix(input_vcf, gff, out_prefix):
	genes = []

	with open(gff[0], "r") as gff_file:
		for line in gff_file:
			if line.startswith(("#", "A", "C", "T", "G", ">")):
				continue
			ann = line.split("\t")[8]
			id = ann.split(";")[0]
			genes.append(id.replace('ID=',''))

	genes_std_uniq = sorted(set(genes))

	#this is ok with gzipped files	
	vcf_in = VariantFile(",".join(input_vcf))
	data = pd.DataFrame(0, \
		columns = vcf_in.header.samples, \
		index = list(genes_std_uniq))

	for line in vcf_in:
		if 'ANN' in line.info:
			data.loc[gene, line.header.samples] =+ len(line.info['ANN'])

	data.to_csv('empty_table.csv', sep='\t')

def main():
	gen_matrix(args.input_vcf, args.gff, args.out_prefix)

if __name__ == "__main__":
    main()
