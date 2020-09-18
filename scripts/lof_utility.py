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
	help="Input full path to multi-VCF file. ")
parser.add_argument("-g", "--gff", \
	nargs= 1, type=str, \
	help="Reference GFF.")
parser.add_argument("-o", "--out_prefix", \
	nargs = 1, type=str, \
	help="For multiple samples, add an output prefix")

args = parser.parse_args()

def parse_gene(gene):
	if '-' in gene:
		for i in gene.split("-"):
			if i.startswith("Sm"):
				gene_lof = i
				return gene_lof
	elif '&' in gene:
		for i in gene.split("&"):
			if i.startswith("Sm"):
				gene_lof = i
				return gene_lof
	else:
		gene_lof = gene
		return gene_lof

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

	gene_lof = None
	for line in vcf_in:
		if 'ANN' in line.info:
			for sample in list(line.header.samples):
				if (line.samples[sample]['GT']).count(None):
					gene = line.info['ANN'][0].split("|")[3]
					gene_lof = parse_gene(gene)
					data.loc[gene_lof, sample] == 0
				else:
					gene = line.info['ANN'][0].split("|")[3]
					gene_lof = parse_gene(gene)
					data.loc[gene_lof, sample] =+ len(line.info['ANN'])


	#can't add to cell value (=+) in pd df...
	#something that could solve this would be to generate a per-gene list of LOF counts per sample
	#then adding that to each row by gene

#	pd.set_option('display.max_columns', None)
#	print(data.head())
#	data.to_csv(out_prefix[0]+'.csv', sep='\t')
	return data

#def compare(data, out_prefix, pops):
	#should pops be like pop1: x, y, z? what format?
	#eliminate genes not found for example
	#take the df, find averages, stdvs for each pop
	#difference, significance in LOF mutations between pops?
	
def main():
	gen_matrix(args.input_vcf, args.gff, args.out_prefix)

if __name__ == "__main__":
    main()
