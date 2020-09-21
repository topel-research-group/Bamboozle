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
parser.add_argument("-p", "--population", action="append" \
	help="Text file of newline-separated individuals in a population.\
	This argument can be used multiple times to define different populations")

args = parser.parse_args()

#makes sure that SnpEff annotations don't get in the way of IDing a gene
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
	#find all the genes in a GFF file, creates a list, orders and finds unique genes
	with open(gff[0], "r") as gff_file:
		for line in gff_file:
			if line.startswith(("#", "A", "C", "T", "G", ">")):
				continue
			ann = line.split("\t")[8]
			id = ann.split(";")[0]
			genes.append(id.replace('ID=',''))

	genes_std_uniq = sorted(set(genes))

	#read in vcf (gzipped files are ok, pysam likes them), create pandas df
	#based on genes in gff, samples in VCF
	vcf_in = VariantFile(",".join(input_vcf))
	data = pd.DataFrame(0, \
		columns = vcf_in.header.samples, \
		index = list(genes_std_uniq))

	#for each line in the vcf, check that it's not in a plastid or mit contig
	#check that it's been annotated by SnpEff, check if for that sample there's a genotype
	#meaning sample in the VCF has been annotated, take gene annotation, strip it down
	#with parse_gene(), add, length of line.info['ANN'] tuple, since it annotates number of 
	#modifications for that gene at that position, add to new val, write it to pandas df created
	#return df so compare() can use it
	gene_lof = None
	for line in vcf_in:
		if line.chrom not in ["Sm_plastid", "Sm_mitochondrion"]:
			if 'ANN' in line.info:
				for sample in list(line.header.samples):
					if (line.samples[sample]['GT']).count(None):
						gene = line.info['ANN'][0].split("|")[3]
						gene_lof = parse_gene(gene)
	
						data.at[gene_lof, sample] =+ 0
					else:
						gene = line.info['ANN'][0].split("|")[3]
						gene_lof = parse_gene(gene)

						new_val = int(data.at[gene_lof, sample]) + len(line.info['ANN'])
						data.at[gene_lof, sample] = new_val
	data.to_csv(out_prefix[0]+'.csv', sep='\t')
	return data

def compare(data, out_prefix, pops):
	#should pops be like pop1: x, y, z? what format?
	#one list per pop?
	#take the df, find averages, stdvs for each gene across samples in each pop
	#test differences, significance in LOF mutations between pops?
	for pop in pops:
		with open(pop) as infile:
			pop_s = splitext(pop)[0]
			pop_s+"_list" = []
			for ind in infile:
				pop+"_list".append(ind)
	
				
	
def main():
	gen_matrix(args.input_vcf, args.gff, args.out_prefix)
	compare(data, args.out_prefix, args.population)

if __name__ == "__main__":
    main()
