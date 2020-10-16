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
import scipy.stats as stats

##INPUT:

##OUTPUT:


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
parser.add_argument("-p", "--population", action="append", \
	help="Text file of newline-separated individuals in a population. This argument can be used multiple times to define different populations")

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
			if line.startswith("##FASTA"):
				break
			elif not line.startswith("#"):
				ann = line.split("\t")[8]
				id = ann.split(";")
				if len(id) > 2:
					name = id[2]
					genes.append(name.replace('gene_name=','').rstrip("\n"))

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
		if 'ANN' in line.info:
			for sample in list(line.header.samples):
				if (line.samples[sample]['GT']).count(None):
					continue
#					gene = line.info['ANN'][0].split("|")[3]
#					gene_lof = parse_gene(gene)

#					data.at[gene_lof, sample] =+ 0
				else:
					gene = line.info['ANN'][0].split("|")[3]
					gene_lof = parse_gene(gene)

		#			print(gene_lof)
		#			print(sample)

# 	https://stackoverflow.com/questions/44849868/python-df-loc-not-working-for-variables
		#			print("original val", data.loc[gene_lof, sample])
					new_val = int(data.loc[gene_lof, sample]) + 1
		#			new_val = int(data.loc[gene_lof, sample]) + len(line.info['ANN'])
		#			print("new_val\n", new_val)
					data.loc[gene_lof, sample] = new_val
		#			print("new_val in df\n",  data.loc[gene_lof, sample])
	data.to_csv(out_prefix[0]+'.csv', sep='\t')
	return data, genes

def compare(out_prefix, pops, data, genes):

	sig = pd.DataFrame(index=data.index)
	list_of_pops = []

	for pop in pops:
		with open(pop) as infile:
			pop_s = pop.split(".")[0]
			pop_s = []
			for ind in infile:
				pop_s.append(ind.replace("\n",""))
	
			list_of_pops.append(pop_s)

	for group in list_of_pops:
		for gene in genes:
			print(data[group].T[gene])


	#find a way to call all inds in a pop so stats can be calculated
#			mean_col = pop.split(".")[0]+"_mean"
#			std_col = pop.split(".")[0]+"_std"
#			stat_col = pop.split(".")[0]+"_stat"
#			p_col = pop.split(".")[0]+"_sig"

#			sig[mean_col] = data[pop_s].mean(axis=1)
#			sig[std_col] = data[pop_s].std(axis=1)
			
			pop_id = pop.split(".")[0]
#			print(type(sig))

def main():
	gen_matrix(args.input_vcf, args.gff, args.out_prefix)
#	data, genes = gen_matrix(args.input_vcf, args.gff, args.out_prefix)
#	data, list_of_pops = compare(args.out_prefix, args.population, data, genes)
#	result = pd.DataFrame(compare(args.out_prefix, args.population, data, list_of_pops))
#	print(result.head)
			

if __name__ == "__main__":
    main()
