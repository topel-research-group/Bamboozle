#!/usr/bin/env python3

import argparse
import os
import pandas as pd
from pysam import VariantFile

#read arguments
parser = argparse.ArgumentParser(description='Process, summarize and plot VCF files and corresponding metadata')

#args in
parser.add_argument("-v", "--input_vcf", \
	nargs= '*', type=str, required=True \
	help="Input full path to VCF file or files. \
	A single file, a list of comma-separated files, or a file with paragraphed file paths")
parser.add_argument("-m", "--metadata", \
	nargs= 1, type=str, \
	help="Input metadata in a TSV format.")
parser.add_argument("-o", "--out_prefix", \
	nargs = 1, type=str, \
	help="For multiple samples, add an output prefix")
#saving these for later
parser.add_argument("--no_plots", \
	help="SV plots will not be generated for the sample(s). \
		Returns an HTML report with Circos plots and a summarized table")
parser.add_argument("--just_the_table", \
	help="Will only produce a summarized table from input VCF(s) and metadata")

args = parser.parse_args()

#detect if only one file, several comma-delimited ones or a .txt file with names of files
#return a state to inform other functions on how to work
def check_input(input_vcf):
	#if it is indeed a string verify if a path to a file or to a .txt with data
	if isinstance(input_vcf, list):
		if len(input_vcf) == 1:
			if input_vcf[0].endswith(".vcf"):
				return "single vcf"
			if input_vcf[0].endswith(".txt"):
				return "vcfs in file" 
		if len(input_vcf) > 1:
			return "vcfs in list"

#summarize vcfs per sample (with metadata?)	-	sorted
#what should this do?...
#
#read input metadata table
#       with open(metadata) as infile:
#               md = csv.reader(infile, delimiter="\t")
#               for row in md:
                        ##do XXX

#output a bed-formatted output with metadata in it to input into the html plotly
#plots there:
#	horizontal 'barplot' with SVs per chrm with annotations highlighted
#		1. needs chromosome layout per sample, that is
#			proportional to contig/chrom size
#			"mappable"? to locs in BED input
#		2. needs to plot bars? at right height too
#		3. maybe dendrogram for n. of samples plotted
#	point plot? comparing populations by hierarchical clustering?

def summarize(input_vcf, state, out_prefix):
	#taking care of vcf first according to the nature of the input
	if state == "single vcf":
		#read in vcf input, take its name
		vcf_in = VariantFile(",".join(input_vcf))
		data = pd.DataFrame(0, \
			columns = ['DEL','INS','DUP','INV','CTX','UNC'], \
			index = list(vcf_in.header.contigs))
		#for all the chromosomes found in the vcf keep as row names
		for line in vcf_in:
			if 'SIMPLE_TYPE' in line.info:
				data.loc[line.chrom, line.info['SIMPLE_TYPE']] += 1
			else:
				data.loc[line.chrom, 'UNC'] += 1
		data_out = ",".join(input_vcf)[:-4] + ".tsv"
		data.to_csv(data_out, sep='\t')
	#if more than one file as comma-sep inputs
	if state == "vcfs in list":
		#create folder if it doesn't exist
		if not os.path.exists(str(out_prefix).strip('[]')[1:-1]):
			os.makedirs(str(out_prefix).strip('[]')[1:-1])
		for vcf in input_vcf:
			vcf_in = VariantFile(vcf.strip("`b,"))
			data_multi = pd.DataFrame(0, \
				columns = ['DEL','INS','DUP','INV','CTX','UNC'], \
				index = list(vcf_in.header.contigs))
			#for all the chromosomes found in the vcf keep as row names
			for line in vcf_in:
				if 'SIMPLE_TYPE' in line.info:
					data_multi.loc[line.chrom, line.info['SIMPLE_TYPE']] += 1
				else:
					data_multi.loc[line.chrom, 'UNC'] += 1
			#name output per sample
			data_out = vcf[:-4] + ".tsv"
			data_multi.to_csv(str(out_prefix).strip('[]')[1:-1] + "/" + data_out, sep='\t')

	#if more than one file as .txt with \n-sep inputs
	if state == "vcfs in file":
		if not os.path.exists(str(out_prefix).strip('[]')[1:-1]):
			os.makedirs(str(out_prefix).strip('[]')[1:-1])
		open_txt = open(input_vcf[0])
		vcf_open = open_txt.readlines()
		for vcf in vcf_open:
			vcf_in = VariantFile(vcf.strip("`b,").rstrip("\n"))
			data_multi = pd.DataFrame(0, \
				columns = ['DEL','INS','DUP','INV','CTX','UNC'], \
				index = list(vcf_in.header.contigs))
			#for all the chromosomes found in the vcf keep as row names
			for line in vcf_in:
				if 'SIMPLE_TYPE' in line.info:
					data_multi.loc[line.chrom, line.info['SIMPLE_TYPE']] += 1
				else:
					data_multi.loc[line.chrom, 'UNC'] += 1
		#there goes the output
			data_out = vcf.strip("`b,").rstrip("\n")[:-4] + ".tsv"
			data_multi.to_csv(str(out_prefix).strip('[]')[1:-1] + "/" + os.path.basename(data_out), sep='\t')

#python dashboard
#interactive plots?
#maybe like a heatmap?

def main():
	#start things out
	state = check_input(args.input_vcf)
	#follow arguments if any has been given
	if args.no_plots:
		summarize(args.input_vcf, state)
		#circos()
	if args.just_the_table:
		summarize(args.input_vcf, state)
	#otherwise run everything
	else:
		summarize(args.input_vcf, state, args.out_prefix)
#		plots()
#		circos()

if __name__ == "__main__":
    main()
