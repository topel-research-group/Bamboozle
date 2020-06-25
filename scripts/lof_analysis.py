#!/usr/bin/env python3

import argparse
from pysam import VariantFile

#read arguments
parser = argparse.ArgumentParser(description='Process, summarize and plot VCF files and corresponding metadata')
parser.add_argument("-v", "--input-vcf", nargs= '*', type=str, \
			help="Input full path to VCF file or files. A single file, a list of comma-separated files, or a file with paragraphed file paths")
parser.add_argument("-m", "--metadata", nargs= 1, type=str, \
			help="Input metadata in a TSV format. First column should be named ID and consist of sample identifiers")
#saving these for later
parser.add_argument("--no-circos", \
			help="Circos plots will not be generated for the sample(s). Returns an HTML report with SV plots and a summarized table")
parser.add_argument("--no-plots", \
			help="SV plots will not be generated for the sample(s). Returns an HTML report with Circos plots and a summarized table")
parser.add_argument("--just-the-table", \
			help="Will only produce a summarized table from input VCF(s) and metadata")

args = parser.parse_args()

#detect if only one file, several comma-delimited ones or a .txt file with names of files
#return a state to inform other functions on how to work
def check_input(input-vcf):
	#if it is indeed a string verify if a path to a file or to a .txt with data
	if isinstance(input-vcf, str):
		if input-vcf.endswith(".vcf"):
			state="single vcf"
			return state
		if input-vcf.endswith(".txt"):
			state="vcfs in file"	
			return state
	if isinstance(input-vcf, list):
		state="vcfs in list"
		return state

#summarize vcfs with metadata
def summarize(input-vcf, metadata):
	#taking care of vcf first according to the nature of the input
	if state == "single vcf":
		#read in vcf input, take its name
		vcf_in = VariantFile(input-vcf)
		
	elif state == "vcfs in file":
		#
	elif state == "vcfs in list":
		#
#read input metadata table
#
#produce full output
#main
def main():
	check_input(args.input-vcf)
	#follow arguments if any has been given
	if args.command == "no-circos":
		summarize(args.input-vcf, args.metadata, state)
		plots()
	elif args.command == "no-plots":
		summarize(args.input-vcf, args.metadata, state)
		circos()
	elif args.command == "just-the-table":
		summarize(args.input-vcf, args.metadata, state)
	#otherwise run everything
	else:
		summarize(args.input-vcf, args.metadata, state)
		plots()
		circos()
