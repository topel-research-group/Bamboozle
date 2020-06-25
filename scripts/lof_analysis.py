#!/usr/bin/env python3

import argparse
from pysam import VariantFile

#read arguments
parser = argparse.ArgumentParser(description='Process, summarize and plot VCF files and corresponding metadata')
parser.add_argument("-v", "--input-vcf", nargs= '*', \
			help="Input full path to VCF file or files. A single file, a list of comma-separated files, or a file with paragraphed file paths")
parser.add_argument("-m", "--metadata", nargs= 1, \
			help="Input metadata in a TSV format. First column should be named ID and consist of sample identifiers")
#saving these for later
#parser.add_argument("--no-circos", \
#			help="Circos plots will not be generated for the sample(s). Returns an HTML report with SV plots and a summarized table")
#parser.add_argument("--no-plots", \
#			help="SV plots will not be generated for the sample(s). Returns an HTML report with Circos plots and a summarized table")
#parser.add_argument("--just-the-table", \
#			help="Will only produce a summarized table from input VCF(s) and metadata")

args = parser.parse_args()

#read input vcf(s)
#detect if only one file, several comma-delimited ones or a .txt file with names of files
bcf_in = VariantFile("test.bcf")  # auto-detect input format
#read input metadata table
#
#produce full output
