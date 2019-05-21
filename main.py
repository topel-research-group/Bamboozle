#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess
import fnmatch
from functools import reduce
from functools import wraps
from time import time
import datetime

#######################################################################

parser = argparse.ArgumentParser(prog="Bamboozle")
parser.add_argument("-f", "--ref", \
                        help="Reference")
parser.add_argument("-F", "--forward", \
                        nargs='*', \
                        help="Forward reads")
parser.add_argument("-R", "--reverse", \
                        nargs='*', \
                        help="Reverse reads")
parser.add_argument("-b", "--bamfile", \
                        help="BAM infile")
parser.add_argument("--sortbam", \
                        help="Sorted BAM infile")
parser.add_argument("--gff", \
                        help="gff infile")
parser.add_argument("--contigsizes", \
                        help="Contig sizes for gff parser")
parser.add_argument("--feature", \
                        help="Feature for gff parser")
parser.add_argument("-t", "--threads", \
                        default=1, \
                        help="Threads")
parser.add_argument("-e", "--snpeff", \
                        nargs='*', \
                        help="Input options for snpeff, without the '-' before")
parser.add_argument("-s", "--snpsift", \
                        action="store_true", \
                        help="Run snpSift")
parser.add_argument("-r", "--clean", \
                        action="store_true", \
                        help="Removes the SAM and BAM files")
parser.add_argument("-p", "--done", \
                        action="store_true", \
                        help="Add an empty file to mark the directory as done")

#######################################################################

group = parser.add_argument_group('Bamparser')
group.add_argument('--coverage', \
			action="store_true", \
			help='Print a statistic for what percentage of bases in an assembly have >=Nx coverage')
group.add_argument('--consensus', \
			action="store_true", \
			help='Extract the consensus sequence of aligned reads from a specific region of the reference sequence (WIP)')
group.add_argument('--zero', \
			action="store_true", \
			help='Find areas of zero coverage and print the reference sequence, along with a GC percentage')
group.add_argument('--deletion1', \
			action="store_true", \
			help='Find deletions')
group.add_argument('--deletion2', \
			action="store_true", \
			help='Find deletion events')
group.add_argument('--deletion3', \
			action="store_true", \
			help='Find frameshift deletion events')
group.add_argument('--deletionx', \
			action="store_true", \
			help='Find deletions occurring within exons')
group.add_argument('--homohetero', \
			action="store_true", \
			help='Attempt to determine whether a deletion is homozygous or heterozygous')
group.add_argument('--median', \
			action="store_true", \
			help='Find regions differing from contig median by +/- 50%%, or just contig medians')
group.add_argument('--long_coverage', \
			action="store_true", \
			help='Find the longest region between given coverage limits for a given contig')
group.add_argument('--complex', \
			action="store_true", \
			help='Print full bed output for median')
group.add_argument('--simple', \
			action="store_true", \
			help='Print median coverage only for median')
group.add_argument('-c', '--contig', \
			help='Gives per-contig coverage stats')
group.add_argument('-d', '--threshold', \
			type=int, \
			nargs='?', \
			const='1', \
			default='20', \
			help='Threshold for calculating coverage percentage; default 20')
group.add_argument("-a", "--range", \
			help="somethingsomsing")
group.add_argument("-m", "--mutations", \
			help="List of mutation events; currently requires output from bamboozle deletion function")
group.add_argument("-x", "--exons", \
			help="Bed file containing exon coordinates (0-based); -m also required")
group.add_argument("-l", "--limits", \
			type=int, \
			nargs=2, \
			help="Specify lower and upper limits for long_coverage function; two arguments required")
group.add_argument("-v", "--verbose", \
			action="store_true", \
			help="Be more verbose")
group.add_argument('--dev', \
			help=argparse.SUPPRESS, action="store_true")

args = parser.parse_args()

if args.feature and args.gff is None:
        parser.error("--feature requires --gff")
elif args.gff and args.feature is None:
        parser.error("--feature requires --gff")

if args.forward and args.reverse and args.ref is None:
	parser.error("--ref [Reference is required]")

#######################################################################

# For determining how long functions take to run
if args.dev == True:
	start_time = time()

#######################################################################

sortbam = args.sortbam

#######################################################################

# Determine whether the function derives from `bamparser.py`
# Avoids the need for a `--bamparse` flag
BamparseList = ["--coverage","--consensus","--zero","--deletion1","--deletion2","--deletion3",\
		"--deletionx","--homohetero","--median","--long_coverage"]

for item1 in BamparseList:
	for item2 in sys.argv:
		if item1 == item2:
			bamparse = True

#######################################################################

# Run pipeline from beginning and if --bamparser run bamparser.py
if args.ref and args.forward and args.reverse:
	if bamparse:
#		from pipeline import main
#		main()

		from pipeline import bowtie2,samtools_view,bam_input,samtools_index
#		from pipeline import bowtie2,samtools_view,samtools_sort,samtools_index
		bowtie2(args)
		samtools_view()
#		samtools_sort()
		bam_input(args)
		samtools_index()

		sortbam = 'Bowtie2/*_sorted.bam' 
		from bamparser import main
		main(sortbam)
	else:
		from pipeline import main
		main()
	
# Run pipeline from bam_input (skips aligning steps) and if --bamparse run bamparser.py
if args.bamfile:
	if bamparse:
		from pipeline import bam_input,samtools_index
		bam_input(args)
		samtools_index()

		sortbam = 'Bowtie2/*_sorted.bam' 
		from bamparser import main 
		main(sortbam)
	else:
		from pipeline import input_files,snpEff_test
		snpEff_test()
		input_files()	

# If input is a sorted bam file and reference run pipeline, if not and --bamparse run bamparser.py
if args.sortbam and bamparse:
	if args.ref:
		from pipeline import input_files,snpEff_test
		snpEff_test()
		input_files()

		from bamparser import main
		main(sortbam)
	else:
		from bamparser import main
		main(sortbam)

#######################################################################

if args.dev == True:
	print("Time taken =",(time() - start_time),"seconds.")
