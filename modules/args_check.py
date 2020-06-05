#!/usr/bin/env python3

import argparse
import sys

# Can information be obtained from parser._subparsers about the parser names?
# subparsers.choices is a dictionary whose keys are the subprocesses
# for i in subparsers.choices.keys():
#	print(i)
# Incorporate this into metavars?

def metavar_mod(metavar):
	x = list(metavar.values())
	y = list(metavar.keys())
	for i in range(0,len(x)):
		yield(str(y[i]) + "\t" + str(x[i].description))


parser = argparse.ArgumentParser(usage="bamboozle.py <command> <args>")

subparsers = parser.add_subparsers(title="Commands", \
					metavar=metavar_mod(subparsers.choices))

#					metavar="""
#	pipeline	[Vilma's pipeline]
#	coverage	Print a statistic for what percentage of bases in an assembly have >=Nx coverage
#	consensus	Extract the consensus sequence of aligned reads from a region of the reference (WIP)
#	zero		Find areas of zero coverage and print the reference sequence, along with a GC percentage
#	deletion1	Find deletions
#	deletion2	Find deletion events
#	deletion3	Find frameshift deletion events
#	deletionx	Find deletions occurring within exons
#	homohetero	Attempt to determine whether a deletion is homozygous or heterozygous
#	median		Find regions differing from contig median by +/- 50%%, or just contig medians
#	long_coverage	Find the longest region between given coverage limits for a given contig
#	barcode		Search the input (sorted) BAM files for suitable barcode regions
#	
#						""")


# Arguments for all commands

all_commands = argparse.ArgumentParser()

all_commands.add_argument("-F", "--forward", \
				nargs='*', \
				help="Forward reads")



pipeline = subparsers.add_parser("pipeline")
coverage = subparsers.add_parser("coverage")
consensus = subparsers.add_parser("consensus")
zero = subparsers.add_parser("zero")
deletion1 = subparsers.add_parser("deletion1")
deletion2 = subparsers.add_parser("deletion2")
deletion3 = subparsers.add_parser("deletion3")
deletionx = subparsers.add_parser("deletionx")
homohetero = subparsers.add_parser("homohetero")
median = subparsers.add_parser("median")
long_coverage = subparsers.add_parser("long_coverage")
barcode = subparsers.add_parser("barcode")

pipeline.add_argument("--gff", \
			help="gff infile")
pipeline.add_argument("--contigsizes", \
			help="Contig sizes for gff parser")
pipeline.add_argument("--feature", \
			help="Feature for gff parser")
pipeline.add_argument("-e", "--snpeff", \
			nargs='*', \
			help="Input options for snpeff, without the '-' before")
pipeline.add_argument("-s", "--snpsift", \
			action="store_true", \
			help="Run snpSift")
pipeline.add_argument("-r", "--clean", \
			action="store_true", \
			help="Removes the SAM and BAM files")
pipeline.add_argument("-p", "--done", \
			action="store_true", \
			help="Add an empty file to mark the directory as done")

# 'coverage' has no unique arguments

consensus.add_argument("-a", "--range", \
			help="somethingsomsing")

# 'zero' has no unique arguments

# 'deletion1' has no unique arguments

# 'deletion2' has no unique arguments

# 'deletion3' has no unique arguments

deletionx.add_argument("-x", "--exons", \
			help="Bed file containing exon coordinates (0-based)")

# 'homohetero' has no unique arguments

median.add_argument("--complex", \
			action="store_true", \
			help="Print full bed output for median")
median.add_argument("--simple", \
			action="store_true", \
			help="Print median coverage only for median")

long_coverage.add_argument("-l", "--limits", \
				type=int, \
				nargs=2, \
				help="Specify lower and upper limits for long_coverage function")

barcode.add_argument("-q", "--quality", \
			type=int, \
			default="20", \
			help="Quality threshold for filtering variants")
barcode.add_argument("--window_size", \
			type=int, \
			default="5000", \
			help="Window size for barcode search")
barcode.add_argument("--primer_size", \
			type=int, \
			default="21", \
			help="Desired size of conserved regions at beginning and end of barcode")

#forward
#reverse
#bamfile
#sortbam
#threads
#outprefix
#verbose
#dev


#ref
#contig
#threshold
#coverage


args = parser.parse_args()


#for i in subparsers.choices.keys(), j in subparsers.choices.values():
#	print(i + "\t" + j.description)

#for i in subparsers.choices.keys():
#	print(i)
#
#print("")

#for j in subparsers.choices.values():
#	print(j.description)

#x = list(subparsers.choices.values())
#y = list(subparsers.choices.keys())

#for i in range(0,len(x)):
#	print(str(y[i]) + "\t" + str(x[i].description))

#print(len(x))
