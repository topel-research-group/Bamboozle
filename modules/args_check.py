#!/usr/bin/env python3

import argparse
import sys

###################################################################################################
# CAN AVOID HAVING TO HARDCODE THE DESCRIPTIONS INTO METAVAR BY HAVING
#	add_help=False 
# ON THE PARENT PARSERS, THEREBY ONLY SHOWING THE HELP TEXT FOR THE SUBPARSERS
###################################################################################################

# Can information be obtained from parser._subparsers about the parser names?
# subparsers.choices is a dictionary whose keys are the subprocesses
# for i in subparsers.choices.keys():
#	print(i)
# Incorporate this into metavars?

#def metavar_mod(metavar):
#	x = list(metavar.values())
#	y = list(metavar.keys())
#	for i in range(0,len(x)):
#		yield(str(y[i]) + "\t" + str(x[i].description))

#subparsers = parser.add_subparsers(title="Commands", \
#					metavar=metavar_mod(subparsers.choices))

###################################################################################################

parser = argparse.ArgumentParser(usage="bamboozle.py <command> <args>")

subparsers = parser.add_subparsers(title="Commands", metavar="")

all_commands = argparse.ArgumentParser(add_help=False)

all_commands.add_argument("-F", "--forward", nargs='*', \
				help="Forward reads")
all_commands.add_argument("-R", "--reverse", nargs='*', \
				help="Reverse reads")
all_commands.add_argument("-b", "--bamfile", \
				help="BAM infile")
all_commands.add_argument("--sortbam", nargs='*', \
				help="Sorted BAM infile (N.B. only the BarcodeSearch function accepts multiple inputs)")
all_commands.add_argument("-t", "--threads", default=1, \
				help="Number of threads to use (note: not currently implemented for all functions)")
all_commands.add_argument("-o", "--outprefix", help="Output file prefix")
all_commands.add_argument("-v", "--verbose", action="store_true", \
				help="Be more verbose")
all_commands.add_argument('--dev', \
				help=argparse.SUPPRESS, action="store_true")

ref_command = argparse.ArgumentParser(add_help=False)
ref_command.add_argument("-f", "--ref", \
				help="Reference")

contig_command = argparse.ArgumentParser(add_help=False)
contig_command.add_argument("-c", "--contig", \
				help="Gives per-contig coverage stats")

threshold_command = argparse.ArgumentParser(add_help=False)
threshold_command.add_argument('-d', '--threshold', type=int, nargs='?', const='1', default='20', \
				help='Threshold for calculating the coverage percentage; default 20')

pipeline = subparsers.add_parser("pipeline", parents=[all_commands, ref_command], \
					usage="bamboozle.py pipeline <args>", \
					help="[Vilma's pipeline]")
pipeline.add_argument("--gff", \
			help="gff infile")
pipeline.add_argument("--contigsizes", \
			help="Contig sizes for gff parser")
pipeline.add_argument("--feature", \
			help="Feature for gff parser")
pipeline.add_argument("-e", "--snpeff", nargs='*', \
			help="Input options for snpeff, without the '-' before")
pipeline.add_argument("-s", "--snpsift", action="store_true", \
			help="Run snpSift")
pipeline.add_argument("-r", "--clean", action="store_true", \
			help="Removes the SAM and BAM files")
pipeline.add_argument("-p", "--done", action="store_true", \
			help="Add an empty file to mark the directory as done")

coverage = subparsers.add_parser("coverage", parents=[all_commands, contig_command], \
					usage="bamboozle.py coverage <args>", \
					help="Print a statistic for what percentage of bases in an assembly have >=Nx coverage")

consensus = subparsers.add_parser("consensus", parents=[all_commands, ref_command, contig_command], \
					usage="bamboozle.py consensus <args>", \
					help="Extract the consensus sequence of aligned reads from a region of the reference (WIP)")
consensus.add_argument("-a", "--range", \
			help="somethingsomsing")

zero = subparsers.add_parser("zero", parents=[all_commands, ref_command, contig_command], \
				usage="bamboozle.py zero <args>", \
				help="Find areas of zero coverage and print the reference sequence, along with a GC percentage")

deletion1 = subparsers.add_parser("deletion1", parents=[all_commands, contig_command, threshold_command], \
					usage="bamboozle.py deletion1 <args>", \
					help="Find deletions")

deletion2 = subparsers.add_parser("deletion2", parents=[all_commands, contig_command, threshold_command], \
					usage="bamboozle.py deletion2 <args>", \
					help="Find deletion events")

deletion3 = subparsers.add_parser("deletion3", parents=[all_commands, contig_command, threshold_command], \
					usage="bamboozle.py deletion3 <args>", \
					help="Find frameshift deletion events")

deletionx = subparsers.add_parser("deletionx", parents=[all_commands, contig_command, threshold_command], \
					usage="bamboozle.py deletionx <args>", \
					help="Find deletions occurring within exons")
deletionx.add_argument("-x", "--exons", \
			help="Bed file containing exon coordinates (0-based)")

homohetero = subparsers.add_parser("homohetero", parents=[all_commands, contig_command, threshold_command], \
					usage="bamboozle.py homohetero <args>", \
					help="Attempt to determine whether a deletion is homozygous or heterozygous")

median = subparsers.add_parser("median", parents=[all_commands, contig_command], \
					usage="bamboozle.py median <args>", \
					help="Find regions differing from contig median by +/- 50%%, or just contig medians")
median.add_argument("--complex", action="store_true", \
			help="Print full bed output for median")
median.add_argument("--simple", action="store_true", \
			help="Print median coverage only for median")

long_coverage = subparsers.add_parser("long_coverage", parents=[all_commands, contig_command], \
					usage="bamboozle.py long_coverage <args>",
					help="Find the longest region between given coverage limits for a given contig")
long_coverage.add_argument("-l", "--limits", type=int, nargs=2, \
			help="Specify lower and upper limits for long_coverage function")

barcode = subparsers.add_parser("barcode", parents=[all_commands, ref_command], \
					usage="bamboozle.py barcode <args>", \
					help="Search the input (sorted) BAM files for suitable barcode regions")
barcode.add_argument("-q", "--quality", type=int, default="20", \
			help="Quality threshold for filtering variants")
barcode.add_argument("--window_size", type=int, default="5000", \
			help="Window size for barcode search")
barcode.add_argument("--primer_size", type=int, default="21", \
			help="Desired size of conserved regions at beginning and end of barcode")


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
