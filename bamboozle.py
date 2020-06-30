#!/usr/bin/env python3


#       Bamboozle v1.1
#       Pipeline that performs bioinformatic analysis including SNP calling 
#       and effect prediction of fastq files or BAM file. 
#
#       Copyright (C) 2018 Vilma Canfjorden. vilma.canfjorden@gmail.com
#       Copyright (C) 2018 Matthew Pinder. matt_pinder13@hotmail.com
#	Copyright (C) 2018 Mats Topel. mats.topel@marine.gu.se
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

#######################################################################
# IMPORTS
#######################################################################

import sys
import os
import argparse
import subprocess
import fnmatch
import glob
from functools import reduce, wraps
from time import time
import datetime

#######################################################################
# ARGUMENTS
#######################################################################

parser = argparse.ArgumentParser(usage="bamboozle.py <command> <args>")

subparsers = parser.add_subparsers(title="Commands", dest="command", metavar="")

# Arguments common to all commands
## Inputs

## DevNote: The block below causes either a help conflict, or the inability to call help when an input file is required...
#input = parser.add_mutually_exclusive_group()
#
#fwd_rev = input.add_argument_group("FASTQ files")
#fwd_rev.add_argument("-F", "--forward", nargs='*', \
#				help="Forward reads")
#fwd_rev.add_argument("-R", "--reverse", nargs='*', \
#				help="Reverse reads")
#
#bam_input = input.add_argument_group("BAM file(s)")
#bam_input.add_argument("-b", "--bamfile", \
#				help="BAM infile")

input_commands = argparse.ArgumentParser(add_help=False)
input_commands.add_argument("-F", "--forward", nargs='*', \
				help="Forward reads; can be a single file or a space-separated list")
input_commands.add_argument("-R", "--reverse", nargs='*', \
				help="Reverse reads; can be a single file or a space-separated list")
input_commands.add_argument("-b", "--bamfile", nargs="*", \
				help="BAM infile")

## Other
other_commands = argparse.ArgumentParser(add_help=False)
other_commands.add_argument("-t", "--threads", default=1, \
				help="Number of threads to use (note: not currently implemented for all functions) (default: 1)")
other_commands.add_argument("-o", "--outprefix", \
				help="Output file prefix")
other_commands.add_argument("-v", "--verbose", action="store_true", \
				help="Be more verbose")
other_commands.add_argument('--dev', \
				help=argparse.SUPPRESS, action="store_true")

# Argument included only in pipeline, consensus, zero, barcode, and lof
ref_command = argparse.ArgumentParser(add_help=False)
ref_command.add_argument("-f", "--ref", \
				help="Reference")

# Argument included only in coverage, consensus, zero, deletion1, deletion2, deletion3, deletionx, homohetero, median, and long_coverage
contig_command = argparse.ArgumentParser(add_help=False)
contig_command.add_argument("-c", "--contig", \
				help="Gives per-contig coverage stats")

# Arguments included only in coverage, deletion1, deletion2, deletion3, deletionx, homohetero, and barcode
threshold_command = argparse.ArgumentParser(add_help=False)
threshold_command.add_argument('-d', '--threshold', type=int, nargs='?', const='1', default='20', \
				help='Threshold for calculating the coverage percentage (default: 20)')

# Arguments included only in pipeline and coverage 
gff_command = argparse.ArgumentParser(add_help=False)
gff_command.add_argument("--gff", \
				help="gff infile")

# Pipeline command
pipeline = subparsers.add_parser("pipeline", parents=[input_commands, other_commands, ref_command, gff_command], \
				usage="bamboozle.py pipeline <args>", \
				help="[Vilma's pipeline]")
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

# Coverage command
coverage = subparsers.add_parser("coverage", parents=[input_commands, other_commands, contig_command, threshold_command, gff_command], \
				usage="bamboozle.py coverage <args>", \
				help="Print a statistic for what percentage of bases in an assembly have >=Nx coverage")

# Consensus command
consensus = subparsers.add_parser("consensus", parents=[input_commands, other_commands, ref_command, contig_command], \
				usage="bamboozle.py consensus <args>", \
				help="Extract the consensus sequence of aligned reads from a region of the reference (WIP)")
consensus.add_argument("-a", "--range", \
			help="somethingsomsing")

# Zero command
zero = subparsers.add_parser("zero", parents=[input_commands, other_commands, ref_command, contig_command], \
				usage="bamboozle.py zero <args>", \
				help="Find areas of zero coverage and print the reference sequence, along with a GC percentage")

# Deletion1 command
deletion1 = subparsers.add_parser("deletion1", parents=[input_commands, other_commands, contig_command, threshold_command], \
				usage="bamboozle.py deletion1 <args>", \
				help="Find deletions")

# Deletion2 command
deletion2 = subparsers.add_parser("deletion2", parents=[input_commands, other_commands, contig_command, threshold_command], \
				usage="bamboozle.py deletion2 <args>", \
				help="Find deletion events")

# Deletion3 command
deletion3 = subparsers.add_parser("deletion3", parents=[input_commands, other_commands, contig_command, threshold_command], \
				usage="bamboozle.py deletion3 <args>", \
				help="Find frameshift deletion events")

# Deletionx command
deletionx = subparsers.add_parser("deletionx", parents=[input_commands, other_commands, contig_command, threshold_command], \
				usage="bamboozle.py deletionx <args>", \
				help="Find deletions occurring within exons")
deletionx.add_argument("-x", "--exons", \
			help="Bed file containing exon coordinates (0-based)")

# Homohetero command
homohetero = subparsers.add_parser("homohetero", parents=[input_commands, other_commands, contig_command, threshold_command], \
				usage="bamboozle.py homohetero <args>", \
				help="Attempt to determine whether a deletion is homozygous or heterozygous")

# Median command
median = subparsers.add_parser("median", parents=[input_commands, other_commands, contig_command], \
				usage="bamboozle.py median <args>", \
				help="Find regions differing from contig median by +/- 50%%, or just contig medians")
median.add_argument("--complex", action="store_true", \
			help="Print full bed output for median")
median.add_argument("--simple", action="store_true", \
			help="Print median coverage only for median")

# Long_coverage command
long_coverage = subparsers.add_parser("long_coverage", parents=[input_commands, other_commands, contig_command], \
				usage="bamboozle.py long_coverage <args>",
				help="Find the longest region between given coverage limits for a given contig")
long_coverage.add_argument("-l", "--limits", type=int, nargs=2, \
			help="Specify lower and upper limits for long_coverage function")

# Barcode command
barcode = subparsers.add_parser("barcode", parents=[input_commands, other_commands, ref_command, threshold_command], \
				usage="bamboozle.py barcode <args>", \
				help="Search the input (sorted) BAM files for suitable barcode regions")
barcode.add_argument("-q", "--quality", type=int, default="20", \
			help="Quality threshold for filtering variants (default: 20)")
barcode.add_argument("--window_size", type=int, default="5000", \
			help="Window size for barcode search (default: 5000)")
barcode.add_argument("--primer_size", type=int, default="21", \
			help="Desired size of conserved regions at beginning and end of barcode (default: 21)")

# SV caller command
sv = subparsers.add_parser("lof", parents=[input_commands, other_commands, ref_command], \
				usage="bamboozle.py lof <args>", \
				help="Run loss-of-function pipeline")
sv.add_argument("--snpeffdb", \
			nargs=1, \
			help="SnpEff database to use")
sv.add_argument("-GFF", "--reference_gff", \
			required = True, \
			help='Reference GFF with gene models for the reference genome')
sv.add_argument("-M", "--masking", \
			required = False, \
			help="Masks found SVs by subtracting SVs found for e.g. additional reads used to correct an assembly")

args = parser.parse_args()

# Get absolute path to Bamboozle directory and assign to args.bamboozledir,
# so that R scripts can be called relative to Bamboozle regardless of cwd

args.bamboozledir = os.path.dirname(os.path.realpath(__file__))

#######################################################################
# DEFINE FUNCTIONS USING BAMPARSER MODULE
#######################################################################

BamparseList = ["coverage","consensus","zero","deletion1","deletion2","deletion3",\
		"deletionx","homohetero","median","long_coverage"]

bamparse = None
if args.command in BamparseList:
	bamparse = True

#######################################################################
# HANDLING BAM FILES
#	First, ensure that all BAM input files are sorted
#	Otherwise, sort them
#	Then, assign all sorted BAMs to args.sortbam
#	Finally, ensure that if multiple BAMs are specified,
#		the barcode command is being run
#		Else warn the user and exit
#######################################################################

#extracting sample name from input BAM, checking if sorted or not

# DevNote - ensure that there is also a .bai file present

#def bam_check(threads, bam_list):
#	args.sortbam = []
#	for bamfile in bam_list:
#		bam_name = os.path.basename(bamfile[:-4])
#		bam_sorted = "%s_sorted.bam" % (bam_name)
#		bam_index = "%s_sorted.bai" % (bam_name)
#
#		#command to check out first line of BAM header and look for "coordinate" (= sorted)
#		cmd1 = "samtools view -H %s | head -n1 | cut -f3 | cut -f2 -d$':'" % (bamfile)
#		proc_1 = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True)
#
#		#if coordinate is present in bam header, bam is sorted
#		std_out, std_error = proc_1.communicate()
#		if std_out.rstrip('\n') == "coordinate":
#			if args.verbose:
#				print("Input BAM " + bamfile + " is already sorted")
#			args.sortbam.append(bamfile)
#		else:
#			print("Input BAM " + bamfile + " is unsorted. Sorting...")
#			cmd2 = "samtools sort -@ %s %s -o %s" % (threads, bamfile,bam_sorted)
#			proc_2 = subprocess.Popen(cmd2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
#			std_out, std_error = proc_2.communicate()
#			cmd3 = "samtools index %s %s" % (bam_sorted, bam_index)
#			proc_3 = subprocess.Popen(cmd3, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
#			std_out, std_error = proc_3.communicate()
#			print("Input BAM " + bamfile + " has been sorted")
#			args.sortbam.append(bam_sorted)
#
#	if len(args.sortbam) == 1:
#		args.sortbam = args.sortbam[0]

#######################################################################
# CHECK DEPENDENCIES
#######################################################################

# DevNote - add checks for snpEff/Java?

def check_samtools():
	try:
		subprocess.check_output('samtools --help', stderr=subprocess.PIPE, shell=True)
	except subprocess.CalledProcessError:
		sys.exit("[Error] Please ensure that Samtools is in your PATH.")
	try:
		subprocess.check_output('samtools depth 2>&1 | grep -- "-aa"', stderr=subprocess.PIPE, shell=True)
	except subprocess.CalledProcessError:
		sys.exit("[Error] This version of samtools does not support the `depth -aa` option; please update samtools.")

def check_bcftools():
	try:
		subprocess.check_output('bcftools --help', stderr=subprocess.PIPE, shell=True)
	except subprocess.CalledProcessError:
		sys.exit("[Error] Please ensure that BCFtools is in your PATH.")

def check_bedtools():
	try:
		subprocess.check_output('bedtools --help', stderr=subprocess.PIPE, shell=True)
	except subprocess.CalledProcessError:
		sys.exit("[Error] Please ensure that bedtools is in your PATH.")

#######################################################################
# TIME DECORATOR
#######################################################################

def timing(function):
	@wraps(function)
	def wrapper(*args, **kwargs):
		now = datetime.datetime.now()
		start = time()
		result = function(*args, **kwargs)
		end = time()
		fh = open("time.log", "a")
		lines_of_text = now.strftime("%Y-%m-%d %H:%M") \
			+ ' Function: ' \
			+ function.__name__ \
			+ ' Elapsed time: {}'.format(end-start) \
			+ ' seconds \n'
		fh.writelines(lines_of_text)
		fh.close()
		return result
	return wrapper

#######################################################################
# CLEANUP STEP
#	Remove SAM and BAM files.
#######################################################################

# DevNote - Add an 'else' statement if --clean has been used with BAM input?
def clean():
	if args.forward and args.reverse:
		for samfile in os.listdir('Bowtie2'):
			if fnmatch.fnmatch(samfile, '*.sam'):
				os.remove('Bowtie2/' + samfile)

		name = os.path.basename(os.getcwd())
		for bamfile in os.listdir('Bowtie2'):
			if fnmatch.fnmatch(bamfile, name + '.bam'):
				os.remove('Bowtie2/' + bamfile)

#######################################################################
# DONE
#	Add empty file when the pipeline is done.
#######################################################################

def done():
	open("pipeline.done", 'a').close()

#######################################################################
# EXIT
#	Exit program.
#######################################################################

def exit():
	sys.exit()

######################################################################


######################################################################
# BAMPARSER
#	Define pipeline based on which Bamparser function is called
######################################################################

def bamparse_func():

#	if not args.sortbam:
#		args.sortbam = "Bowtie2/*sorted.bam" 

	if args.command == "coverage":
		import modules.coverage_stats as cs
		if args.gff and not args.outprefix:
			sys.exit("[Error] If --gff is specified, please ensure that -o is also specified.")
		if args.dev:
			import cProfile
			cProfile.runctx('cs.main(args)', globals(), locals())
		else:
			cs.main(args)

	elif args.command == "consensus":
		import modules.consensus as con
		if args.ref and args.contig and args.range:
			if args.dev:
				import cProfile
				cProfile.runctx('con.main(args)', globals(), locals())
			else:
				con.main(args)
		else:
			sys.exit("[Error] Please ensure that a reference [-f], contig [-c] and range [-a] are given.")

	elif args.command == "zero":
		import modules.zero_regions as zr
		check_bedtools()
		if args.ref and args.contig:
			if args.dev:
				import cProfile
				cProfile.runctx('zr.main(args)', globals(), locals())
			else:
				zr.main(args)
		else:
			sys.exit("[Error] Please ensure that a reference [-f] and contig [-c] are given.")

	elif args.command in ["deletion1", "deletion2", "deletion3", "deletionx", "homohetero"]:
		import modules.deletion as dl
		if args.command == "deletionx" and not args.exons:
			sys.exit("[Error] Please ensure that a bed file of exons [-x] is given.")
		elif args.dev:
			import cProfile
			cProfile.runctx('dl.main(args)', globals(), locals())
		else:
			dl.main(args)

	elif args.command == "median":
		import modules.median_deviation as md
		if args.simple or args.complex:
			if args.dev:
				import cProfile
				cProfile.runctx('md.main(args)', globals(), locals())
			else:
				md.main(args)
		else:
			sys.exit("[Error] Please specify --simple for medians only or --complex for full output")

	elif args.command == "long_coverage":
		import modules.coverage_limits as cl
		if args.dev:
			import cProfile
			cProfile.runctx('cl.main(args)', globals(), locals())
		else:
			cl.main(args)

	else:
		parser.print_help(sys.stderr)
		exit()

def main():
	check_samtools()

	# Ensure that the input into the main pipeline is in sorted BAM format
	import modules.input_files as infiles
	infiles.main(args)

	if args.command == "lof":
		import modules.sv_caller as sv
		sv.main(args, bam_name)

	if bamparse:
		bamparse_func()

	if args.command == "barcode":
		import modules.barcodesearch as bcs
		check_bcftools()
		check_bedtools()
		if args.dev:
			import cProfile
			cProfile.runctx('bcs.barcode(args)', globals(), locals())
		else:
			bcs.barcode(args)

	if args.command == "pipeline":
		if bool(args.feature) != bool(args.gff):
			sys.exit("[Error] --feature requires --gff, and vice versa.")

		import modules.pipeline as pl
		pl.snpEff_test(args)
		pl.bcftools(args)
		pl.annotation(args)

		if args.snpsift:
			pl.snpsift(args)
		if args.clean:
			clean()
		if args.done:
			done()

#######################################################################

if __name__ == "__main__":
	main()
