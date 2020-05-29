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
from functools import reduce
from functools import wraps
import datetime

#######################################################################
# ARGUMENTS
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
			nargs='*', \
			help="Sorted BAM infile (N.B. only the BarcodeSearch function accepts multiple inputs)")
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
parser.add_argument("-o", "--outprefix", \
			help="Output file prefix")

#subparsers = parser.add_subparsers(help="sub-command help")

#group = subparsers.add_parser("Bamparser", parents=[parser], add_help=False, help='Bamparser')
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
			help='Threshold for calculating the coverage percentage; default 20')
group.add_argument("-a", "--range", \
			help="somethingsomsing")
group.add_argument("-x", "--exons", \
			help="Bed file containing exon coordinates (0-based); requires -m")
group.add_argument("-l", "--limits", \
			type=int, \
			nargs=2, \
			help="Specify lower and upper limits for long_coverage function; two arguments required")
group.add_argument("-v", "--verbose", \
			action="store_true", \
			help="Be more verbose")
group.add_argument('--dev', \
			help=argparse.SUPPRESS, action="store_true")

#barcode = subparsers.add_parser("barcodesearch", parents=[parser], add_help=False, help='BarcodeSearch')
barcode = parser.add_argument_group('BarcodeSearch')

barcode.add_argument("--barcode", \
			action="store_true", \
			help="Search the input (sorted) BAM files for suitable barcode regions")
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

args = parser.parse_args()

if not args.coverage:
	if args.feature and args.gff is None:
	        parser.error("--feature requires --gff")
	elif args.gff and args.feature is None:
	        parser.error("--feature requires --gff")

if args.forward and args.reverse and args.ref is None:
	parser.error("--ref [Reference is required]")

#######################################################################
# DEFINE FUNCTIONS USING BAMPARSER MODULE
# 	As the 'bamparse' flags are not tied to a position, must cycle
# 	through all flags
#######################################################################

BamparseList = ["--coverage","--consensus","--zero","--deletion1","--deletion2","--deletion3",\
		"--deletionx","--homohetero","--median","--long_coverage"]

bamparse = False

for item1 in BamparseList:
	for item2 in sys.argv:
		if item1 == item2:
			bamparse = True

#######################################################################
# HANDLING OF MULTIPLE SORTED BAM INPUTS
#	Ensure that if multiple sorted BAM inputs are specified,
#	BarcodeSearch is the function being run
#	Else warn the user and exit
#######################################################################

if args.sortbam:
	if len(args.sortbam) == 1:
		args.sortbam = args.sortbam[0]

	elif len(args.sortbam) > 1 and not args.barcode:
		print("Please note that only BarcodeSearch currently accepts multiple BAM inputs.")
		exit()

#######################################################################
# CHECK DEPENDENCIES
#######################################################################

# DevNote - add checks for snpEff/Java?

def check_samtools():
	try:
		subprocess.check_output('samtools --help', stderr=subprocess.PIPE, shell=True)
	except subprocess.CalledProcessError:
		print("Please ensure that Samtools is in your PATH.")
		exit()
	try:
		subprocess.check_output('samtools depth 2>&1 | grep -- "-aa"', stderr=subprocess.PIPE, shell=True)
	except subprocess.CalledProcessError:
		print("This version of samtools does not support the `depth -aa` option; please update samtools.")
		exit()

def check_bcftools():
	try:
		subprocess.check_output('bcftools --help', stderr=subprocess.PIPE, shell=True)
	except subprocess.CalledProcessError:
		print("Please ensure that BCFtools is in your PATH.")
		exit()

def check_bedtools():
	try:
		subprocess.check_output('bedtools --help', stderr=subprocess.PIPE, shell=True)
	except subprocess.CalledProcessError:
		print("Please ensure that bedtools is in your PATH.")
		exit()

#######################################################################

current_directory = os.getcwd()
name = os.path.basename(current_directory)
add = '../'
add2 = '../Bowtie2/'
threads = str(args.threads)
base = name + '.contigs'
sam = name + '.sam'
bam = name + '.bam'

sorted_bam_out = ""
if args.sortbam:
	sorted_bam_out = add + str(args.sortbam)
else:
	sorted_bam_out = add2 + name + '_sorted.bam'

sorted_bam_bai = name + '_sorted.bam.bai'

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
# BOWTIE2
#	Running bowtie2-build to index reference genome and bowtie2 to align.
#######################################################################

@timing
def bowtie2():
	try:
		subprocess.check_output('bowtie2 --help', stderr=subprocess.PIPE, shell=True)
	except subprocess.CalledProcessError:
		print("Please ensure that Bowtie2 is in your path.")
		exit()

	# Makes new directory 'Bowtie2' if it doesn't exists.
	if not os.path.exists('Bowtie2'):
		os.makedirs('Bowtie2')

	log_file=open('pipeline.log','a')
	# Selected input files using forward and reverse flags,
	# the flags can take several input files.
	file1 = ''
	file2 = ''
	if args.forward:
		f1 = []
		for name in args.forward:
			f1.append(add+name)
		file1 += ','.join(map(str, f1))

	if args.reverse:
		f2 = []
		for name2 in args.reverse:
			f2.append(add+name2)
		file2 += ','.join(map(str, f2))

	# Bowtie2-build, inputs are reference in fasta format and
	# base name for index files, the output are the index files.

	cmd1 = ['bowtie2-build', add+args.ref, base]
	process1 = subprocess.Popen(cmd1, \
		stdout=subprocess.PIPE, \
		stderr = log_file, \
		cwd='Bowtie2')
	while process1.wait() is None:
		pass
	process1.stdout.close()

	# Bowtie2 align step, input are the index files from Bowtie2-build,
	# fastq files (forward and reverse) the output is a SAM file.
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*.rev.1.bt2'):
			cmd2 = ['bowtie2', \
				'-p', threads, \
				'--no-unal', \
				'--very-sensitive', \
				'-x', base, \
				'-1', file1, \
				'-2', file2, \
				'-S', sam]
			process2 = subprocess.Popen(cmd2, \
				stdout=subprocess.PIPE, \
				stderr = log_file, \
				cwd='Bowtie2')
			while process2.wait() is None:
				pass
			process2.stdout.close()
	log_file.close()

#######################################################################
# SAMTOOLS VIEW
#	Converting SAM to BAM using samtools view.
#######################################################################

@timing
def samtools_view():
	check_samtools()
	log_file=open('pipeline.log','a')
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*.sam'):
			cmd3 = ('samtools view -@ %s \
				-b \
				-o %s \
				%s') \
				% (args.threads, bam, sam)
			process3 = subprocess.Popen(cmd3, \
				stdout=subprocess.PIPE, \
				stderr = log_file, \
				shell=True, \
				cwd='Bowtie2')
			while process3.wait() is None:
				pass
			process3.stdout.close()
	log_file.close()

#######################################################################
# SAMTOOLS SORT - FROM BOWTIE PIPELINE
#	Sort BAM files.
#######################################################################

@timing
def samtools_sort():
	check_samtools()
	log_file=open('pipeline.log','a')
#	if glob.glob("Bowtie2/*sorted.bam"):
#		print("Please remove bam files from the Bowtie2 directory before retrying.")
#		exit()
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*.bam'):
			cmd4 = ['samtools', 'sort', \
				'-@', threads, \
				bam, \
				'-o', sorted_bam_out]
			process4 = subprocess.Popen(cmd4, \
				stdout=subprocess.PIPE, \
				stderr = log_file, \
				cwd='Bowtie2')
			while process4.wait() is None:
				pass
			process4.stdout.close()
	log_file.close()

#######################################################################
# SAMTOOLS SORT - FROM BAM INPUT
#	BAM input file by using the '-b' flag.
#######################################################################

@timing
def bam_input():
	check_samtools()
	log_file=open('pipeline.log','a')
#	if glob.glob("Bowtie2/*sorted.bam"):
#		print("Please remove bam files from the Bowtie2 directory before retrying.")
#		exit()
	cmd5 = ['samtools', 'sort', \
		'-@', threads, \
		add+args.bamfile, \
		'-o', sorted_bam_out]
	process5 = subprocess.Popen(cmd5, \
		stdout=subprocess.PIPE, \
		stderr = log_file, \
		cwd='Bowtie2')
	while process5.wait() is None:
		pass
	process5.stdout.close()
	log_file.close()

#######################################################################
# SAMTOOLS INDEX
#	Index sorted BAM files.
#######################################################################

@timing
def samtools_index():
	check_samtools()
	log_file=open('pipeline.log','a')
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*_sorted.bam'):
			cmd6 = ['samtools','index', \
				'-@', threads, \
				sorted_bam_out, \
				sorted_bam_bai]
			process6 = subprocess.Popen(cmd6, \
				stdout=subprocess.PIPE, \
				stderr = log_file, \
				cwd='Bowtie2')
			while process6.wait() is None:
				pass
			process6.stdout.close()
	log_file.close()

#######################################################################
# CLEANUP STEP
#	Remove SAM and BAM files.
#######################################################################

def clean():
	if args.clean:
		for samfile in os.listdir('Bowtie2'):
			if fnmatch.fnmatch(samfile, '*.sam'):
				os.remove('Bowtie2/' + samfile)

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

#######################################################################
# INPUT FILES
#	Define pipeline based on the type of input file
#	Mainly calls the Pipeline module
#######################################################################

def input_files():
	import modules.pipeline as pl
	pl.snpEff_test(args)

	if args.sortbam:
		pl.bcftools(args,threads,sorted_bam_out)
		pl.annotation(args)

	elif args.bamfile:
		bam_input()
		samtools_index()
		pl.bcftools(args,threads,sorted_bam_out)
		pl.annotation(args)
	else:
		bowtie2()
		samtools_view()
		samtools_sort()
		samtools_index()
		pl.bcftools(args,threads,sorted_bam_out)
		pl.annotation(args)

	if args.snpsift:
		pl.snpsift(args)

	if args.clean:
		clean()

	if args.done:
		done()

######################################################################

# DevNote - ensure that there is also a .bai file present

######################################################################
# BAMPARSER
#	Define pipeline based on which Bamparser function is called
######################################################################

def bamparse_func():
	import modules.bamparser as bp

#	input_files()

	if not args.sortbam:
		args.sortbam = "Bowtie2/*sorted.bam" 

	if args.coverage:
		import modules.coverage_stats as cs
		check_samtools()
		if args.gff and not args.outprefix:
			print("If --gff is specified, please ensure that -o is also specified.")
			exit()
		if args.dev:
			import cProfile
			cProfile.runctx('cs.main(args)', globals(), locals())
		else:
			cs.main(args)

	elif args.consensus:
		if args.ref and args.contig and args.range:
			bp.extract_sequence(args)
		else:
			print("Please ensure that a reference [-f], contig [-c] and range [-a] are given.")
			exit()
	elif args.zero:
		import modules.zero_regions as zr
		check_bedtools()
		if args.ref and args.contig:
			if args.dev:
				import cProfile
				cProfile.runctx('zr.main(args)', globals(), locals())
			else:
				zr.main(args)
		else:
			print("Please ensure that a reference [-f] and contig [-c] are given.")
			exit()

	elif args.deletion1 or args.deletion2 or args.deletion3 or args.deletionx or args.homohetero:
		import modules.deletion as dl
		check_samtools()
		if args.deletionx and not args.exons:
			print("Please ensure that a bed file of exons [-x] is given.")
			exit()
		elif args.dev:
			import cProfile
			cProfile.runctx('dl.main(args)', globals(), locals())
		else:
			dl.main(args)

	elif args.median:
		import modules.median_deviation as md
		check_samtools()
		if args.simple or args.complex:
			if args.dev:
				import cProfile
				cProfile.runctx('md.main(args)', globals(), locals())
			else:
				md.main(args)
		else:
			print("Please specify --simple for medians only or --complex for full output")
			exit()

	elif args.long_coverage:
		import modules.coverage_limits as cl
		check_samtools()
		if args.dev:
			import cProfile
			cProfile.runctx('cl.main(args)', globals(), locals())
		else:
			cl.main(args)

	else:
		parser.print_help(sys.stderr)
		exit()

def main():
	if bamparse:
		bamparse_func()

	if args.barcode:
		import modules.barcodesearch as bcs
		check_samtools()
		check_bcftools()
		if args.dev:
			import cProfile
			cProfile.runctx('bcs.barcode(args)', globals(), locals())
		else:
			bcs.barcode(args)
		
	if args.gff and args.feature:
		import modules.pipeline as pl
		try:
			pl.annotation(args)
		except:
			input_files()

	if not bamparse and not args.barcode:
		input_files()

#######################################################################

if __name__ == "__main__":
	main()
