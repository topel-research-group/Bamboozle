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

parser = argparse.ArgumentParser(usage="bamboozle.py <command> <args>")

subparsers = parser.add_subparsers(title="Commands", dest="command", metavar="")

# Arguments common to all commands
## Inputs

## DevNote: The block below causes either a help conflict, or the inability to call help when an input file is required...
#infiles = parser.add_mutually_exclusive_group()
#
#fwd_rev = infiles.add_argument_group("FASTQ files")
#fwd_rev.add_argument("-F", "--forward", nargs='*', \
#				help="Forward reads")
#fwd_rev.add_argument("-R", "--reverse", nargs='*', \
#				help="Reverse reads")
#
#bam_input = infiles.add_argument_group("BAM file(s)")
#bam_input.add_argument("-b", "--bamfile", \
#				help="BAM infile")
#bam_input.add_argument("--sortbam", nargs='*', \
#				help="Sorted BAM infile (N.B. only the BarcodeSearch function accepts multiple inputs)")

input_commands = argparse.ArgumentParser(add_help=False)
input_commands.add_argument("-F", "--forward", nargs='*', \
				help="Forward reads")
input_commands.add_argument("-R", "--reverse", nargs='*', \
				help="Reverse reads")
input_commands.add_argument("-b", "--bamfile", \
				help="BAM infile")
input_commands.add_argument("--sortbam", nargs='*', \
				help="Sorted BAM infile (N.B. only the BarcodeSearch function accepts multiple inputs)")

## Other
other_commands = argparse.ArgumentParser(add_help=False)
other_commands.add_argument("-t", "--threads", default=1, \
				help="Number of threads to use (note: not currently implemented for all functions)")
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

# Arguments included only in deletion1, deletion2, deletion3, deletionx, and homohetero
threshold_command = argparse.ArgumentParser(add_help=False)
threshold_command.add_argument('-d', '--threshold', type=int, nargs='?', const='1', default='20', \
				help='Threshold for calculating the coverage percentage; default 20')

# Pipeline command
pipeline = subparsers.add_parser("pipeline", parents=[input_commands, other_commands, ref_command], \
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

# Coverage command
coverage = subparsers.add_parser("coverage", parents=[input_commands, other_commands, contig_command], \
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
barcode = subparsers.add_parser("barcode", parents=[input_commands, other_commands, ref_command], \
				usage="bamboozle.py barcode <args>", \
				help="Search the input (sorted) BAM files for suitable barcode regions")
barcode.add_argument("-q", "--quality", type=int, default="20", \
			help="Quality threshold for filtering variants")
barcode.add_argument("--window_size", type=int, default="5000", \
			help="Window size for barcode search")
barcode.add_argument("--primer_size", type=int, default="21", \
			help="Desired size of conserved regions at beginning and end of barcode")

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

if args.command == "pipeline":
	if args.feature and args.gff is None:
	        parser.error("--feature requires --gff")
	elif args.gff and args.feature is None:
	        parser.error("--feature requires --gff")

if args.forward and args.reverse and args.ref is None:
	parser.error("--ref [Reference is required]")

#######################################################################
# DEFINE FUNCTIONS USING BAMPARSER MODULE
#######################################################################

BamparseList = ["--coverage","--consensus","--zero","--deletion1","--deletion2","--deletion3",\
		"--deletionx","--homohetero","--median","--long_coverage"]

bamparse = None
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

# DevNote - this will need changing after the removal of the sortbam option

#if args.sortbam:
#	if len(args.sortbam) == 1:
#		args.sortbam = args.sortbam[0]
#	elif len(args.sortbam) > 1 and args.command != "barcode":
#		print("Please note that only BarcodeSearch currently accepts multiple BAM inputs.")
#		exit()

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
#if args.sortbam:
#	sorted_bam_out = add + str(args.sortbam)
#else:
#	sorted_bam_out = add2 + name + '_sorted.bam'
#
#sorted_bam_bai = name + '_sorted.bam.bai'

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

#extracting sample name from input BAM, checking if sorted or not
bam_name = args.bamfile[:-4]
bam_sorted = "%s_sorted.bam" % (bam_name)
bam_index = "%s_sorted.bai" % (bam_name)

def bam_check(bamfile,bam_sorted,bam_index):
        #command to check out first line of BAM header and look for "coordinate" (= sorted)
        cmd1 = "samtools view -H %s | head -n1 | cut -f3 | cut -f2 -d$':'" % (bamfile)
        proc_1 = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True)

        #if coordinate is present in bam header, bam is sorted
        std_out, std_error = proc_1.communicate()
        if std_out.rstrip('\n') == "coordinate":
                print("Input BAM was already sorted")
                args.sortbam = args.bamfile
        else:
                cmd2 = "samtools sort %s -o %s" % (bamfile,bam_sorted)
                proc_2 = subprocess.Popen(cmd2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                std_out, std_error = proc_2.communicate()
                cmd3 = "samtools index %s %s" % (bam_sorted, bam_index)
                proc_3 = subprocess.Popen(cmd3, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                std_out, std_error = proc_3.communicate()
                print("Input BAM has been sorted")
                args.sortbam = bam_sorted

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
		bam_check(args.bamfile, bam_sorted, bam_index)
		pl.bcftools(args,threads,sorted_bam_out)
		pl.annotation(args)

#	elif args.bamfile:
#		bam_input()
#		samtools_index()
#		pl.bcftools(args,threads,sorted_bam_out)
#		pl.annotation(args)
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

#	if not args.sortbam:
#		args.sortbam = "Bowtie2/*sorted.bam" 

	if args.command == "coverage":
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

	elif args.command == "consensus":
		import modules.consensus as con
		if args.ref and args.contig and args.range:
			if args.dev:
				import cProfile
				cProfile.runctx('con.main(args)', globals(), locals())
			else:
				con.main(args)
		else:
			print("Please ensure that a reference [-f], contig [-c] and range [-a] are given.")
			exit()

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
			print("Please ensure that a reference [-f] and contig [-c] are given.")
			exit()

	elif args.command in ["deletion1", "deletion2", "deletion3", "deletionx", "homohetero"]:
		import modules.deletion as dl
		check_samtools()
		if args.command == "deletionx" and not args.exons:
			print("Please ensure that a bed file of exons [-x] is given.")
			exit()
		elif args.dev:
			import cProfile
			cProfile.runctx('dl.main(args)', globals(), locals())
		else:
			dl.main(args)

	elif args.command == "median":
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

	elif args.command == "long_coverage":
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
	if args.command == "lof":
		import modules.sv_caller as sv
		check_samtools()
		sv.main(args, bam_name)

	if bamparse:
		bam_check(args.bamfile, bam_sorted, bam_index)
		bamparse_func()

	if args.command == "barcode":
		import modules.barcodesearch as bcs
		check_samtools()
		check_bcftools()
		bam_check(args.bamfile, bam_sorted, bam_index)
		if args.dev:
			import cProfile
			cProfile.runctx('bcs.barcode(args)', globals(), locals())
		else:
			bcs.barcode(args)

	if args.command == "pipeline":		
		import modules.pipeline as pl
		try:
			pl.annotation(args)
		except:
			input_files()

	# If the command is not a 'bamparse', barcode or lof command, run [Vilma's pipeline]
	if not bamparse and args.command not in ["barcode", "lof"]:
		input_files()
		#not sure if this is where it should go?
		bam_check(args.bamfile, bam_sorted, bam_index)

#######################################################################

if __name__ == "__main__":
	main()
