#!/usr/bin/env python3


#       Bamboozle v1.1
#       Pipeline that performs bioinformatic analysis including SNP calling 
#       and effect prediction of fastq files or BAM file. 
#
#       Copyright (C) 2018 Vilma Canfjorden. vilma.canfjorden@gmail.com
#       Copyright (C) 2018 Matthew Pinder. matt_pinder13@hotmail.com
#       Copyright (C) 2020 André Soares. andre.soares@bioenv.gu.se
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
import subprocess
import fnmatch
from functools import reduce, wraps
from time import time
import datetime

#######################################################################
# GLOBAL VARIABLES
#######################################################################

current_directory = os.getcwd()
name = os.path.basename(current_directory)
add = '../'
add2 = '../Bowtie2/'
base = name + '.contigs'
sam = name + '.sam'
bam = name + '.bam'

sorted_bam_out = add2 + name + '_sorted.bam'
sorted_bam_bai = add2 + name + '_sorted.bam.bai'

#######################################################################
# HANDLING MULTIPLE INPUT FILES
#	If input is given in a comma-separated list, convert it to
#	a list rather than a string
#######################################################################

# DevNote - this appears not to work currently...

def list_format(args):
	for name in ["args.forward", "args.reverse", "args.bamfile"]:
		if name and "," in name:
			sys.exit("[Error] FASTQ/BAM input must be space-separated, not comma-separated")

#######################################################################
# HANDLING BAM FILES
#	First, ensure that all BAM input files are sorted
#	Otherwise, sort them
#	Then, assign all sorted BAMs to args.sortbam
#	Finally, ensure that if multiple BAMs are specified,
#		the barcode command is being run
#		Else warn the user and exit
#######################################################################

# Extracting sample name from input BAM, checking if sorted or not

# DevNote - ensure that there is also a .bai file present

# DevNote - combine this with the samtools sort/index commands below, rather than having the same function run twice

def bam_check(args):
	args.sortbam = []
	for bamfile in args.bamfile:
		bam_name = os.path.basename(bamfile[:-4])
		bam_sorted = "Bowtie2/%s_sorted.bam" % (bam_name)
		bam_index = "Bowtie2/%s_sorted.bai" % (bam_name)

		#command to check out first line of BAM header and look for "coordinate" (= sorted)
		cmd1 = ["samtools", "view", "-H", bamfile]
		proc_1 = subprocess.Popen(cmd1, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True)

		#if coordinate is present in bam header, bam is sorted
		std_out, std_error = proc_1.communicate()
		if "coordinate" in std_out.split("\n")[0]:
			if args.verbose:
				print("Input BAM " + bamfile + " is already sorted")
			args.sortbam.append(bamfile)
		else:
			print("Input BAM " + bamfile + " is unsorted. Sorting...")

			# Makes new directory 'Bowtie2' if it doesn't exists.
			if not os.path.exists('Bowtie2'):
				os.makedirs('Bowtie2')

			cmd2 = ["samtools", "sort", "-@", str(args.threads), bamfile, "-o", bam_sorted]
			proc_2 = subprocess.Popen(cmd2, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
			std_out, std_error = proc_2.communicate()
			cmd3 = ["samtools", "index", bam_sorted, bam_index]
			proc_3 = subprocess.Popen(cmd3, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
			std_out, std_error = proc_3.communicate()
			print("Input BAM " + bamfile + " has been sorted")
			args.sortbam.append(bam_sorted)

	if len(args.sortbam) == 1:
		args.sortbam = args.sortbam[0]

	return(args.sortbam)

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
def bowtie2(args):
	try:
		subprocess.check_output(['bowtie2', '--help'])
	except FileNotFoundError:
		sys.exit("[Error] Please ensure that Bowtie2 is in your path.")

	# Ensure that a reference is provided
	if not args.ref:
		sys.exit("[Error] A reference fasta (--ref) must be provided.")

	# Makes new directory 'Bowtie2' if it doesn't exists.
	if not os.path.exists('Bowtie2'):
		os.makedirs('Bowtie2')

	log_file=open('pipeline.log','a')
	# Selected input files using forward and reverse flags,
	# the flags can take several input files.

	file1 = ''
	file2 = ''
	if args.forward:
		if isinstance(args.forward, str):
			file1 = str(add + args.forward)
		else:
			f1 = []
			for name in args.forward:
				f1.append(add+name)
			file1 += ','.join(map(str, f1))

	if args.reverse:
		if isinstance(args.reverse, str):
			file2 = str(add + args.reverse)
		else:
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
				'-p', str(args.threads), \
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
def samtools_view(args):
	log_file=open('pipeline.log','a')
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*.sam'):
			cmd3 = ["samtools", "view", "-@", str(args.threads), "-b", "-o", bam, sam]
			process3 = subprocess.Popen(cmd3, \
				stdout=subprocess.PIPE, \
				stderr = log_file, \
				shell=False, \
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
def samtools_sort(args):
	log_file=open('pipeline.log','a')
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*.bam'):
			cmd4 = ['samtools', 'sort', \
				'-@', str(args.threads), \
				bam, \
				'-o', sorted_bam_out]
			process4 = subprocess.Popen(cmd4, \
				stdout=subprocess.PIPE, \
				stderr = log_file, \
				cwd='Bowtie2')
			while process4.wait() is None:
				pass
			process4.stdout.close()
		args.sortbam = str('Bowtie2/' + name + '_sorted.bam')
	log_file.close()

#######################################################################
# SAMTOOLS INDEX
#	Index sorted BAM files.
#######################################################################

@timing
def samtools_index(args):
	log_file=open('pipeline.log','a')
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*_sorted.bam'):
			cmd6 = ['samtools','index', \
				'-@', str(args.threads), \
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
# SAMTOOLS PHASE
#	Phase alleles for barcoding purposes
#	DevNote - THIS IS FREEZING UP! Error Error
#######################################################################

@timing
def phasing(bamfile, fileprefix, threads):
	log_file=open('pipeline.log','a')
	cmd7 = ['samtools', 'phase', \
		'-b', fileprefix, \
		bamfile]
	process7 = subprocess.Popen(cmd7, \
		stdout=subprocess.DEVNULL, \
		stderr=log_file, \
		shell=False)
	while process7.wait() is None:
		pass

	for phasefile in [fileprefix + ".0", fileprefix + ".1"]:
		cmd8 = ['samtools', 'index', \
			'-@', str(threads), \
			str(phasefile + ".bam"), \
			str(phasefile + ".bai")]
		process8 = subprocess.Popen(cmd8, \
			stdout=subprocess.DEVNULL, \
			stderr=log_file, \
			shell=False)
		while process8.wait() is None:
			pass

	log_file.close()

#######################################################################
# BCFTOOLS
#	Generate VCF files for phased files
#	DevNote - should be combined with the version in barcodesearch
#######################################################################

@timing
def phasevcf(mainfile, allele0, allele1, reference, threads, qual):
	log_file=open('pipeline.log','a')
	quality = "%QUAL>" + str(qual)
	for infile in [mainfile, allele0, allele1]:
		invcf = infile.replace(".bam", ".vcf.gz")
		if not os.path.isfile(invcf):
			cmdA = ["bcftools", "mpileup", "--threads", str(threads), "--fasta-ref", reference, infile]
			procA = subprocess.Popen(cmdA, stdout=subprocess.PIPE, stderr=log_file, shell=False)

			cmdB = ["bcftools", "call", "--threads", str(threads), "-mv"]
			procB = subprocess.Popen(cmdB, stdin=procA.stdout, stdout=subprocess.PIPE, stderr=log_file, shell=False)

			cmdC = ["bcftools", "filter", "--threads", str(threads), "-i", quality, "-Oz", "-o", invcf]
			procC = subprocess.Popen(cmdC, stdin=procB.stdout, stdout=subprocess.PIPE, stderr=log_file, shell=False)
			while procC.wait() is None:
				pass
			procC.stdout.close()

			cmdD = ["bcftools", "index", "-f", "--threads", str(threads), invcf]
			procD = subprocess.Popen(cmdD, stdout=subprocess.PIPE, stderr=log_file, shell=False)
			while procD.wait() is None:
				pass
			procD.stdout.close()

	log_file.close()

#######################################################################
# MAIN
#	Check the input arguments:
#		If FASTQ - map, sort and index
#		If BAM - sort and index if unsorted
#	Then add to and return args.sortbam
#######################################################################

def main(args):
	# Ensure list inputs are in the correct format
	list_format(args)

	# If FASTQ input is given, generate a sorted BAM file
	if args.forward and args.reverse:
		if args.command == "barcode":
			sys.exit("[Error] BarcodeSearch currently doesn't accept FASTQ input; please give BAM input instead.")
		else:
			bowtie2(args)
			samtools_view(args)
			samtools_sort(args)
			samtools_index(args)

	# If BAM input is given, ensure it's sorted; otherwise sort and index it
	elif args.bamfile:
		if len(args.bamfile) > 1 and args.command != "barcode":
			sys.exit("[Error] Please note that only BarcodeSearch currently accepts multiple BAM inputs.")
		elif len(args.bamfile) == 1 and args.command == "barcode":
			sys.exit("[Error] Please note that BarcodeSearch requires multiple BAM inputs.")
		else:
			bam_check(args)

	# If BarcodeSearch is being run, check whether samtools phase has been run; if not, run it
	# The functions below don't work! Pregenerate files for testing

	if args.command == "barcode":
		for infile in args.sortbam:
			noext = os.path.splitext(infile)[0]
			phase0 = noext + ".0.bam"
			phase1 = noext + ".1.bam"
			mainvcf = noext + ".vcf.gz.csi"
			phase0vcf = noext + ".0.vcf.gz.csi"
			phase1vcf = noext + ".1.vcf.gz.csi"

			if args.ploidy == "haploid":
				if not os.path.isfile(mainvcf):
					sys.exit("[Error] Please ensure that you run 'bcftools index' on all input files. This will be automated in future versions of Bamboozle.")
			if args.ploidy == "diploid":
				if not (os.path.isfile(phase0) and os.path.isfile(phase1)):
					print("Generating phased BAM files")
					phasing(infile, noext, args.threads)
				if not (os.path.isfile(mainvcf) and os.path.isfile(phase0vcf) and os.path.isfile(phase1vcf)):
					print("Generating phased VCF files")
					phasevcf(infile, phase0, phase1, args.ref, args.threads, args.quality)
	return(args.sortbam)

#######################################################################

if __name__ == "__main__":
	main()
