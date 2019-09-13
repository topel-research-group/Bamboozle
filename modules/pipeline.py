#!/usr/bin/env python3


#	Pipeline that performs bioinformatic analysis including SNP calling 
#	and effect prediction of fastq files or BAM file. 
#
#	Copyright (C) 2018 Vilma Canfjorden. vilma.canfjorden@gmail.com
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <https://www.gnu.org/licenses/>.


import sys
import os 
import argparse
import subprocess
import fnmatch
import glob
from functools import reduce 
from functools import wraps
from time import time
import datetime

#######################################################################

current_directory = os.getcwd()
name = os.path.basename(current_directory)
add = '../'
add2 = '../Bowtie2/'

sorted_bam_bai = name + '_sorted.bam.bai'
bcftools_out = name + '.bcftools_filtered.vcf.gz'
annotated_vcf = name + '.snpeff_annotated.vcf'
annotated_table = name + '.snpsift_table.txt'

#######################################################################

# Time decorator.
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


# Variant calling using bcftools mpileup, input is sorted BAM file, 
# output file is a gzipped vcf file, 
# makes new directory 'Bcftools' if it doesn't exists.
@timing
def bcftools(args,threads,sorted_bam_out):
	log_file=open('pipeline.log','a')
	ref_path = add + str(args.ref)
	if not os.path.exists('Bcftools'):
		os.makedirs('Bcftools')

	cmd7 = ("bcftools mpileup --threads %s -Ou -f %s %s \
		| bcftools call --threads %s -Ou -mv \
		| bcftools filter -s LowQual -e 'QUAL<20' -Oz -o %s") \
	% (threads, ref_path, sorted_bam_out, threads, bcftools_out)
#	% (threads, add+args.ref, sorted_bam_out, threads, bcftools_out)
	process7 = subprocess.Popen(cmd7, \
		stdout=subprocess.PIPE, \
		stderr = log_file, \
		shell=True, \
		cwd='Bcftools')
	while process7.wait() is None:
		pass
	process7.stdout.close()
	log_file.close()
				
# Checks for dependencies required for snpEff. 
def snpEff_test(args):
	# Checks if there is a Skeletonema database, 
	# if it doesn't exists the program will exit 
	# and it has to be created using 'snpEff build'.
	try:
		cmdx = ('snpEff databases | grep "Skeletonema"')
		processx = subprocess.check_output(cmdx, shell=True)
	
	except subprocess.CalledProcessError as e:
		if e.returncode >= 1:
			print('snpEff: Skeletonema database not found, exit program...')
			exit()

	# Try to import gffutils if gff and feature flag is used.	
	if args.gff and args.feature:
		try:
			import gffutils

		except ImportError:
			sys.stderr.write("[Error] The python module \"gffutils\" \
					is not installed\n") 
			sys.stderr.write("[--] Would you like to install it now using \
					'pip install gffutils' [Y/N]?\n")
			answer = sys.stdin.readline()
			if answer[0].lower() == "y":
				sys.stderr.write("[--] Running \"pip install gffutils\"\n")
				from subprocess import call
				call(["pip", "install", "gffutils"])
			else:
				sys.exit("[Error] Exiting due to missing dependency \"gffutils\"")

# Annotating variant calling output using snpEff, output is a vcf, 
# the vcf file is bgzipped to work as an input file to the Fst analysis,
# the original vcf file is kept by using the -c flag. 
@timing
def annotation(args):					
	log_file=open('pipeline.log','a')
	for file in os.listdir('Bcftools'):
		if fnmatch.fnmatch(file, '*.bcftools_filtered.vcf.gz'):
			my_output = annotated_vcf
			my_interval = ""
			if args.gff and args.feature:
				from modules.parse_gff_2 import main as parse
				parse(args.gff, args.feature, args.contigsizes)
				out = add + 'out.gff'	
				my_interval = "-fi %s" % out
				my_output = name + '_' + args.feature + '_snpeff_annotated.vcf'

			# If you want to specify options yourself.
			if args.snpeff:
				opt = '-'+' -'.join(args.snpeff)
				my_args = my_interval + " " + opt \
					+ " Skeletonema_marinoi_v1.1.1.1 \
					 -stats snpEff_summary.html"
			else:
				my_args = my_interval + \
					" Skeletonema_marinoi_v1.1.1.1 \
					-stats snpEff_summary.html"

			cmd8 = ("snpEff	%s %s > %s") \
				% (my_args, bcftools_out, my_output)
			process8 = subprocess.Popen(cmd8, \
				stdout=subprocess.PIPE, \
				stderr = log_file, \
				shell=True, \
				cwd='Bcftools')
			while process8.wait() is None:
				pass
			process8.stdout.close()

			cmd9 = ('bgzip -c %s > %s') \
				% (my_output, my_output + '.gz')
			process9 = subprocess.Popen(cmd9, \
				stdout=subprocess.PIPE, \
				stderr = log_file, \
				shell=True, \
				cwd='Bcftools')
			while process9.wait() is None:
				pass
			process9.stdout.close()

			# Add headers to gff parsed vcf file for fst statistics.
			if args.gff and args.feature:
				my_output_hdr = name + '_' + args.feature + '_hdr_snpeff_annotated.vcf'
				cmd_a = ("bcftools view -h %s > hdr.txt") \
					% (my_output)
				process_a = subprocess.Popen(cmd_a, \
					stdout=subprocess.PIPE, \
					stderr = log_file, \
					shell=True, \
					cwd='Bcftools')
				while process_a.wait() is None:
					pass
				process_a.stdout.close()

				cmd_b = ('''sed -i '/##INFO=<ID=MQ/a##INFO=<ID=out_ID,Number=1,Type=String,Description="none">\
					\\n##INFO=<ID=out_Parent,Number=1,Type=String,Description="none">\
					\\n##INFO=<ID=out_type,Number=1,Type=String,Description="none">\
					\\n##INFO=<ID=out_source,Number=1,Type=String,Description="none">' \
					hdr.txt''') 
				process_b = subprocess.Popen(cmd_b, \
					stdout=subprocess.PIPE, \
					stderr = log_file, \
					shell=True, \
					cwd='Bcftools')
				while process_b.wait() is None:
					pass
				process_b.stdout.close()
				
				cmd_c = ("bcftools reheader -h hdr.txt %s > %s") \
					% (my_output, my_output_hdr)
				process_c = subprocess.Popen(cmd_c, \
					stdout=subprocess.PIPE, \
					stderr = log_file, \
					shell=True, \
					cwd='Bcftools')
				while process_c.wait() is None:
					pass
				process_c.stdout.close()

				cmd_d = ('bgzip -c %s > %s') \
					% (my_output_hdr, my_output_hdr + '.gz')
				process_d = subprocess.Popen(cmd_d, \
					stdout=subprocess.PIPE, \
					stderr = log_file, \
					shell=True, \
					cwd='Bcftools')
				while process_d.wait() is None:
					pass
				process_d.stdout.close()
				os.remove('out.gff')
				os.remove('Bcftools/hdr.txt')
				for vcffile in os.listdir('Bcftools'):
					if fnmatch.fnmatch(vcffile, '*'+args.feature+'_snpeff_annotated.vcf'+'*'):
						os.remove('Bcftools/' + vcffile)
			else:
				pass
	log_file.close()
					

# Filtering and making a summary of annotated files using 
# the vcf (not bgzipped) output file from snpEff, 
# the summary will be in table format, tab separated.
@timing
def snpsift(args):
	for file in os.listdir('Bcftools'):
		if fnmatch.fnmatch(file, '*_annotated.vcf'):
			my_output = annotated_vcf
			if args.gff and args.feature:
				my_output = name + '_' + args.feature + '_hdr_snpeff_annotated.vcf'

			cmd10 = ('java -jar /usr/local/packages/snpEff/SnpSift.jar \
			extractFields -e "." -s "," %s CHROM POS "EFF[*].GENE" REF ALT QUAL DP AF \
			"EFF[*].EFFECT" "EFF[*].AA" "EFF[*].FUNCLASS" > %s') \
			% (my_output, annotated_table) 
			process10 = subprocess.Popen(cmd10, \
				stdout=subprocess.PIPE, \
				shell=True, \
				cwd='Bcftools')
			while process10.wait() is None:
				pass
			process10.stdout.close()

