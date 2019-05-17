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
from functools import reduce 
from functools import wraps
from time import time
import datetime

#######################################################################

parser = argparse.ArgumentParser(prog="ADD-SCRIPT-NAME-HERE")
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

parser.add_argument('--bamparse', action="store_true", help ="Run bamparser")
parser.add_argument('--coverage', action="store_true", help='Print a statistic for what percentage of bases in an assembly have >=Nx coverage')
parser.add_argument('--consensus', action="store_true", help='Extract the consensus sequence of aligned reads from a specific region of the reference sequence (WIP)')
parser.add_argument('--zero', action="store_true", help='Find areas of zero coverage and print the reference sequence, along with a GC percentage')
parser.add_argument('--deletion1', action="store_true", help='Find deletions')
parser.add_argument('--deletion2', action="store_true", help='Find deletion events')
parser.add_argument('--deletion3', action="store_true", help='Find frameshift deletion events')
parser.add_argument('--deletionx', action="store_true", help='Find deletions occurring within exons')
parser.add_argument('--homohetero', action="store_true", help='Attempt to determine whether a deletion is homozygous or heterozygous')
parser.add_argument('--median', action="store_true", help='Find regions differing from contig median by +/- 50%%, or just contig medians')
parser.add_argument('--long_coverage', action="store_true", help='Find the longest region between given coverage limits for a given contig')

parser.add_argument('--complex', action="store_true", help='Print full bed output for median')
parser.add_argument('--simple', action="store_true", help='Print median coverage only for median')
parser.add_argument('-c', '--contig', help='Gives per-contig coverage stats')
parser.add_argument('-d', '--threshold', type=int, nargs='?', const='1', default='20', help='Threshold for calculating coverage percentage; default 20')
parser.add_argument("-a", "--range", help="somethingsomsing")
parser.add_argument("-m", "--mutations", help="List of mutation events; currently requires output from bamboozle deletion function")
parser.add_argument("-x", "--exons", help="Bed file containing exon coordinates (0-based); -m also required")
parser.add_argument("-l", "--limits", type=int, nargs=2, help="Specify lower and upper limits for long_coverage function; two arguments required")
parser.add_argument("-v", "--verbose", action="store_true", help="Be more verbose")
parser.add_argument('--dev', help=argparse.SUPPRESS, action="store_true")

args = parser.parse_args()

if args.feature and args.gff is None:
        parser.error("--feature requires --gff")
elif args.gff and args.feature is None:
        parser.error("--feature requires --gff")

#######################################################################

current_directory = os.getcwd()
name = os.path.basename(current_directory)
add = '../'
add2 = '../Bowtie2/'
threads = str(args.threads)
base = name + '.contigs'
sam = name + '.sam'
bam = name + '.bam' 

if args.sortbam:
	sorted_bam_out = add + str(args.sortbam)
else:
	sorted_bam_out = add2 + name + '_sorted.bam'

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

# Makes new directory 'Bowtie2' if it doesn't exists.
if not os.path.exists('Bowtie2'):
	os.makedirs('Bowtie2')

# Running bowtie2-build to index reference genome and bowtie2 to align. 
@timing
def bowtie2(args):
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
			
# Converting SAM to BAM using samtools view.
@timing
def samtools_view():
	log_file=open('pipeline.log','a')
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*.sam'):
			cmd3 = ('samtools view \
				-@ %s \
				-Sb %s \
				> %s') \
				% (args.threads, sam, bam)
			process3 = subprocess.Popen(cmd3, \
				stdout=subprocess.PIPE, \
				stderr = log_file, \
				shell=True, \
				cwd='Bowtie2')
			while process3.wait() is None:
				pass
			process3.stdout.close()
	log_file.close()

# Sort BAM files.
@timing
def samtools_sort():
	log_file=open('pipeline.log','a')
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*.bam'):
			cmd4 = ['samtools', 'sort', \
				'-@', '$NSLOTS', \
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
		
# BAM input file by using the '-b' flag.
@timing
def bam_input(args):
	log_file=open('pipeline.log','a')
	cmd5 = ['samtools', 'sort', \
		'-@', '$NSLOTS', \
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

# Index sorted BAM files.
@timing
def samtools_index():
	log_file=open('pipeline.log','a')
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*_sorted.bam'):
			cmd6 = ['samtools','index', \
				sorted_bam_out, sorted_bam_bai]
			process6 = subprocess.Popen(cmd6, \
				stdout=subprocess.PIPE, \
				stderr = log_file, \
				cwd='Bowtie2')
			while process6.wait() is None:
				pass
			process6.stdout.close()
	log_file.close()

# Remove SAM and BAM files.
def clean():
	if args.clean:
		for samfile in os.listdir('Bowtie2'):
			if fnmatch.fnmatch(samfile, '*.sam'):
				os.remove('Bowtie2/' + samfile)
 
		for bamfile in os.listdir('Bowtie2'):
			if fnmatch.fnmatch(bamfile, name + '.bam'):
				os.remove('Bowtie2/' + bamfile)

# Variant calling using bcftools mpileup, input is sorted BAM file, 
# output file is a gzipped vcf file, 
# makes new directory 'Bcftools' if it doesn't exists.
@timing
def bcftools(args):
	log_file=open('pipeline.log','a')
	if not os.path.exists('Bcftools'):
		os.makedirs('Bcftools')

	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*_sorted.bam'):
			cmd7 = ("bcftools mpileup --threads %s -Ou -f %s %s \
				| bcftools call --threads %s -Ou -mv \
	 			| bcftools filter -s LowQual -e 'QUAL<20' -Oz -o %s") \
			% (threads, add+args.ref, sorted_bam_out, threads, bcftools_out)
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
def snpEff_test():
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
				from parse_gff_2 import main as parse
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
def snpsift():
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

# Add empty file when the pipeline is done.
def done():
	open("pipeline.done", 'a').close()
	

# If the '-b' flag is used this function will run, 
# excluding the first steps of the program. 
def input_files():
	if args.sortbam:
		bcftools(args)
		annotation(args)
	else:
		bam_input(args)
		samtools_index()
		bcftools(args)
		annotation(args)
	if args.snpsift:
		snpsift()

	if args.clean:
		clean()

	if args.done:
		done()

# Exit program.
def exit():
	sys.exit()

def main():
	snpEff_test()

	if args.bamfile or args.sortbam:
		input_files()
		exit()	

	bowtie2(args)
	samtools_view()
	samtools_sort()
	samtools_index()
	bcftools(args)
	annotation(args)

	if args.snpsift:
		snpsift()
	
	if args.clean:
		clean()

	if args.done:
		done()

if __name__ == "__main__":
	# Use gff parser without running whole pipeline.
	if args.gff and args.feature:
		try:
			annotation(args)
		except:
			main()
	else:
		main()

