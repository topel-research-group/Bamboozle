#!/usr/bin/env python3

#	Pipeline that performs bioinformatic analysis including SNP calling and effect prediction of fastq files or BAM file 
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
import subprocess
import argparse
import fnmatch
import os 

##################################################################################
parser = argparse.ArgumentParser(prog="ADD-SCRIPT-NAME-HERE")
parser.add_argument("-f", "--ref", required=True, help="Reference")
parser.add_argument("-F", "--forward", nargs='*', help="Forward reads")
parser.add_argument("-R", "--reverse", nargs='*', help="Reverse reads")
parser.add_argument("-b", "--bamfile", help="BAM infile")  
parser.add_argument("-t", "--threads", default=1, help="Threads")
parser.add_argument("-s", "--snpsift", action="store_true", help="Run snpSift")
parser.add_argument("-r", "--clean", action="store_true", help="Removes the SAM and BAM files")
parser.add_argument("-p", "--done", action="store_true", help="Add an empty file to mark the directory as done")
args = parser.parse_args()
##################################################################################

current_directory = os.getcwd()
name = os.path.basename(current_directory)
threads = str(args.threads)
base = name + '.contigs'
sam = name + '.sam'
bam = name + '.bam' 
sorted_bam_out = name + '_sorted.bam'
sorted_bam_bai = name + '_sorted.bam.bai'
bcftools_out = name + '.bcftools_filtered.vcf.gz'
annotated_vcf = name + '.snpeff_annotated.vcf'
annotated_vcf_gz = name + '.snpeff_annotated.vcf.gz'
annotated_table = name + '.snpsift_table.txt'
add = '../'
add2 = '../Bowtie2/'

##################################################################################
# Makes new directory 'Bowtie2' if it doesn't exists
bowtie2_directory = os.path.join(current_directory, r'Bowtie2')
if not os.path.exists(bowtie2_directory):
	os.makedirs(bowtie2_directory)

# Running bowtie2-build to index reference genome and bowtie2 to align 
def bowtie2(args):
	# Selected input files using forward and reverse flags, the flags can take several input files
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

	# Bowtie2-build step
	cmd1 = ['bowtie2-build', add+args.ref, base]
	process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, cwd='Bowtie2')	
	while process1.wait() is None:
		pass

	# Bowtie2 align step
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*.rev.1.bt2'):
			cmd2 = ['bowtie2', '-p', threads, '--no-unal', '--very-sensitive', '-x', base, '-1', file1, '-2', file2, '-S', sam]	
			process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, cwd='Bowtie2')
			while process2.wait() is None:
				pass

# Converting SAM to BAM using samtools view
def samtools_view(args):
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*.sam'):
			cmd3 = ('samtools view -@ %s -Sb %s > %s') % (args.threads, sam, bam)
			process3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE, cwd='Bowtie2', shell=True)
			while process3.wait() is None:
				pass

# Sort BAM files
def samtools_sort(args):
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*.bam'):
			cmd4 = ['samtools', 'sort', '-@', '$NSLOTS', bam, '-o', sorted_bam_out]
			process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, cwd='Bowtie2')
			while process4.wait() is None:
				pass	
		
# BAM input file by using the '-b' flag
def bam_input(args):
	cmd5 = ['samtools', 'sort', '-@', '$NSLOTS', add+args.bamfile, '-o', sorted_bam_out]
	process5 = subprocess.Popen(cmd5, stdout=subprocess.PIPE, cwd='Bowtie2')
	while process5.wait() is None:
		pass	

# Index sorted BAM files
def samtools_index(args):
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*_sorted.bam'):
			cmd7 = ['samtools','index', sorted_bam_out, sorted_bam_bai]
			process7 = subprocess.Popen(cmd7, stdout=subprocess.PIPE, cwd='Bowtie2')
			while process7.wait() is None:
				pass

# Remove SAM and BAM files
def clean():
	if args.clean:
		for samfile in os.listdir('Bowtie2'):
			if fnmatch.fnmatch(samfile, '*.sam'):
				os.remove(current_directory + '/Bowtie2/' + samfile)

		for bamfile in os.listdir('Bowtie2'):
			if fnmatch.fnmatch(bamfile, name + '.bam'):
				os.remove(current_directory + '/Bowtie2/' + bamfile)

# Variant calling using bcftools mpileup, makes new directory 'Bcftools' if it doesn't exists
def bcftools(args):
	bcftools_directory = os.path.join(current_directory, r'Bcftools')
	if not os.path.exists(bcftools_directory):
		os.makedirs(bcftools_directory)

	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*_sorted.bam'):
			cmd9 = ("bcftools mpileup -Ou -f %s %s | bcftools call -Ou -mv | bcftools filter -s LowQual \
			-e 'QUAL<20 || DP>100' -Oz -o %s") % (add+args.ref, add2+sorted_bam_out, bcftools_out)
			process9 = subprocess.Popen(cmd9, stdout=subprocess.PIPE, shell=True, cwd='Bcftools')
			while process9.wait() is None:
				pass

# Variant annotation and effect prediction, the snpEff_test checks if there is a Skeletonema database, 
# if it doesn't exists the program will exit and you have to create it by using 'snpEff build'
def snpEff_test():
	try:
		cmdx = ('snpEff databases | grep "Skeletonema"')
		processx = subprocess.check_output(cmdx, shell=True)
	
	except subprocess.CalledProcessError as e:
		if e.returncode >= 1:
			print('snpEff: Skeletonema database not found, exit program...')
			exit()

# Annotating variant calling output using snpEff, output is a vcf, the vcf file is bgzipped to work as an input file to the Fst analysis,
# the original vcf file is kept by using the -c flag 
def annotation(args):					
	for file in os.listdir('Bcftools'):
		if fnmatch.fnmatch(file, '*.bcftools_filtered.vcf.gz'):
			cmd11 = ("snpEff -no-downstream -no-upstream -no-intron -no-intergenic -classic Skeletonema_marinoi_v1.1.1.1 \
			-stats snpEff_summary.html %s > %s") % (bcftools_out, annotated_vcf)
			process11 = subprocess.Popen(cmd11, stdout=subprocess.PIPE, shell=True, cwd='Bcftools')
			while process11.wait() is None:
				pass

			cmd12 = ('bgzip -c %s > %s') % (annotated_vcf, annotated_vcf_gz)
			process12 = subprocess.Popen(cmd12, stdout=subprocess.PIPE, shell=True, cwd='Bcftools')
			while process12.wait() is None:
				pass

# Filtering and making a summary of annotated files using the vcf (not bgzipped) output file from snpEff, the summary will be in table format
def snpsift(args):
	for file in os.listdir('Bcftools'):
		if fnmatch.fnmatch(file, '*_annotated.vcf'):
			cmd13 = ('java -jar /usr/local/packages/snpEff/SnpSift.jar extractFields -e "." -s "," %s \
			CHROM POS "EFF[*].GENE" REF ALT QUAL DP AF "EFF[*].EFFECT" "EFF[*].AA" "EFF[*].FUNCLASS" > %s') % (annotated_vcf, annotated_table) 
			process13 = subprocess.Popen(cmd13, stdout=subprocess.PIPE, shell=True, cwd='Bcftools')
			while process13.wait() is None:
				pass

# Add empty file when the pipeline is done
def done():
	open("pipeline.done", 'a').close()
	

# If the '-b' flag is used this function will run, excluding the first steps of the program 
def input_files():
	bam_input(args)
	samtools_index(args)
	bcftools(args)
	annotation(args)
	if args.snpsift:
		snpsift(args)

	if args.clean:
		clean()

	if args.done:
		done()

# Exit program
def exit():
	sys.exit()

def main():
	snpEff_test()

	if args.bamfile:
		input_files()
		exit()	

	bowtie2(args)
	samtools_view(args)
	samtools_sort(args)
	samtools_index(args)
	bcftools(args)
	annotation(args)

	if args.snpsift:
		snpsift(args)
	
	if args.clean:
		clean()

	if args.done:
		done()

if __name__ == "__main__":
	main()

