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
import subprocess
import argparse
import fnmatch
import os 

#######################################################################
parser = argparse.ArgumentParser(prog="ADD-SCRIPT-NAME-HERE")
parser.add_argument("-f", "--ref", required=True, help="Reference")
parser.add_argument("-F", "--forward", nargs='*', help="Forward reads")
parser.add_argument("-R", "--reverse", nargs='*', help="Reverse reads")
parser.add_argument("-b", "--bamfile", help="BAM infile")  
parser.add_argument("--gff", help="gff infile")  
parser.add_argument("--feature", help="Feature for gff parser")  
parser.add_argument("-t", "--threads", default=1, help="Threads")
parser.add_argument("-s", "--snpsift", action="store_true", \
			help="Run snpSift")
parser.add_argument("-r", "--clean", action="store_true", \
			help="Removes the SAM and BAM files")
parser.add_argument("-p", "--done", action="store_true", \
			help="Add an empty file to mark the directory as done")
args = parser.parse_args()

if args.feature and args.gff is None:
	parser.error("--feature requires --gff")
#######################################################################

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
annotated_table = name + '.snpsift_table.txt'
add = '../'
add2 = '../Bowtie2/'

#######################################################################
# Makes new directory 'Bowtie2' if it doesn't exists.
if not os.path.exists('Bowtie2'):
	os.makedirs('Bowtie2')

# Running bowtie2-build to index reference genome and bowtie2 to align. 
def bowtie2(args):
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
	process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, cwd='Bowtie2')	
	while process1.wait() is None:
		pass
	process1.stdout.close()

	# Bowtie2 align step, input are the index files from Bowtie2-build, 
	# fastq files (forward and reverse) the output is a SAM file.
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*.rev.1.bt2'):
			cmd2 = ['bowtie2', '-p', threads, '--no-unal', '--very-sensitive', \
			'-x', base, '-1', file1, '-2', file2, '-S', sam]	
			process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, cwd='Bowtie2')
			while process2.wait() is None:
				pass
			process2.stdout.close()

# Converting SAM to BAM using samtools view.
def samtools_view():
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*.sam'):
			cmd3 = ('samtools view -@ %s -Sb %s > %s') % (args.threads, sam, bam)
			process3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE, cwd='Bowtie2', shell=True)
			while process3.wait() is None:
				pass
			process3.stdout.close()

# Sort BAM files.
def samtools_sort():
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*.bam'):
			cmd4 = ['samtools', 'sort', '-@', '$NSLOTS', bam, '-o', sorted_bam_out]
			process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, cwd='Bowtie2')
			while process4.wait() is None:
				pass	
			process4.stdout.close()
		
# BAM input file by using the '-b' flag.
def bam_input(args):
	cmd5 = ['samtools', 'sort', '-@', '$NSLOTS', add+args.bamfile, '-o', sorted_bam_out]
	process5 = subprocess.Popen(cmd5, stdout=subprocess.PIPE, cwd='Bowtie2')
	while process5.wait() is None:
		pass	
	process5.stdout.close()

# Index sorted BAM files.
def samtools_index():
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*_sorted.bam'):
			cmd6 = ['samtools','index', sorted_bam_out, sorted_bam_bai]
			process6 = subprocess.Popen(cmd6, stdout=subprocess.PIPE, cwd='Bowtie2')
			while process6.wait() is None:
				pass
			process6.stdout.close()

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
def bcftools(args):
	if not os.path.exists('Bcftools'):
		os.makedirs('Bcftools')

	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*_sorted.bam'):
			cmd7 = ("bcftools mpileup --threads %s -Ou -f %s %s | bcftools call -Ou -mv \
	 		| bcftools filter -s LowQual -e 'QUAL<20 || DP>100' -Oz -o %s") \
			% (threads, add+args.ref, add2+sorted_bam_out, bcftools_out)
			process7 = subprocess.Popen(cmd7, stdout=subprocess.PIPE, shell=True, cwd='Bcftools')
			while process7.wait() is None:
				pass
			process7.stdout.close()

# Checks for dependencies required for snpEff 
def snpEff_test():
	# Checks if there is a Skeletonema database, if it doesn't exists the program 
	# will exit and it has to be created using 'snpEff build'.
	try:
		cmdx = ('snpEff databases | grep "Skeletonema"')
		processx = subprocess.check_output(cmdx, shell=True)
	
	except subprocess.CalledProcessError as e:
		if e.returncode >= 1:
			print('snpEff: Skeletonema database not found, exit program...')
			exit()

	# Try to import gffutils if gff and feature flag is used	
	if args.gff and args.feature:
		try:
			import gffutils

		except ImportError:
			sys.stderr.write("[Error] The python module \"gffutils\" is not installed\n") 
			sys.stderr.write("[--] Would you like to install it now using 'pip install gffutils' [Y/N]?\n")
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
def annotation():					
	for file in os.listdir('Bcftools'):
		if fnmatch.fnmatch(file, '*.bcftools_filtered.vcf.gz'):
			my_output = annotated_vcf
			my_interval = ""
			if args.gff and args.feature:
				from modules.parse_gff import main as parse
				parse(args.gff, args.feature)
				out = add+'out.gff'	
				my_interval = "-interval %s" % out
				my_output = name + '_' + args.feature + '.snpeff_annotated.vcf'
			my_args = my_interval + " -no-downstream -no-upstream -no-intron -no-intergenic \
			-classic Skeletonema_marinoi_v1.1.1.1 -stats snpEff_summary.html"
			cmd8 = ("snpEff	%s %s > %s") % (my_args, bcftools_out, my_output)
			process8 = subprocess.Popen(cmd8, stdout=subprocess.PIPE, shell=True, cwd='Bcftools')
			while process8.wait() is None:
				pass
			process8.stdout.close()

			cmd9 = ('bgzip -c %s > %s') % (my_output, my_output +'.gz')
			process9 = subprocess.Popen(cmd9, stdout=subprocess.PIPE, shell=True, cwd='Bcftools')
			while process9.wait() is None:
				pass
			process9.stdout.close()

# Filtering and making a summary of annotated files using 
# the vcf (not bgzipped) output file from snpEff, 
# the summary will be in table format, tab separated.
def snpsift():
	for file in os.listdir('Bcftools'):
		if fnmatch.fnmatch(file, '*_annotated.vcf'):
			cmd10 = ('java -jar /usr/local/packages/snpEff/SnpSift.jar \
			extractFields -e "." -s "," %s CHROM POS "EFF[*].GENE" REF ALT QUAL DP AF \
			"EFF[*].EFFECT" "EFF[*].AA" "EFF[*].FUNCLASS" > %s') \
			% (annotated_vcf, annotated_table) 
			process10 = subprocess.Popen(cmd10, stdout=subprocess.PIPE, shell=True, cwd='Bcftools')
			while process10.wait() is None:
				pass
			process10.stdout.close()

# Add empty file when the pipeline is done.
def done():
	open("pipeline.done", 'a').close()
	

# If the '-b' flag is used this function will run, 
# excluding the first steps of the program. 
def input_files():
	bam_input(args)
	samtools_index()
	bcftools(args)
	annotation()
	if args.snpsift:
		snpsift()

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
	samtools_view()
	samtools_sort()
	samtools_index()
	bcftools(args)
	annotation()

	if args.snpsift:
		snpsift()
	
	if args.clean:
		clean()

	if args.done:
		done()

if __name__ == "__main__":
	main()

