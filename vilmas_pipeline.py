#!/usr/bin/env python

import sys
import subprocess
import argparse
import fnmatch
import os 
													        
##################################################################################
parser = argparse.ArgumentParser(prog="ADD-SCRIPT-NAME-HERE")
parser.add_argument("-v", "--verbose", action="store_true", help="Be more verbose")
parser.add_argument("-b", "--bowtie2", action="store_true", help="Run Bowtie2")
parser.add_argument("-s", "--samtools", action="store_true", help="Run Samtools")
parser.add_argument("-c", "--vcalling", action="store_true", help="Run mpileup")
parser.add_argument("-a", "--annotation", action="store_true", help="Run snpEff")
parser.add_argument("-f", "--filtering", action="store_true", help="Run snpSift")
args = parser.parse_args()
##################################################################################

ref = 'Skeletonema_marinoi_Ref_v1.1.1.fst'
base = 'P8352_102.contigs'
file1 = 'P8352_102_S1_L001_R1_001.fastq.gz'
file2 = 'P8352_102_S1_L001_R2_001.fastq.gz'
sam = 'P8352_102.sam'
bam = 'P8352_102.bam'
sorted_bam = 'P8352_102_sorted.bam'
sorted_bam_bai = 'P8352_102_sorted.bam.bai'
mpileup_out = 'outputvcffile.mpileup.vcf.gz'
annotated_vcf = 'snpeff_annotated.vcf'
annotated_filtered_vcf = 'snpeff_annotated_filtered.vcf'

##################################################################################

# Running bowtie2-build to index reference genome and bowtie2 to align
def bowtie2():
	cmd1 = ['bowtie2-build', ref, base]
	process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
	for line in process1.stdout:
		print line
	
	for file in os.listdir('.'):
		if fnmatch.fnmatch(file, '*.rev.2.bt2'):
			cmd2 = ['bowtie2', '--no-unal', '--very-sensitive', '-x', base, '-1', file1, '-2', file2, '-S', sam]	
			process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE)

# Converting SAM to BAM using samtools view
def samtools():
	for file in os.listdir('.'):
		if fnmatch.fnmatch(file, '*.sam'):
			with open(bam, 'wb',0) as output_file: 
				cmd3 = ['samtools', 'view', '-Sb', sam]
				process3 = subprocess.Popen(cmd3, stdout=output_file)
				while process3.wait() is None:
					pass
	# Sort BAM files
	for file in os.listdir('.'):
       		if fnmatch.fnmatch(file, '*.bam'):
			cmd4 = ['samtools', 'sort', bam, '-o', sorted_bam]
			process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE)
			while process4.wait() is None:
                        	pass	
	# Index sorted BAM files
	for file in os.listdir('.'):
		if fnmatch.fnmatch(file, '*_sorted.bam'):
			cmd5 = ['samtools','index', sorted_bam, sorted_bam_bai]
			process5 = subprocess.Popen(cmd5, stdout=subprocess.PIPE)

# Variant calling using mpileup
def vcalling():
	for file in os.listdir('.'):
                if fnmatch.fnmatch(file, '*_sorted.bam'):
			cmd6 = ("samtools mpileup -u -g -f Skeletonema_marinoi_Ref_v1.1.1.fst %s | bcftools call -v -m -O z -o %s") \
			% (sorted_bam, mpileup_out)
			process6 = subprocess.Popen(cmd6, stdout=subprocess.PIPE, shell=True)
			while process6.wait() is None:
                                pass

# Variant annotation and effect prediction
def annotation():
	for file in os.listdir('.'):
                if fnmatch.fnmatch(file, '*.mpileup.vcf.gz'):
			cmd7 = ("snpEff -no-downstream -no-upstream -no-intron -no-intergenic -classic Skeletonema_marinoi_v1.1.1.1 \
			-stats snpEff_summary.html %s > %s") % (mpileup_out, annotated_vcf)
			process7 = subprocess.Popen(cmd7, stdout=subprocess.PIPE, shell=True)
			while process7.wait() is None:
                                pass

# Filtering and manipulation of annotated files
def filtering():
	for file in os.listdir('.'):
                if fnmatch.fnmatch(file, '*_annotated.vcf'):
			cmd8 = ('cat %s | java -jar /usr/local/packages/snpEff/SnpSift.jar \
			filter "(QUAL>=10)&(DP>=10)" > %s') % (annotated_vcf, annotated_filtered_vcf) 
			process8 = subprocess.Popen(cmd8, stdout=subprocess.PIPE, shell=True)
			while process8.wait() is None:
				pass

#def main():

if __name__ == "__main__":

	if args.bowtie2:
		bowtie2()

	if args.samtools:
		samtools()

	if args.vcalling:
		vcalling()

	if args.annotation:
		annotation()

	if args.filtering:
		filtering()


