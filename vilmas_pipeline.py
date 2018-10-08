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
parser.add_argument("-t", "--bcftools", action="store_true", help="Run Bcftools")
args = parser.parse_args()
##################################################################################

current_directory = os.getcwd()
name = os.path.basename(current_directory)
ref = '/proj/data11/vilma/Pipeline_vilma/P8352_102/Skeletonema_marinoi_Ref_v1.1.1.fst'
base = name + '.contigs'
#file1 = '/proj/data11/vilma/Pipeline_vilma/P8352_102/P8352_102_S1_L001_R1_001.fastq.gz'
#file2 = '/proj/data11/vilma/Pipeline_vilma/P8352_102/P8352_102_S1_L001_R2_001.fastq.gz'
sam = name + '.sam'
bam = name + '.bam'
sorted_bam = name + '_sorted.bam'
sorted_bam_input = current_directory + '/Bowtie2/' + name + '_sorted.bam'
sorted_bam_bai = name + '_sorted.bam.bai'
mpileup_out = name + '.mpileup.vcf.gz'
bcftools_out = name + '.bcftools_filtered.vcf'
annotated_vcf = name + '.snpeff_annotated.vcf'
annotated_filtered_vcf = name + '_annotated_filtered.vcf'

# Find the files in current working directory
file1 = current_directory + '/' 
for f1 in os.listdir('.'):
        if fnmatch.fnmatch(f1, '*_R1_*f*q.gz'):
                file1+=str(f1)

file2 = current_directory + '/'
for f2 in os.listdir('.'):
        if fnmatch.fnmatch(f2, '*_R2_*f*q.gz'):
                file2+=str(f2)

##################################################################################

# Running bowtie2-build to index reference genome and bowtie2 to align
def bowtie2():
	bowtie2_directory = os.path.join(current_directory, r'Bowtie2')
	if not os.path.exists(bowtie2_directory):
   		os.makedirs(bowtie2_directory)
	cmd1 = ['bowtie2-build', ref, base]
	process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, cwd=('Bowtie2'))	
	for line in process1.stdout:
		print line
	
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*.rev.2.bt2'):
			cmd2 = ['bowtie2', '--no-unal', '--very-sensitive', '-x', base, '-1', file1, '-2', file2, '-S', sam]	
			process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, cwd=('Bowtie2'))
			while process2.wait() is None:
                                pass

# Converting SAM to BAM using samtools view
def samtools():
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*.sam'):
			cmd3 = ('samtools view -Sb %s > %s') % (sam, bam)
			process3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE, cwd=('Bowtie2'), shell=True)
			while process3.wait() is None:
				pass

	# Sort BAM files
	for file in os.listdir('Bowtie2'):
       		if fnmatch.fnmatch(file, '*.bam'):
			cmd4 = ['samtools', 'sort', bam, '-o', sorted_bam]
			process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, cwd=('Bowtie2'))
			while process4.wait() is None:
                        	pass	

	# Index sorted BAM files
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*_sorted.bam'):
			cmd5 = ['samtools','index', sorted_bam, sorted_bam_bai]
			process5 = subprocess.Popen(cmd5, stdout=subprocess.PIPE, cwd=('Bowtie2'))

# Variant calling using samtools mpileup
def vcalling():
	vcalling_directory = os.path.join(current_directory, r'Vcalling')
        if not os.path.exists(vcalling_directory):
                os.makedirs(vcalling_directory)
	for file in os.listdir('Bowtie2'):
                if fnmatch.fnmatch(file, '*_sorted.bam'):
			cmd6 = ("samtools mpileup -u -g -f %s %s | bcftools call -v -m -O z -o %s") \
			% (ref, sorted_bam_input, mpileup_out)
			process6 = subprocess.Popen(cmd6, stdout=subprocess.PIPE, shell=True, cwd='Vcalling')
			while process6.wait() is None:
                                pass

# Variant calling using bcftools mpileup
def bcftools():
	bcftools_directory = os.path.join(current_directory, r'Bcftools')
        if not os.path.exists(bcftools_directory):
                os.makedirs(bcftools_directory)
        for file in os.listdir('Bowtie2'):
                if fnmatch.fnmatch(file, '*_sorted.bam'):
                        cmd6 = ("bcftools mpileup -Oz -f %s %s | bcftools call -v -m -O z | bcftools filter -s LowQual -e 'QUAL<20 || DP>100' > %s") \
                        % (ref, sorted_bam_input, bcftools_out)
                        process6 = subprocess.Popen(cmd6, stdout=subprocess.PIPE, shell=True, cwd='Bcftools')
                        while process6.wait() is None:
                                pass

# Variant annotation and effect prediction
def annotation():
	for file in os.listdir('Vcalling'):
                if fnmatch.fnmatch(file, '*.mpileup.vcf.gz'):
			cmd7 = ("snpEff -no-downstream -no-upstream -no-intron -no-intergenic -classic Skeletonema_marinoi_v1.1.1.1 \
			-stats snpEff_summary.html %s > %s") % (mpileup_out, annotated_vcf)
			process7 = subprocess.Popen(cmd7, stdout=subprocess.PIPE, shell=True, cwd='Vcalling')
			while process7.wait() is None:
                                pass

# Filtering and manipulation of annotated files
def filtering():
	for file in os.listdir('Vcalling'):
                if fnmatch.fnmatch(file, '*_annotated.vcf'):
			cmd8 = ('cat %s | java -jar /usr/local/packages/snpEff/SnpSift.jar \
			filter "(QUAL>=10)&(DP>=10)" > %s') % (annotated_vcf, annotated_filtered_vcf) 
			process8 = subprocess.Popen(cmd8, stdout=subprocess.PIPE, shell=True, cwd='Vcalling')
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

	if args.bcftools:
		bcftools()
