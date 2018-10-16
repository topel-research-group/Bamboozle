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
parser.add_argument("-l", "--bcftools", action="store_true", help="Run Bcftools")
parser.add_argument("-t", "--threads", default=1, help="Threads")
parser.add_argument("-r", "--clean", action="store_true", help="Removes the SAM and BAM files")
parser.add_argument("-p", "--done", action="store_true", help="Add an empty file to mark the directory as done")
args = parser.parse_args()
##################################################################################

current_directory = os.getcwd()
name = os.path.basename(current_directory)
ref = '/proj/data11/vilma/Pipeline_vilma/Skeletonema_marinoi_Ref_v1.1.1.fst'
base = name + '.contigs'
sam = name + '.sam'
bam = name + '.bam'
sorted_bam = name + '_sorted.bam'
sorted_bam_input = current_directory + '/Bowtie2/' + name + '_sorted.bam'
sorted_bam_input_2 = current_directory + '/' + name + '_sorted.bam'
sorted_bam_bai = name + '_sorted.bam.bai'
bcftools_out = name + '.bcftools_filtered.vcf.gz'
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
	process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, cwd='Bowtie2')	
	while process1.wait() is None:
		pass

	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*.rev.1.bt2'):
			cmd2 = ['bowtie2', '-p', args.threads, '--no-unal', '--very-sensitive', '-x', base, '-1', file1, '-2', file2, '-S', sam]	
			process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, cwd='Bowtie2')
			while process2.wait() is None:
                                pass

	# Converting SAM to BAM using samtools view
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*.sam'):
			cmd3 = ('samtools view -@ %s -Sb %s > %s') % (args.threads, sam, bam)
			process3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE, cwd='Bowtie2', shell=True)
			while process3.wait() is None:
				pass
		else:
			continue	

# Sort BAM files
def samtools():
	for file in os.listdir('Bowtie2'):
       		if fnmatch.fnmatch(file, '*.bam'):
			cmd4 = ['samtools', 'sort', '-@', '$NSLOTS', bam, '-o', sorted_bam]
			process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, cwd='Bowtie2')
			while process4.wait() is None:
                    		pass	
		else: 
			continue	

	# If the input files are BAM files and previous steps haven't been made the program will look 
	# for the BAM files in cwd instead of the Bowtie2 directory
	for file in os.listdir('.'):
		if fnmatch.fnmatch(file, '*.bam'):
			input_bam = file
            		cmd4 = ['samtools', 'sort', '-@', '$NSLOTS', input_bam, '-o', sorted_bam] 
              		process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, cwd='.')
               		while process4.wait() is None:
               			pass  

		# Index sorted BAM files
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*_sorted.bam'):
			cmd5 = ['samtools','index', sorted_bam, sorted_bam_bai]
			process5 = subprocess.Popen(cmd5, stdout=subprocess.PIPE, cwd='Bowtie2')
			while process5.wait() is None:
				pass
		else:
			continue	

	for file in os.listdir('.'):
		if fnmatch.fnmatch(file, '*_sorted.bam'):
                    	cmd5 = ['samtools','index', sorted_bam, sorted_bam_bai]
                      	process5 = subprocess.Popen(cmd5, stdout=subprocess.PIPE, cwd=('.'))

        # Remove SAM and BAM files
	if args.clean:
		for samfile in os.listdir('Bowtie2'):
			if fnmatch.fnmatch(samfile, '*.sam'):
				os.remove(current_directory + '/Bowtie2/' + samfile)

		for bamfile in os.listdir('Bowtie2'):
			if fnmatch.fnmatch(bamfile, name + '.bam'):
				os.remove(current_directory + '/Bowtie2/' + bamfile)

# Variant calling using bcftools mpileup
def bcftools():
	bcftools_directory = os.path.join(current_directory, r'Bcftools')
        if not os.path.exists(bcftools_directory):
                os.makedirs(bcftools_directory)

        for file in os.listdir('Bowtie2'):
                if fnmatch.fnmatch(file, '*_sorted.bam'):
                        cmd6 = ("bcftools mpileup -Ou -f %s %s | bcftools call -Ou -mv | bcftools filter -s LowQual \
			-e 'QUAL<20 || DP>100' -Oz -o %s") % (ref, sorted_bam_input, bcftools_out)
                        process6 = subprocess.Popen(cmd6, stdout=subprocess.PIPE, shell=True, cwd='Bcftools')
                        while process6.wait() is None:
                                pass
		else:
			continue	

	# If input files were BAM files
	for file in os.listdir('.'):
                if fnmatch.fnmatch(file, '*_sorted.bam'):
                        cmd6 = ("bcftools mpileup -Ou -f %s %s | bcftools call -Ou -mv | bcftools filter -s LowQual \
                        -e 'QUAL<20 || DP>100' -Oz -o %s") % (ref, sorted_bam_input_2, bcftools_out, bcftools_out)
                        process6 = subprocess.Popen(cmd6, stdout=subprocess.PIPE, shell=True, cwd='Bcftools')
                        while process6.wait() is None:
                                pass

# Variant annotation and effect prediction
def annotation():
	for file in os.listdir('Bcftools'):
                if fnmatch.fnmatch(file, '*.bcftools_filtered.vcf.gz'):
			cmd7 = ("snpEff -no-downstream -no-upstream -no-intron -no-intergenic -classic Skeletonema_marinoi_v1.1.1.1 \
			-stats snpEff_summary.html %s > %s") % (bcftools_out, annotated_vcf)
			process7 = subprocess.Popen(cmd7, stdout=subprocess.PIPE, shell=True, cwd='Bcftools')
			while process7.wait() is None:
                                pass

# Filtering and manipulation of annotated files
def filtering():
	for file in os.listdir('Bcftools'):
                if fnmatch.fnmatch(file, '*_annotated.vcf'):
			cmd8 = ('cat %s | java -jar /usr/local/packages/snpEff/SnpSift.jar \
			filter "(QUAL>=10)&(DP>=10)" > %s') % (annotated_vcf, annotated_filtered_vcf) 
			process8 = subprocess.Popen(cmd8, stdout=subprocess.PIPE, shell=True, cwd='Bcftools')
			while process8.wait() is None:
				pass

# Add empty file when the pipeline is done
def done():
	open("pipeline.done", 'a').close()
	
def main():
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

	if args.done:
		done()

if __name__ == "__main__":
	main()

