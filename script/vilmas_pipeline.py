#!/usr/bin/env python

import sys
import subprocess
import argparse
import fnmatch
import os 
													        
##################################################################################
parser = argparse.ArgumentParser(prog="ADD-SCRIPT-NAME-HERE")
#parser.add_argument("-f", "--ref", required=True, help="Reference")
parser.add_argument("-s", "--snpsift", action="store_true", help="Run snpSift")
#parser.add_argument("-i", "--infile", help="BAM infile")  
parser.add_argument("-t", "--threads", default=1, help="Threads")
parser.add_argument("-r", "--clean", action="store_true", help="Removes the SAM and BAM files")
parser.add_argument("-p", "--done", action="store_true", help="Add an empty file to mark the directory as done")
args = parser.parse_args()
##################################################################################

current_directory = os.getcwd()
name = os.path.basename(current_directory)
ref = '/proj/data11/Skeletonema_marinoi_Ref_v1.1_Primary.all.fst'
base = name + '.contigs'
sam = name + '.sam'
bam = current_directory + '/Bowtie2/' + name + '.bam'
sorted_bam = name + '_sorted.bam'
sorted_bam_out = current_directory + '/Bowtie2/' + name + '_sorted.bam'
sorted_bam_bai = name + '_sorted.bam.bai'
bcftools_out = name + '.bcftools_filtered.vcf.gz'
annotated_vcf = name + '.snpeff_annotated.vcf'
annotated_filtered_vcf = name + '.snpsift_filtered.vcf'

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

def input_files():

# Running bowtie2-build to index reference genome and bowtie2 to align
def bowtie2_build():
	bowtie2_directory = os.path.join(current_directory, r'Bowtie2')
	if not os.path.exists(bowtie2_directory):
   		os.makedirs(bowtie2_directory)
	cmd1 = ['bowtie2-build', ref, base]
	process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, cwd='Bowtie2')	
	while process1.wait() is None:
		pass
def bowtie2():
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*.rev.1.bt2'):
			cmd2 = ['bowtie2', '-p', args.threads, '--no-unal', '--very-sensitive', '-x', base, '-1', file1, '-2', file2, '-S', sam]	
			process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, cwd='Bowtie2')
			while process2.wait() is None:
                                pass

def samtools_view():
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
def samtools_sort():
	cmd4 = ['samtools', 'sort', '-@', '$NSLOTS', bam, '-o', sorted_bam_out]
	for file in os.listdir('Bowtie2'):
       		if fnmatch.fnmatch(file, '*.bam'):
			process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, cwd='Bowtie2')
			while process4.wait() is None:
                    		pass	

	# If the input files are BAM files and previous steps haven't been made the program will look 
	# for the BAM files in cwd instead of the Bowtie2 directory
			for file in os.listdir('.'):
				if fnmatch.fnmatch(file, '*.bam'):
					process6 = subprocess.Popen(cmd4, stdout=subprocess.PIPE)
					while process6.wait() is None:
						pass  
		
	# BAM infile
	if args.infile:
		cmd5 = ['samtools', 'sort', '-@', '$NSLOTS', args.infile, '-o', sorted_bam_out]
                process5 = subprocess.Popen(cmd5, stdout=subprocess.PIPE)
                while process5.wait() is None:
                	pass	

def samtools_index():
	# Index sorted BAM files
	cmd7 = ['samtools','index', sorted_bam_out, sorted_bam_bai]
	for file in os.listdir('Bowtie2'):
		if fnmatch.fnmatch(file, '*_sorted.bam'):
			process7 = subprocess.Popen(cmd7, stdout=subprocess.PIPE, cwd='Bowtie2')
			while process7.wait() is None:
				pass
		else:
			for file in os.listdir('.'):
				if fnmatch.fnmatch(file, '*_sorted.bam'):
               	   			process8 = subprocess.Popen(cmd7, stdout=subprocess.PIPE )
					while process8.wait() is None:
                   	       			pass

def clean():
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
                        cmd9 = ("bcftools mpileup -Ou -f %s %s | bcftools call -Ou -mv | bcftools filter -s LowQual \
			-e 'QUAL<20 || DP>100' -Oz -o %s") % (ref, sorted_bam_out, bcftools_out)
                        process9 = subprocess.Popen(cmd9, stdout=subprocess.PIPE, shell=True, cwd='Bcftools')
                        while process9.wait() is None:
                                pass

# Variant annotation and effect prediction
def annotation():
	for file in os.listdir('Bcftools'):
                if fnmatch.fnmatch(file, '*.bcftools_filtered.vcf.gz'):
			cmd11 = ("snpEff -no-downstream -no-upstream -no-intron -no-intergenic -classic Skeletonema_marinoi_v1.1.1.1 \
			-stats snpEff_summary.html %s > %s") % (bcftools_out, annotated_vcf)
			process11 = subprocess.Popen(cmd11, stdout=subprocess.PIPE, shell=True, cwd='Bcftools')
			while process11.wait() is None:
                                pass

# Filtering and manipulation of annotated files
def filtering():
	for file in os.listdir('Bcftools'):
                if fnmatch.fnmatch(file, '*_annotated.vcf'):
			cmd12 = ('cat %s | java -jar /usr/local/packages/snpEff/SnpSift.jar \
			filter "(QUAL>=10)&(DP>=10)" > %s') % (annotated_vcf, annotated_filtered_vcf) 
			process12 = subprocess.Popen(cmd12, stdout=subprocess.PIPE, shell=True, cwd='Bcftools')
			while process12.wait() is None:
				pass

# Add empty file when the pipeline is done
def done():
	open("pipeline.done", 'a').close()
	
def main():
	input_files()
	bowtie2_build()
	bowtie2()
	samtools_view()
	samtools_sort()
	samtools_index()
	bcftools()
	annotation()

	if args.filtering:
 		filtering()
	
	if args.clean:
		clean()

	if args.done:
		done()

if __name__ == "__main__":
	main()

