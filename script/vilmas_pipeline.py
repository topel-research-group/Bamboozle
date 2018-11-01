#!/usr/bin/env python

import sys
import subprocess
import argparse
import fnmatch
import os 
													        
##################################################################################
parser = argparse.ArgumentParser(prog="ADD-SCRIPT-NAME-HERE")
parser.add_argument("-f", "--ref", required=True, help="Reference")
parser.add_argument("-s", "--snpsift", action="store_true", help="Run snpSift")
parser.add_argument("-c", "--cwd", action="store_true", help="Find the files in current working directory")
parser.add_argument("-F", "--forward", nargs='*', help="Forward reads")
parser.add_argument("-R", "--reverse", nargs='*', help="Reverse reads")
parser.add_argument("-b", "--bamfile", help="BAM infile")  
parser.add_argument("-t", "--threads", default=1, help="Threads")
parser.add_argument("-r", "--clean", action="store_true", help="Removes the SAM and BAM files")
parser.add_argument("-p", "--done", action="store_true", help="Add an empty file to mark the directory as done")
args = parser.parse_args()
##################################################################################

current_directory = os.getcwd()
name = os.path.basename(current_directory)
base = name + '.contigs'
sam = name + '.sam'
bam = current_directory + '/Bowtie2/' + name + '.bam'
sorted_bam = name + '_sorted.bam'
sorted_bam_out = current_directory + '/Bowtie2/' + name + '_sorted.bam'
sorted_bam_bai = name + '_sorted.bam.bai'
bcftools_out = name + '.bcftools_filtered.vcf.gz'
annotated_vcf = name + '.snpeff_annotated.vcf'
annotated_filtered_vcf = name + '.snpsift_filtered.vcf'

# Selected input files using forward and reverse flags, the flags can take several input files
f1 = [] 
if args.forward:
	for name in args.forward:
		f1.append(current_directory + '/' + name)
	file1 = ','.join(map(str, f1)) 

f2 = [] 
if args.reverse:
	for name2 in args.reverse:  
		f2.append(current_directory + '/' + name2)
	file2 = ','.join(map(str, f2))

# Find the files in current working directory
if args.cwd:
	for fname1 in os.listdir('.'):
		if fnmatch.fnmatch(fname1, '*_R1_*f*q.gz'):
			n1 = os.path.abspath(fname1)
			f1.append(n1)
	file1 = ','.join(map(str, f1))	

	for fname2 in os.listdir('.'):
		if fnmatch.fnmatch(fname2, '*_R2_*f*q.gz'):
			n2 = os.path.abspath(fname2)
			f2.append(n2) 
        file2 = ','.join(map(str, f2)) 

##################################################################################

def input_files():
	if args.bamfile:
                bam_input()
                samtools_index()
                bcftools()
		snpEff_test()
                annotation()
	exit()	

# Running bowtie2-build to index reference genome and bowtie2 to elign
def bowtie2():
	bowtie2_directory = os.path.join(current_directory, r'Bowtie2')
	if not os.path.exists(bowtie2_directory):
   		os.makedirs(bowtie2_directory)
	cmd1 = ['bowtie2-build', args.ref, base]
	process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, cwd='Bowtie2')	
	while process1.wait() is None:
		pass

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
		
def bam_input():
	# BAM infile
	if args.bamfile:
		cmd5 = ['samtools', 'sort', '-@', '$NSLOTS', args.bamfile, '-o', sorted_bam_out]
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
			-e 'QUAL<20 || DP>100' -Oz -o %s") % (args.ref, sorted_bam_out, bcftools_out)
                        process9 = subprocess.Popen(cmd9, stdout=subprocess.PIPE, shell=True, cwd='Bcftools')
                        while process9.wait() is None:
                                pass

# Variant annotation and effect prediction
def snpEff_test():
	try:
		cmdx = ('snpEff databases | grep "Skeletonema"')
		processx = subprocess.check_output(cmdx, shell=True)
	except OSError:
		raise RuntimeError('Skeletonema database not found; building new database...')
		try:
			cmdy = ['snpEff', 'build', '-gff3', '-v', 'Skeletonema_marinoi_v1.1.1.1'] 
			processy = subprocess.Popen(cmdy, stdout=subprocess.PIPE)
		except:
			print('Skeletonema database not found; printing manual...')

			cmdz = ['snpEff', 'build', '-h']
			processz = subprocess.Popen(cmdz, stdout=subprocess.PIPE)
			while processz.wait() is None:
			        pass
			print ('''
			Build Skeletonema database:     
				Step 1: 
				To 'snpEFF.contig' Add this:
				#---
				# Skeletonema marinoi, version 1.1.1.1
				#---
				Skeletonema_marinoi_v1.1.1.1.genome : Skeletonema_marinoi_v1.1.1.1
				
				Step 2:
				mkdir path/to/snpEff/data/Skeletonema_marinoi_v1.1.1.1
				cd path/to/snpEff/data/Skeletonema_marinoi_v1.1.1.1
				rsync -hav path/to/Skeletonema.gff.gz . 
				mv Skeletonema.gff.gz genes.gff.gz 
				
				Additional:
				GFF3 files can have sequence information either in the same file or in a separate fasta file. 
				In order to add sequence information in the GFF file, you can do this:

				cat Skeletonema.gff > genes.gff
				echo "###"  >> genes.gff
				echo "##FASTA"  >> genes.gff
				cat sequence.fa  >> genes.gff
			''')

def annotation():					
	for file in os.listdir('Bcftools'):
                if fnmatch.fnmatch(file, '*.bcftools_filtered.vcf.gz'):
			cmd11 = ("snpEff -no-downstream -no-upstream -no-intron -no-intergenic -classic Skeletonema_marinoi_v1.1.1.1 \
			-stats snpEff_summary.html %s > %s") % (bcftools_out, annotated_vcf)
			process11 = subprocess.Popen(cmd11, stdout=subprocess.PIPE, shell=True, cwd='Bcftools')
			while process11.wait() is None:
                                pass

# Filtering and manipulation of annotated files
def snpsift():
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
	bowtie2()
	samtools_view()
	samtools_sort()
	samtools_index()
	bcftools()
	snpEff_test()
	annotation()

	if args.snpsift:
 		snpsift()
	
	if args.clean:
		clean()

	if args.done:
		done()

if __name__ == "__main__":
	main()

