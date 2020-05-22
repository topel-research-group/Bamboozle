#!/usr/bin/env python3
# coding=utf-8

#André Soares - 15/04/2020
#Pseudocode for a script calling SVs and inferring their LOF potential from BAM/SAM files

#very much based on Matt's(?) pipeline.py in Bamboozle/modules

#Input and arguments:
#	1) sorted BAM alignments of diatom (S. marinoi and others in the future?) genomes against a reference
#	2) genome reference
#	3) Pilon-corrected genome reference
#	4) GFF gene model reference
#	?) option to run submodules individually?

#Outputs: 
#	1) VCF files of alterations in gene models, with LOF prediction.
#	2) summary tables of alterations that should include type (deletion, insertion, inversion, duplication, other genomic rearrangmements),
#size in bp, quality of call, LOF potential, gene product, genomic location, zygosity, ?...


#Steps:
#0. Check if reference has been indexed with bwa, if not, index it
#1. Take in sorted alignments, verify they exist, are sorted BAM files
#1a. Take in fasta, GFF, Pilon-corrected assembly references as arguments (check if BAMs were aligned to those refs?)
#2. Run GRIDSS
#3. Mask VCF files of SV calls with bedtools using Pilon assembly
#4. Run snpEff on masked SV calls to obtain predicted alterations in gene models
#4a. (Modify output headers?)
#5. Run snpEff snpSift to obtain summary tables

#Help messages:
#1. Use a bwa-mem-indexed reference genome
#2. Even though ref has been indexed with bwa-mem it is ok to have bowtie2-aligned reads since from the same genome
#3. Input alignments have to be sorted BAM
#4. 

#  NOT SURE about this license? Secondary though

#       Pipeline that performs bioinformatic analysis including SNP calling
#       and effect prediction of fastq files or BAM file.
#
#       Copyright (C) 2020 André Soares, andre.soares@bioenv.gu.se
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

# Time decorator.
# returns time elapsed during processing of each function
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

#First things first - taking in all arguments needed
parser = argparse.ArgumentParser(description='SV caller - title says it all dunnit')
parser.add_argument('-B', '--sorted_BAM', required=True,
                   help='Sorted BAM file (can be bowtie2- or bwa-mem-aligned). Assuming sample name from first two "_" delimited fields of input name!')
#parser.add_argument('-n', '--sample_name',
#		   help='Sample name')
parser.add_argument('-F', '--reference_FASTA', required=True,
		   help='Reference genome in FASTA format (needs to be bwa-mem-indexed, this will checked automatically though)')
parser.add_argument('-G', '--reference_GFF', required=True,
		   help='Reference GFF with gene models for the reference genome')
parser.add_argument('-P', '--ref_Pilon',
		   help="Reference Pilon-corrected assembly in FASTA format (doesn't need indexing)")
parser.add_argument('-t', '--threads', default='8', type=int,
		   help='Number of threads to run sv_caller.py with')
#arguments to variables
args = parser.parse_args()
bam = args.one
bam_name = bam.split('_')[:1]
reffa = args.three
refgff = args.four
refpil = args.five
threads = args.six

#create output folder if it doesn't exist
if not os.path.exists('sv_caller_output'):
                os.makedirs('sv_caller_output')

# - function to make sure input is as needed!
#       1 - input alignment is sorted BAM
@timing
def bam_check(bam):
	#command to check out first line of BAM header and look for "coordinate" (= sorted)
	cmd1 = "samtools view -H %s \
		| head -n1 \
		| cut -f3 \
		| cut -f2 -d$':'" % (bam)
	proc_1 = subprocess.Popen(cmd1,
		cwd='sv_caller_output')

	if proc1.stdout.read() == "coordinate":
		print("Input BAM is sorted")
	else:
		cmd2 = "samtools sort %s \
        	        -o %s" % (bam,bam_name+"_sorted.bam")

		proc_2 = subprocess.Popen(cmd2,
			cwd='sv_caller_output')
		print("Input BAM has been sorted")

#       2 - input reference genome has associated bwa-mem index
@timing
def ref_check(reffa):
	#assuming a path string was provided, remove file?
	for file in os.listdir(os.path.dirname(reffa)):
		if file.endswith(('.fai','.amb','.ann','.bwt','.pac','.sa')):
			print("bwa-mem indices exist")
		else:
			cmd3 = "bwa index %s" % (reffa)
			proc_3 = subprocess.Popen(cmd3,
				cwd='sv_caller_output')
			print("bwa-mem indices didn't exist but they sure do now")
#
# makes new directory 'gridss' if it doesn't exist.
@timing
def gridss(bam, reffa, threads, java_gridss, assembly_bam_out, vcf_out):
        log_file=open('pipeline.log','a')
	vcf_out = bam_name+"_sorted.vcf"
	assembly_bam_out = bam_name+"_assembly.bam"
	
        cmd4 = "gridss.sh \
		%s, \
		-r %s, \
		-a %s, \
		-o %s, \
		-t %s, \
		-j %s" % (bam, reffa, assembly_bam_out, vcf_out, threads, java_gridss)

        process4 = subprocess.Popen(cmd4, \
                stdout=subprocess.PIPE, \
                stderr = log_file, \
                shell=True, \
                cwd='sv_caller_output')
        while process4.wait() is None:
                pass
        process4.stdout.close()
        log_file.close()

#
# HERE GOES R SCRIPT TO ANNOTATE DEL, INS, ETC
#

# BEDTOOLS masking of SV calls goes here
@timing
def masking(vcf_out, refpil, masked_vcf_out):
	masked_vcf_out = bam_name+"_sorted_masked.vcf"

	cmd5 = "bedtools intersect \
		-v \
		-b %s \
		-a %s \
		-sorted \
		> %s" % (refpil, vcf_out, masked_vcf_out)
	proc_5 = subprocess.Popen(cmd5, \
		cwd='sv_caller_output')
# Checks for dependencies required for snpEff.
def snpeff(masked_vcf_out, masked_vcf_out_csv, masked_vcf_out_ann):
	masked_vcf_out_csv = bam_name+"_sorted_masked.csv"
	masked_vcf_out_ann = bam_name+"_sorted_masked_ann.vcf"
        # Checks if there is a Skeletonema database,
        # if it doesn't exists the program will exit
        # and it has to be created using 'snpEff build'.
        try:
                cmd6 = "snpEff databases | grep Smarinoi.v112"
		################## this didn't work last I tried, check snpeff database!
                proc_6 = subprocess.check_output(cmd6, shell=True)

        except subprocess.CalledProcessError as e:
                if e.returncode >= 1:
                        print('snpEff: Skeletonema database not found, exit program...')
                        exit()
	cmd7 = "snpEff eff Smarinoi.v112 \
		%s \
		-c /home/andre/snpEff.config \
		-csvStats %s \
		> %s" % (masked_vcf_out, masked_vcf_out_csv, masked_vcf_out_ann)
	proc_7 = subprocess.Popen(cmd7, \
                cwd='sv_caller_output')

def main():
	bam_check()
	ref_check()
	gridss()

main()
