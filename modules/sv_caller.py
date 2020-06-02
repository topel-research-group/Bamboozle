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

#First things first - taking in all arguments needed
#arguments to variables
#args = parser.parse_args()
#bam = args.bamfile
#ref = args.ref
#refgff = args.reference_gff
#refpil = args.masking
#threads = args.threads
#extracting sample name from input BAM
#bam_name = bamfile[:-4]


def ref_check(args):
	#list of bwa index files
	bwa_suf = (".amb",".ann",".bwt",".pac",".sa")

	#creating a list of bwa indices if they exist
	suf_list = []
	for file in os.listdir(os.path.dirname(args.ref)):
		if file.endswith(bwa_suf):
			print(file + " is already in the directory")
			suf_list.append(file)
	#if list of bwa indices in folder is <1, create indices
	if len(suf_list) < 1:
		#this doesn't overwrite the reference if called
		cmd3 = "bwa index %s" % (args.ref)
		proc_3 = subprocess.Popen(cmd3, shell=True, universal_newlines=True)
		std_out, std_error = proc_3.communicate()
		print("bwa-mem indices didn't exist but they sure do now")
#ref_check(args.ref)


#GRIDSS

def gridss(bamfile, ref, threads, java_gridss, assembly_bam_out, vcf_out):

	#create output folder if it doesn't exist
	if not os.path.exists('sv_caller_output'):
		os.makedirs('sv_caller_output')

	#naming outputs
	vcf_out = "sv_caller_output/%s_sorted.vcf" % (bam_name)
	assembly_bam_out = "sv_caller_output/%s_sorted_assembly.vcf" % (bam_name)

	java_gridss="/usr/local/packages/gridss-2.8.3/gridss-2.8.3-gridss-jar-with-dependencies.jar"

	log_file=open('sv_caller_run.log','a')

	cmd4 = "gridss.sh %s,-r %s, -a %s, -o %s, -t %s, -j %s" % (bamfile, args.ref, assembly_bam_out, vcf_out, threads, java_gridss)
	proc_4 = subprocess.Popen(cmd4, shell=True)
	std_out, std_error = proc_4.communicate()
	
	while proc_4.wait() is None:
		pass
	proc_4.stdout.close()
	log_file.close()

	print("GRIDSS finished calling SVs")

#gridss(bamfile, args.ref, threads, java_gridss, assembly_bam_out, vcf_out)
#exit()
#
# HERE GOES R SCRIPT TO ANNOTATE DEL, INS, ETC
#

# BEDTOOLS masking of SV calls goes here
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

def main(args, bam_name):
	ref_check(args.ref)
	gridss(args.bamfile, args.ref, args.threads, java_gridss, assembly_bam_out, vcf_out)

#main(args, bam_name)
