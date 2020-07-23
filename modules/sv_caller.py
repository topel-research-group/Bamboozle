#!/usr/bin/env python3
# coding=utf-8

#André Soares - 15/04/2020
#Pseudocode for a script calling SVs and inferring their LOF potential from BAM/SAM files

#	SVCaller provides a pipeline to call and annotate structural variants
#	from reads mapped to a reference genome.
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

#checks for the presence of indices for the FASTA reference
#if non existing, uses bwa to create new indices
def ref_check(reference):
	#list of bwa index files
	bwa_suf = (".amb",".ann",".bwt",".pac",".sa")

	#creating a list of bwa indices if they exist
	suf_list = []
	for file in os.listdir(os.path.dirname(reference)):
		if file.endswith(bwa_suf):
			print(file + " is already in the directory")
			suf_list.append(file)
	#if list of bwa indices in folder is <1, create indices
	if len(suf_list) < 1:
		#this doesn't overwrite the reference if called
		cmd3 = ['bwa','index', reference]
		proc_3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE,
			shell=False, universal_newlines=True)
		std_out, std_error = proc_3.communicate()
		print("bwa-mem indices didn't exist but they sure do now")

#
#GATK goes here
#

# BEDTOOLS masking of SV calls goes here

def masking(vcf_out, refpil, masked_vcf_out):
	cmd5 = ['bedtools', 'intersect', '-v','-b', refpil, '-a', vcf_out, '-sorted', '-header']
	with open(masked_vcf_out, "w+") as f:
		proc_5 = subprocess.Popen(cmd5, stdout=f, shell=False)
	std_error = proc_5.communicate()

#
#R script was removed from here
#

# Checks for dependencies required for snpEff.
def snpeff(snpeffdb1, masked_ann_vcf_out, bamboozledir1, masked_vcf_out_lof_csv, masked_vcf_out_lof_ann):
	cmd7 = ['snpEff', 'eff', snpeffdb1.replace("'", ""), masked_ann_vcf_out, '-c', bamboozledir1+'/data/snpeff/snpEff.config', '-csvStats', masked_vcf_out_lof_csv]
	with open(masked_vcf_out_lof_ann, "w+") as f:
		proc_7 = subprocess.Popen(cmd7, stdout=f, shell=False)
	std_error = proc_7.communicate()

# Filters SnpEff (and GRIDSS) annotations and tidies headers
def filter(masked_vcf_out_lof_ann, masked_vcf_out_lof_ann_filt, masked_vcf_out_lof_ann_filt_clean):
	#removes FORMAT, INFO fields
	cmd8 = ['bcftools', 'annotate', '-x', 'FORMAT,INFO', masked_vcf_out_lof_ann, '-Oz', '-o', masked_vcf_out_lof_ann_filt+'.gz']
	proc_8 = subprocess.Popen(cmd8, shell=False)
	std_out, std_error = proc_8.communicate()

	cmd8_5 = ['tabix', '-p', 'vcf', masked_vcf_out_lof_ann_filt+'.gz']
	proc_8_5 = subprocess.Popen(cmd8_5, shell=False)
	std_out, std_error = proc_8_5.communicate()

	#bgzips, indexes filt file
	cmd9 = ['bgzip', masked_vcf_out_lof_ann]
	proc_9 = subprocess.Popen(cmd9, shell=False)
	std_out, std_error = proc_9.communicate()

	cmd9_5 = ['tabix', '-p', 'vcf', masked_vcf_out_lof_ann+'.gz']
	proc_9_5 = subprocess.Popen(cmd9_5, shell=False)
	std_out, std_error = proc_9_5.communicate()

	#adds only relevant header columns from filt file
	cmd10 = ['bcftools', 'annotate', '-c', 'FORMAT/GT,INFO/EVENT,INFO/REF,INFO/RP,INFO/RPQ,INFO/SVLEN,INFO/SVTYPE,INFO/SIMPLE_TYPE,INFO/ANN,INFO/LOF,INFO/NMD', '-a', masked_vcf_out_lof_ann+'.gz', masked_vcf_out_lof_ann_filt+'.gz', '-Oz', '-o', masked_vcf_out_lof_ann_filt_clean+'.gz']
	proc_10 =  subprocess.Popen(cmd10, shell=False)
	std_out, std_error = proc_10.communicate()

def main(args, bam_name):
	#gridss java
	java_gridss="/usr/local/packages/gridss-2.8.3/gridss-2.8.3-gridss-jar-with-dependencies.jar"
	#outputs for gridss
	vcf_out = "sv_caller_output/%s_svcalls.vcf" % (bam_name)
	assembly_bam_out = "sv_caller_output/%s_assembly.vcf" % (bam_name)
	#outputs for bedtools
	masked_vcf_out = "sv_caller_output/%s_sorted_masked.vcf" % (bam_name)
	#output for R script
	masked_ann_vcf_out = "sv_caller_output/%s_sv_annotated.vcf" % (bam_name)
	#outputs for snpeff
	masked_vcf_out_lof_csv = "sv_caller_output/%s_sorted_masked_lof.csv" % (bam_name)
	masked_vcf_out_lof_ann = "sv_caller_output/%s_sorted_masked_lof.vcf" % (bam_name)
	#outputs for bcftools
	masked_vcf_out_lof_ann_filt = "sv_caller_output/%s_sorted_masked_lof_filt.vcf" % (bam_name)
	masked_vcf_out_lof_ann_filt_clean = "sv_caller_output/%s_sorted_masked_lof_filt_clean.vcf" % (bam_name)
	#clean database variable
	snpeff_db = str(args.snpeffdb).strip('[]')
	
	#calling functions for sv_caller
	ref_check(args.ref)
	gridss(args.sortbam, args.ref, args.threads, java_gridss, assembly_bam_out, vcf_out)
	#only apply masking() if it's been called
	if args.masking:
		masking(vcf_out, args.masking, masked_vcf_out)
		annotate(masked_vcf_out, bam_name, args.bamboozledir)
		snpeff(snpeff_db, masked_ann_vcf_out, args.bamboozledir,masked_vcf_out_lof_csv, masked_vcf_out_lof_ann)
		filter(masked_vcf_out_lof_ann, masked_vcf_out_lof_ann_filt, masked_vcf_out_lof_ann_filt_clean)
	else:
		annotate(vcf_out, bam_name, args.bamboozledir)
		snpeff(snpeff_db, masked_ann_vcf_out, args.bamboozledir, masked_vcf_out_lof_csv, masked_vcf_out_lof_ann)
		filter(masked_vcf_out_lof_ann, masked_vcf_out_lof_ann_filt, masked_vcf_out_lof_ann_filt_clean)
