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
def ref_check(reference, ref_dict):
	#list of bwa index files
	ref_dict_suf = ".dict"

	#creating a list of bwa indices if they exist
	suf_list = []
	for file in os.listdir(os.path.dirname(reference)):
		if file.endswith(ref_dict_suf):
			print(file + " is already in the directory")
			suf_list.append(file)
	#if list of bwa indices in folder is <1, create indices
	if len(ref_dict_suf) < 1:
		#this doesn't overwrite the reference if called
		cmd3 = ['java','-jar',java_picard,'CreateSequenceDictionary',\
			'REFERENCE='+reference,'OUTPUT='+ref_dict]
		proc_3 = subprocess.Popen(cmd3,
			shell=False)
		std_out, std_error = proc_3.communicate()
		print("GATK refernce index didn't exist but they sure do now")

#GATK
def run_gatk(bam_in, bam_rg, bam_dup, bam_fm, vcf_out, bam_name, reference, java_picard):
	#format input BAMs to make GATK happy
	cmd4a = ['java','-jar',java_picard,'AddOrReplaceReadGroups',\
		'I='+bam_in,'O='+bam_rg,\
		'SORT_ORDER=coordinate','RGID=foo','RGLB=bar',\
		'RGPL=illumina','RGSM='+bam_name,\
		'RGPU=bc1','CREATE_INDEX=True']
	proc_4a = subprocess.Popen(cmd4a,
		shell=False)
	std_out, std_error = proc_4a.communicate()

	#mark duplicate reads in input BAMs
	cmd4b = ['java','-jar',java_picard,'MarkDuplicates',\
		'I='+bam_rg,'O='+bam_dup,\
		'M=sv_caller_output/marked_dup_metrics.txt']
	proc_4b = subprocess.Popen(cmd4b,
		shell=False)
	std_out, std_error = proc_4b.communicate()

	#"ensures that all mate-pair information is in sync between each read and its mate pair"
	cmd4c = ['java','-jar',java_picard,'FixMateInformation',\
ww		'I='+bam_dup,'O='+bam_fm,\
		'ADD_MATE_CIGAR=true']
	proc_4c = subprocess.Popen(cmd4c,
		shell=False)
	std_out, std_error = proc_4c.communicate()

	#index output bam
	cmd4d = ['samtools', 'index', bam_fm]
	proc_4d = subprocess.Popen(cmd4d,
		shell=False)
	std_out, std_error = proc_4d.communicate()

	#validate output bam
	cmd4e = ['ValidateSamFile', '-I', bam_fm, '-MODE', 'SUMMARY']
	proc_4d = subprocess.Popen(cmd4d,
		shell=False)
	std_out, std_error = proc_4d.communicate()

	#run haplotype caller
	cmd4f = ['gatk', '--java-options', '"-Xmx4g"', 'HaplotypeCaller',\
		'-G', 'StandardAnnotation',\
		'-G', 'AS_StandardAnnotation',\
		'-G', 'StandardHCAnnotation',\
		'-R', reference,\
		'-I', bam_fm,\
		'-O', vcf_out]

	proc_4f = subprocess.Popen(cmd4f,
		shell=False)
	std_out, std_error = proc_4f.communicate()

# BEDTOOLS masking of SV calls goes here

def masking(vcf_out, refpil, masked_vcf_out):
	cmd5 = ['bedtools', 'intersect', '-v','-b', \
		refpil, \
		'-a', vcf_out, \
		'-sorted', '-header']
	with open(masked_vcf_out, "w+") as f:
		proc_5 = subprocess.Popen(cmd5, stdout=f, shell=False)
	std_error = proc_5.communicate()

# Checks for dependencies required for snpEff.
def snpeff(snpeffdb1, masked_vcf_out, bamboozledir1, masked_vcf_out_lof_csv, masked_vcf_out_lof_ann):
	cmd7 = ['snpEff', 'eff', snpeffdb1.replace("'", ""),\
		 masked_ann_vcf_out, \
		'-c', bamboozledir1+'/data/snpeff/snpEff.config', \
		'-csvStats', masked_vcf_out_lof_csv]
	with open(masked_vcf_out_lof_ann, "w+") as f:
		proc_7 = subprocess.Popen(cmd7, stdout=f, shell=False)
	std_error = proc_7.communicate()

def main(args, bam_name):
	#picard java
	java_picard="/usr/local/packages/picard-tools-2.18.26/picard.jar"
	#ref_dict
	ref_dict = "sv_caller_output/"+reference.replace('.fasta','.dict')
	#outputs for gatk
	bam_rg = "sv_caller_output/%s_rdgrp.bam" % (bam_name)
	bam_dup = "sv_caller_output/%s_rdgrp_nodups.bam" % (bam_name)
	bam_fm = "sv_caller_output/%s_rdgrp_nodups_fixmate.bam" % (bam_name)
	vcf_out = "sv_caller_output/%s_svcalls.vcf" % (bam_name)
	#outputs for bedtools
	masked_vcf_out = "sv_caller_output/%s_sorted_masked.vcf" % (bam_name)
	#outputs for snpeff
	masked_vcf_out_lof_csv = "sv_caller_output/%s_sorted_masked_lof.csv" % (bam_name)
	masked_vcf_out_lof_ann = "sv_caller_output/%s_sorted_masked_lof.vcf" % (bam_name)
	#clean database variable
	snpeff_db = str(args.snpeffdb).strip('[]')
	
	#calling functions for sv_caller
	ref_check(args.ref, ref_dict)
	run_gatk(args.sortbam, bam_rg, bam_dup, bam_fm, vcf_out, bam_name, args.ref, java_picard)
	#only apply masking() if it's been called
	if args.masking:
		masking(vcf_out, args.masking, masked_vcf_out)
		annotate(masked_vcf_out, bam_name, args.bamboozledir)
		snpeff(snpeff_db, masked_ann_vcf_out, args.bamboozledir,masked_vcf_out_lof_csv, masked_vcf_out_lof_ann)
	else:
		annotate(vcf_out, bam_name, args.bamboozledir)
		snpeff(snpeff_db, masked_ann_vcf_out, args.bamboozledir, masked_vcf_out_lof_csv, masked_vcf_out_lof_ann)
