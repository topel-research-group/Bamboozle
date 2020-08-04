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
	from pathlib import Path
	#list of bwa index files
	java_picard="/usr/local/packages/picard-tools-2.18.26/picard.jar"
	#creating a list of bwa indices if they exist
	suf_list = []
	if Path(ref_dict).exists():
		print("Reference dictionary is already in the directory")
	else:
		#this doesn't overwrite the reference if called
		cmd3 = ['java','-jar',java_picard,'CreateSequenceDictionary',\
			'R='+reference,\
			'O='+ref_dict]
		proc_3 = subprocess.Popen(cmd3, \
			shell=False)
		std_out, std_error = proc_3.communicate()
		print("GATK reference index didn't exist but it sure does now")

#GATK
def run_gatk(bamfile, reference, java_picard, threads):
	#generate variables for the function
	bam_name = os.path.basename(bamfile[:-4])
	out_dir = "sv_caller_output/%s" % (bam_name)
	bam_in = bamfile
	bam_rg = "%s/%s_rdgrp.bam" % (out_dir, bam_name)
	bam_dup = "%s/%s_rdgrp_nodups.bam" % (out_dir, bam_name)
	bam_fm = "%s/%s_rdgrp_nodups_fixmate.bam" % (out_dir, bam_name)
	vcf_out = "%s/%s_svcalls.vcf" % (out_dir, bam_name)
	#format input BAMs to make GATK happy
	cmd4a = ['java','-jar',java_picard,'AddOrReplaceReadGroups',\
		'I=', bam_in, 'O=', bam_rg,\
		'SORT_ORDER=coordinate','RGID=foo','RGLB=bar',\
		'RGPL=illumina','RGSM='+bam_name,\
		'RGPU=bc1','CREATE_INDEX=True']
	proc_4a = subprocess.Popen(cmd4a, \
		shell=False)
	std_out, std_error = proc_4a.communicate()
	#mark duplicate reads in input BAMs
	out_metrics = "%s/marked_dup_metrics.txt" % (out_dir)
	cmd4b = ['java','-jar',java_picard,'MarkDuplicates',\
		'I='+bam_rg,'O=', bam_dup,\
		'M=', out_metrics]
	proc_4b = subprocess.Popen(cmd4b, \
		shell=False)
	std_out, std_error = proc_4b.communicate()
	#"ensures that all mate-pair information is in sync between each read and its mate pair"
	cmd4c = ['java','-jar',java_picard,'FixMateInformation',\
		'I=', bam_dup,'O=', bam_fm,\
		'ADD_MATE_CIGAR=true']
	proc_4c = subprocess.Popen(cmd4c, \
		shell=False)
	std_out, std_error = proc_4c.communicate()
	#index output bam
	cmd4d = ['samtools', 'index', bam_fm]
	proc_4d = subprocess.Popen(cmd4d, \
		shell=False)
	std_out, std_error = proc_4d.communicate()
	#validate output bam
	cmd4e = ['ValidateSamFile', '-I', bam_fm, '-MODE', 'SUMMARY']
	proc_4d = subprocess.Popen(cmd4d, \
		shell=False)
	std_out, std_error = proc_4d.communicate()
	#run haplotype caller
	cmd4f = ['gatk', \
#		'--java-options "-Xmx4g"', \
		'HaplotypeCaller',\
		'-G', 'StandardAnnotation',\
		'-G', 'StandardHCAnnotation',\
		'-R', reference,\
#		'--native-pair-hmm-threads', threads,\
		'-I', bam_fm,\
		'-O', vcf_out]
	proc_4f = subprocess.Popen(cmd4f, \
		shell=False)
	std_out, std_error = proc_4f.communicate()

#masks the output for each called vcf with previous calls using reads generated to correct the assembly

def masking(bamfile, refpil):
	#generating variables for function
	bam_name = os.path.basename(bamfile[:-4])
	out_dir = "sv_caller_output/%s" % (bam_name)
	vcf_out = "%s/%s_svcalls.vcf" % (out_dir, bam_name)
	masked_vcf_out = "%s/%s_sorted_masked.vcf" % (out_dir, bam_name)

	cmd5 = ['bedtools', 'intersect', '-v','-b', \
		refpil, \
		'-a', vcf_out, \
		'-sorted', '-header']
	with open(masked_vcf_out, "w+") as f:
		proc_5 = subprocess.Popen(cmd5, stdout=f, shell=False)
	std_error = proc_5.communicate()

# Checks for dependencies required for snpEff.
def snpeff(snpeffdb1, bamfile, bamboozledir1, threads):
	#generating variables for function
	bam_name = os.path.basename(bamfile[:-4])
	out_dir = "sv_caller_output/%s" % (bam_name)
	masked_vcf_out = "%s/%s_sorted_masked.vcf" % (out_dir, bam_name)
	masked_vcf_out_lof_csv = "%s/%s_sorted_masked_lof.csv" % (out_dir, bam_name)
	masked_vcf_out_lof_ann = "%s/%s_sorted_masked_lof.vcf" % (out_dir, bam_name)

	cmd7 = ['SnpEff', snpeffdb1.replace("'", ""), masked_vcf_out, '-t', threads, '-c', bamboozledir1+'/data/snpeff/snpEff.config', '-lof', '-noStats']
	#	'-csvStats', masked_vcf_out_lof_csv]
	with open(masked_vcf_out_lof_ann, "w+") as f:
		proc_7 = subprocess.Popen(cmd7, stdout=f, shell=False)
	std_error = proc_7.communicate()

def main(args):
	#picard java
	java_picard="/usr/local/packages/picard-tools-2.18.26/picard.jar"
	#ref_dict
	ref_dict = args.ref.replace(".fasta",".dict")
	#clean database variable
	snpeff_db = str(args.snpeffdb).strip('[]')

	#create output directory if it doesn't exist
	if not os.path.exists('sv_caller_output'):
		os.mkdir('sv_caller_output')

	#calling functions for sv_caller
	ref_check(args.ref, ref_dict)
	if isinstance(args.sortbam, list):
		for bamfile in args.sortbam:
			#create sample-specific directories
			bam_name = os.path.basename(bamfile[:-4])
			spl_dir = "sv_caller_output/%s" % (bam_name)
			if not os.path.exists(spl_dir):
				os.mkdir(spl_dir)
			#run pipeline
			run_gatk(bamfile, args.ref, java_picard, args.threads)
			if args.masking:
				masking(bamfile, args.masking)
				snpeff(snpeff_db, bamfile, args.bamboozledir, args.threads)
			else:
				snpeff(snpeff_db, bamfile, args.bamboozledir, args.threads)
	else:
		#create sample-specific directory
		bam_name = os.path.basename(args.sortbam[:-4])
		spl_dir = "sv_caller_output/%s" % (bam_name)
		if not os.path.exists(spl_dir):
			os.mkdir(spl_dir)
		#run pipeline
		run_gatk(args.sortbam, args.ref, java_picard, args.threads)
		if args.masking:
			masking(args.sortbam, args.masking)
			snpeff(snpeff_db, args.sortbam, args.bamboozledir, args.threads)
		else:
			snpeff(snpeff_db, args.sortbam, args.bamboozledir, args.threads)
