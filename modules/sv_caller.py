#!/usr/bin/env python3
# coding=utf-8

#André Soares - 15/04/2020
#Pseudocode for a script calling SVs and inferring their LOF potential from BAM/SAM files

#Outputs: 
#	1) VCF files of alterations in gene models, with LOF prediction.
#	2) summary tables of alterations that should include type (deletion, insertion, inversion, duplication, other genomic rearrangmements),
#size in bp, quality of call, LOF potential, gene product, genomic location, zygosity, ?...

#5. Run snpEff snpSift to obtain summary tables

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
		cmd3 = "bwa index %s" % (reference)
		proc_3 = subprocess.Popen(cmd3, shell=True, universal_newlines=True)
		std_out, std_error = proc_3.communicate()
		print("bwa-mem indices didn't exist but they sure do now")

#GRIDSS

def gridss(bamfile, reference, threads, java_gridss, assembly_bam_out, vcf_out):
	#create output folder if it doesn't exist
	if not os.path.exists('sv_caller_output'):
		os.makedirs('sv_caller_output')

#	log_file=open('sv_caller_run.log','a')

	cmd4 = "gridss.sh %s -r %s -a %s -o %s -t %s -j %s" % (bamfile, reference, assembly_bam_out, vcf_out, threads, java_gridss)
	proc_4 = subprocess.Popen(cmd4, shell=True)
	
	std_out, std_error = proc_4.communicate()
	
	print("GRIDSS finished calling SVs")

# BEDTOOLS masking of SV calls goes here

def masking(vcf_out, refpil, masked_vcf_out):
	cmd5 = "bedtools intersect -v -b %s -a %s -sorted > %s" % (refpil, vcf_out, masked_vcf_out)
	proc_5 = subprocess.Popen(cmd5, shell=True)
	std_out, std_error = proc_5.communicate()

# Use R script provided by GRIDSS authors to annotate SVs as DEl, INS, etc

def annotate(masked_vcf_out, bam_name, bamboozledir):
	cmd6 = "Rscript --vanilla %s/scripts/bamboozle_sv_caller_qc_sum.R %s %s" % (bamboozledir, masked_vcf_out, bam_name)
	proc_6 = subprocess.Popen(cmd6, shell=True)
	std_out, std_error = proc_6.communicate()

# Checks for dependencies required for snpEff.
def snpeff(snpeffdb1, masked_ann_vcf_out, bamboozledir1, bamboozledir2, snpeffdb2, masked_vcf_out_lof_csv, masked_vcf_out_lof_ann):
	cmd7 = "snpEff eff %s %s -c %s/data/snpeff/snpEff.config -dataDir %s/data/snpeff/%s -csvStats %s > %s" % (snpeffdb1, masked_ann_vcf_out, bamboozledir2, snpeffdb2, masked_vcf_out_lof_csv, masked_vcf_out_lof_ann)

	proc_7 = subprocess.Popen(cmd7, shell=True)
	std_out, std_error = proc_7.communicate()

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

	#clean database variable
	snpeff_db = str(args.snpeffdb).strip('[]')
	
	#location of snpeff config file
	config_file = "%s/data/snpeff/snpEff.config" %s (args.bamboozledir)
	#location of genome file
	genome_loc = "%s/data/snpeff/" %s (args.bamboozledir)

#	#modify snpeff config file
#	with open(config_file, "w") as config:
#		for line in lines:
#			if line.strip().startswith("data.dir"):
#				config.write("data.dir  = " + genome_loc)
#	config.close()

	#calling functions for sv_caller
	ref_check(args.ref)
	gridss(args.bamfile, args.ref, args.threads, java_gridss, assembly_bam_out, vcf_out)
	#only apply masking() if it's been called
	if args.masking:
		masking(vcf_out, args.masking, masked_vcf_out)
		annotate(masked_vcf_out, bam_name, args.bamboozledir)
		snpeff(snpeff_db, masked_ann_vcf_out, args.bamboozledir, args.bamboozledir, snpeff_db, masked_vcf_out_lof_csv, masked_vcf_out_lof_ann)
	else:
		annotate(vcf_out, bam_name, args.bamboozledir)
		snpeff(snpeff_db, masked_ann_vcf_out, args.bamboozledir, args.bamboozledir, snpeff_db, masked_vcf_out_lof_csv, masked_vcf_out_lof_ann)
