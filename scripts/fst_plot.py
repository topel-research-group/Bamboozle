#!/usr/bin/env python3
# coding=utf-8

#André Soares - 12/08/2020

#       This script provides a pipeline to calculate Fst statistics and plot them
#       using PCA.
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
import csv
import pandas as pd

parser = argparse.ArgumentParser(usage="fst_plot.py <command> [options]")

parser.add_argument("-v", "--vcf_list", \
		help="Text file of newline-separated list of HaplotypeCaller-run VCF files")
parser.add_argument("-p", "--populations", \
		help="Text file listing newline-separated text files of individuals for each population. \
		Individual names in each listed text file must correspond to sample name in the corresponding VCF.")
parser.add_argument("-r", "--reference", \
		help="Reference")
parser.add_argument("-o", "--out_name", \
		help="Name for output files")
parser.add_argument("-t", "--threads", \
		help="Number of threads")

args = parser.parse_args()

def prep_input(vcf_list):
	vcf_list = open(vcf_list, 'r')
	list_a = []
	for line in csv.reader(vcf_list):
		list_a.append(str(line)[2:-2])
	in_list_gz = [s + '.gz' for s in list_a]

	for vcf in list_a:
		cmd1 = ['bgzip'] + vcf.split(",")
		proc1 = subprocess.Popen(cmd1, shell=False)
		std_out, std_error = proc1.communicate()
	for gzvcf in in_list_gz:
		cmd2 = ['tabix', '-p', 'vcf'] + gzvcf.split(",")
		proc2 = subprocess.Popen(cmd2, shell=False)
		std_out, std_error = proc2.communicate()
	return in_list_gz

def comb_geno(vcf_list, out_name, reference, pops, threads):

	in_list_gz = prep_input(vcf_list)

	cmd3_2 = []
	for gzvcf in in_list_gz:
		cmd3_2.append('--variant '+gzvcf)

	cmd3_2_f = [vars for gzvcf in cmd3_2 for vars in gzvcf.split()]
	print(cmd3_2_f)

	cmd3_1 = ['gatk', 'CombineGVCFs', \
		'-R', reference]
	cmd3_3 = ['-O', out_name+'.vcf.gz']
	cmd3 = cmd3_1 + cmd3_2_f + cmd3_3

	print(cmd3)

	proc3 = subprocess.Popen(cmd3, \
		shell=False)
	std_out, std_error = proc3.communicate()
#	FILTERING?

#	java_opts = "-Xmx4G -XX:ParallelGCThreads=%s" % (threads)
#	cmd4 = ['gatk', \
#		'--java-options', java_opts, \
#		'GenotypeGVCFs', \
#		'-R', reference, \
#		'--variant', out_name+'.vcf.gz', \
#		'-O', out_name+'_geno.vcf.gz']
#	proc4 = subprocess.Popen(cmd4, \
#		shell=False)
#
#	cmd5_2 = []
#	with open(pops) as infile2:
#		for pop in infile2:
#			cmd5_2.append('--weir-fst-pop '+pop.replace('\n',''))
#
#	cmd5_1 = ['vcftools' \
#		'--gzvcf', out_name+'_geno.vcf.gz', \
#		'--out', out_name+'.vcf', \
#		'--012']
#	cmd5 = cmd5_1 + cmd5_2
#
#	proc5 = subprocess.Popen(cmd5, \
#		shell=False)

#def pca():

if __name__ == "__main__":
#	prep_input(args.vcf_list)
	comb_geno(args.vcf_list, args.out_name, args.reference, \
		args.populations, args.threads)
