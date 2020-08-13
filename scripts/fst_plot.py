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

parser = argparse.ArgumentParser(usage="fst_plot.py <command> [options]")

parser.add_argument("-v", "--vcf_list", \
		help="Text file of newline-separated list of HaplotypeCaller-run VCF files")
parser.add_argument("-p", "--populations", action="append" \
		help="Text file listing newline-separated text files of individuals for each population. \
		Individual names in each listed text file must correspond to sample name in the corresponding VCF.")
parser.add_argument("-f", "--ref", \
		help="Reference")

def prep_input(vcf_list):
	with open(vcf_list) as infile:
		in_list = csv.reader(infile, delimiter="\n")
		for vcf in in_list:
			cmd1 = ['bgzip', vcf]
			proc1 = subprocess.Popen(cmd1, \
				shell=False, stdout=subprocess.PIPE)
			cmd2 = ['tabix', '-p', 'vcf']
				proc2 = subprocess.Popen(cmd2, \
				shell=False, stdin=proc1.stdout)
	
			proc1.stdout.close()
			std_out, std_error = proc2.communicate()

def comb_geno(vcf_list, out_name, reference, pops, threads):
	cmd2_2 = []
	with open(vcf_list) as infile:
		for vcf in infile:
			cmd2_2.append('--variant '+vcf.replace('\n','')+'.gz')

	cmd2_1 = ['gatk', 'CombineGVCFs', \
		'-R', reference,
		'-O', out_name+'.vcf.gz']
	cmd2 = cmd2_1 + cmd2_2
	proc2 = subprocess.Popen(cmd2, \
		shell=False)

	java_opts = "-Xmx4G -XX:ParallelGCThreads=%s" % (threads)
	cmd3 = ['gatk', \
		'--java-options', java_opts, \
		'GenotypeGVCFs', \
		'-R', reference, \
		'--variant', out_name+'.vcf.gz', \
		'-O', out_name+'_geno.vcf.gz']
	proc3 = subprocess.Popen(cmd3, \
		shell=False)

	cmd4_2 = []
	with open(pops) as infile2:
		for pop in infile2:
			cmd4_2.append('--weir-fst-pop '+pop.replace('\n',''))

	cmd4_1 = ['vcftools' \
		'--gzvcf', out_name+'_geno.vcf.gz', \
		'--out', out_name+'.vcf', \
		'--012']
	cmd4 = cmd4_1 + cmd4_2

	proc4 = subprocess.Popen(cmd4, \
                shell=False)

def pca():
