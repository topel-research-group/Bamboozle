#!/usr/bin/env python3
# coding=utf-8

#André Soares - 8/10/2020
#Quick pipeline to filter a VCF by quality and zygosity, apply bgzip and tabix

import argparse
import os
import subprocess

#read arguments
parser = argparse.ArgumentParser(description='Provide a VCF, take two filtered, gzipped and indexed ones, one for homologous and another for heterozygous mutations!")

#args in
parser.add_argument("-v", "--input_vcf", \
        nargs= '*', type=str, required=True, \
	help="Input VCF")

def all_things(vcf):

#	INPUT: A VCF file, fresh out of sv_caller.py.
#
#	OUTPUT: So much more. Two FILTERED VCF files, one with homozyguous
#	LOF mutations, another with heterozyguous LOF mutations. Both
#	gzipped and indexed.

	vcf_name = vcf[:-4]
	filt_hom_vcf = vcf_name+'_filt_hom.vcf'
	filt_het_vcf = vcf_name+'_filt_het.vcf'

	cmd1 = ['java', '-jar', '/usr/local/packages/snpEff/SnpSift.jar', 'filter', \
		'"(( QUAL >= 100) && (DP >= 30) | (countHom() > 0 )))"', vcf]
	with open(filt_hom_vcf, "w") as f:
		proc1 = subprocess.Popen(cmd1, stdout=f, shell=False)
		std_error = proc1.communicate()

	cmd1a = ['bgzip', '-f', filt_hom_vcf]
	proc1a = subprocess.Popen(cmd1a, shell=False)
	std_out, std_error = proc1a.communicate()

	cmd1b = ['tabix', '-f', '-p', 'vcf', filt_hom_vcf+'.gz']
	proc1b = subprocess.Popen(cmd1b, shell=False)
	std_out, std_error = proc1b.communicate()

	cmd2 = ['java', '-jar', '/usr/local/packages/snpEff/SnpSift.jar', 'filter', \
		'"(( QUAL >= 100) && (DP >= 30) | (countHet() > 0 )))"', vcf]
        with open(filt_het_vcf, "w") as f:
                proc2 = subprocess.Popen(cmd2, stdout=f, shell=False)
                std_error = proc2.communicate()

	cmd2a = ['bgzip', '-f', filt_het_vcf]
	proc2a = subprocess.Popen(cmd2a, shell=False)
	std_out, std_error = proc2a.communicate()

	cmd2b = ['tabix', '-f', '-p', 'vcf', filt_het_vcf+'.gz']
	proc2b = subprocess.Popen(cmd2b, shell=False)
	std_out, std_error = proc2b.communicate()

def main():
	all_things(args.vcf)
