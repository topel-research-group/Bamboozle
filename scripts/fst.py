#!/usr/bin/env python3

#	Input from pipeline -> this script 
#	-> output Fst statistics table and graph
#	Version: VCFtools/v0.1.13

#	Copyright (C) 2018 Vilma Canfjorden. vilma.canfjorden@gmail.com
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <https://www.gnu.org/licenses/>.


import sys
import subprocess
import argparse
import fnmatch
import os
import glob
import pandas as pd
import numpy as np
from scipy.stats import uniform
from scipy.stats import randint
import matplotlib.pyplot as plt

#######################################################################

parser = argparse.ArgumentParser(prog="fst.py")
parser.add_argument("-c", "--clean", \
		action="store_true", \
		help="Remove tmp files")
parser.add_argument("-1", "--pop1", \
		required=True, \
		help="Population 1 input directory")
parser.add_argument("-2", "--pop2", \
		required=True, \
		help="Population 2 input directory")
args = parser.parse_args()

#######################################################################

current_directory = os.getcwd()
name = os.path.basename(current_directory)
merged_vcf_pop1 = name + '_merged_pop1.vcf.gz'
merged_vcf_pop2 = name + '_merged_pop2.vcf.gz'
names1 = 'name_1_list.txt'
names2 =  'name_2_list.txt'
indv_txt_pop1 = name + '_indv_names_pop1.txt'
indv_txt_pop2 = name + '_indv_names_pop2.txt'
population_list = 'pop_list.txt'
all_pop_merged = 'all_pop_merged.vcf.gz'
fst_out = 'pop1_pop2'
fst_out_in = '../Populations/pop1_pop2.weir.fst' 
fst_out_flt = 'tmp.pop1_pop2_flt.table'
fst_out_flt_results = 'tmp.pop1_pop2_flt_results.table'
fst_out_flt2_results = 'tmp.pop1_pop2_flt2_results.table'
fst_results_sorted = 'pop1_pop2_flt_results_sorted.table'
fst_results_sorted_csv = 'pop1_pop2_flt_results_sorted.csv'
path_for_plot = 'Fst_stats/'
add = '../'

#######################################################################

# Perform Fst-statistics on snpEff results (gzipped vcf-files).
def main():
	directories = '*/*/Bcftools/*.snpeff_annotated.vcf.gz'
	file_list = glob.glob(directories)
	for f in file_list:
		cmd1 = ['bcftools', 'index', '-c', '-f', f]
		process1 = subprocess.Popen(cmd1, \
			stdout=subprocess.PIPE)
		while process1.wait() is None:
			pass
		process1.stdout.close()

	# Make directory for the merged vcf-files 
	# for population1 and population2.
	if not os.path.exists('Populations'):
		os.makedirs('Populations')

	# Making a list of vcf-files that will be input to bcftools merge 
	# and then merge population1.
	directories2 = args.pop1 + '/*/Bcftools/*.snpeff_annotated.vcf.gz'
	name_list1 = glob.glob(directories2)
	myfile = open("name_1_list.txt","w")
	for n1 in name_list1:
		myfile.write(add + "%s\n" % n1)

	myfile.close()
	cmd2 = ['bcftools', 'merge', \
		'-l', add+names1, \
		'-Oz', '-o', merged_vcf_pop1]   
	process2 = subprocess.Popen(cmd2, \
		stdout=subprocess.PIPE, \
		cwd='Populations')
	while process2.wait() is None:
		pass
	process2.stdout.close()

	# Making a list of vcf-files that will be input to bcftools merge 
	# and then merge population2.
	directories3 = args.pop2 + '/*/Bcftools/*.snpeff_annotated.vcf.gz'
	name_list2 = glob.glob(directories3)
	myfile2 = open("name_2_list.txt","w")
	for n2 in name_list2:
		myfile2.write(add + "%s\n" % n2)

	myfile2.close()
	cmd3 = ['bcftools', 'merge', \
		'-l', add+names2, \
		'-Oz', '-o', merged_vcf_pop2]   
	process3 = subprocess.Popen(cmd3, \
		stdout=subprocess.PIPE, \
		cwd='Populations')
	while process3.wait() is None:
		pass
	process3.stdout.close()

	# Making a txt file of the names of the individuals in the populations 
	# that is needed for vcftools --wei-fst-pop and indexing 
	# the merged files for population1 and population2.
	for file in os.listdir('Populations'):
		if fnmatch.fnmatch(file, '*_merged_pop1.vcf.gz'):
			cmd4 = ['bcftools', 'index', '-c', '-f', merged_vcf_pop1]
			process4 = subprocess.Popen(cmd4, \
				stdout=subprocess.PIPE, \
				cwd='Populations')
			while process4.wait() is None:
				pass
			process4.stdout.close()

			cmd5 = ('bcftools query -l %s > %s') \
				% (merged_vcf_pop1, indv_txt_pop1) 
			process5 = subprocess.Popen(cmd5, \
				stdout=subprocess.PIPE, \
				shell=True, \
				cwd='Populations')
			while process5.wait() is None:
				pass
			process5.stdout.close()

		elif fnmatch.fnmatch(file, '*_merged_pop2.vcf.gz'):
			cmd6 = ['bcftools', 'index', '-c', '-f', merged_vcf_pop2]
			process6 = subprocess.Popen(cmd6, \
				stdout=subprocess.PIPE, \
				cwd='Populations')
			while process6.wait() is None:
				pass
			process6.stdout.close()

			cmd7 = ('bcftools query -l %s > %s') \
				% (merged_vcf_pop2, indv_txt_pop2)
			process7 = subprocess.Popen(cmd7, \
				stdout=subprocess.PIPE, \
				shell=True, \
				cwd='Populations')
			while process7.wait() is None:
				pass
			process7.stdout.close()

	
	# Making a list of vcf-files that will be input to bcftools merge 
	# and then merge population1 and population2 to an "all_merged" 
	# vcf file, this file will be the input file to 
	# vcftools --weir-fst-pop. 
	directories4 = 'Populations/*_merged_*.vcf.gz'
	pop_list = glob.glob(directories4)
	myfile3 = open("pop_list.txt","w")
	for p in pop_list:
		myfile3.write(add + "%s\n" % p)

	myfile3.close()
	cmd8 = ['bcftools', 'merge', \
		'-l', add+population_list, \
		'-Oz', '-o', all_pop_merged] 
	process8 = subprocess.Popen(cmd8, \
		stdout=subprocess.PIPE, \
		cwd='Populations')
	while process8.wait() is None:
		pass
	process8.stdout.close()

	# Making directory for Fst-results, 
	# input-files to pandas and matplotlib.
	if not os.path.exists('Fst_stats'):
		os.makedirs('Fst_stats')

	# Fst_statistics using vcftools --weir-fst-pop, input files are a
	# vcf file with all merged populations one txt file with names 
	# of the individuals from population1 and one txt file with 
	# names of the individulas from population2, output is a table 
	# of Fst values and a log file of the results. 
	for file in os.listdir('Populations'):
		if fnmatch.fnmatch(file, 'all_pop_merged.vcf.gz'):
			cmd9 = ['vcftools', \
				'--gzvcf', all_pop_merged, \
				'--weir-fst-pop', indv_txt_pop1, \
				'--weir-fst-pop', indv_txt_pop2, \
				'--out', fst_out]
			process9 = subprocess.Popen(cmd9, \
				stdout=subprocess.PIPE, \
				cwd='Populations')
			while process9.wait() is None:
				pass
			process9.stdout.close()

	# Filtering the resulting files from vcftools and making a new 
	# directory called 'Fst_stats' with the resulting files, 
	# output is a csv file and a tab separated table, the csv file 
	# will be the input file to pandas, matplotlib and highcharts. 
	for file in os.listdir('Populations'):
		if fnmatch.fnmatch(file, '*.weir.fst'):
			cmd10 = ('cat %s | grep -v "nan" > %s') \
				% (fst_out_in, fst_out_flt)
			process10 = subprocess.Popen(cmd10, \
				stdout=subprocess.PIPE, \
				shell=True, \
				cwd='Fst_stats')
			while process10.wait() is None:
				pass
			process10.stdout.close()

	# Removing the results below zero.
	for file in os.listdir('Fst_stats'):
		if fnmatch.fnmatch(file, '*flt.table'):
			cmd11 = ("awk '{if ($3 >0) print}' %s > %s") \
				% (fst_out_flt, fst_out_flt_results)
			process11 = subprocess.Popen(cmd11, \
				stdout=subprocess.PIPE, \
				shell=True, \
				cwd='Fst_stats')
			while process11.wait() is None:
				pass
			process11.stdout.close()

			# Rearrange columns (if needed).
			cmd12 = ('''awk '{print $1 "\\t" $2 "\\t" $3}' %s > %s''') \
			% (fst_out_flt_results, fst_out_flt2_results)
			process12 = subprocess.Popen(cmd12, \
				stdout=subprocess.PIPE, \
				shell=True, \
				cwd='Fst_stats')
			while process12.wait() is None:
				pass
			process12.stdout.close()

			# Sorting the POS column (needed for x-axis in highcharts).
			cmd13 = ("cat %s | sort -n > %s") \
				% (fst_out_flt2_results, fst_results_sorted)
			process13 = subprocess.Popen(cmd13, \
				stdout=subprocess.PIPE, \
				shell=True, \
				cwd='Fst_stats')
			while process13.wait() is None:
				pass
			process13.stdout.close()

			# Making a csv file.
			cmd14 = ('cat %s | tr "\\t" ","  > %s') \
				% (fst_results_sorted, fst_results_sorted_csv)
			process14 = subprocess.Popen(cmd14, \
				stdout=subprocess.PIPE, \
				shell=True, \
				cwd='Fst_stats')
			while process14.wait() is None:
				pass
			process14.stdout.close()

	# Making a plot of the Fst results using pandas and matplotlib, 
	# input is the csv file and the output is a pdf file with the plot.
	for file in os.listdir('Fst_stats'):
		if fnmatch.fnmatch(file, 'pop1_pop2_flt_results_sorted.csv'):
			# Import csv file with Fst results.
			gl = pd.read_csv('Fst_stats/pop1_pop2_flt_results_sorted.csv')

			# Optimize memory usage.
			gl_int = gl.select_dtypes(include=['int'])
			converted_int = gl_int.apply(pd.to_numeric,downcast='unsigned')
			gl_float = gl.select_dtypes(include=['float'])
			converted_float = gl_float.apply(pd.to_numeric,downcast='float')
			optimized_gl = gl.copy()
			optimized_gl[converted_int.columns] = converted_int
			optimized_gl[converted_float.columns] = converted_float

			# Convert CHROM column from object to category.
			gl_obj = gl.select_dtypes(include=['object']).copy()
			chrom = gl_obj.CHROM
			chrom_cat = chrom.astype('category')
			converted_obj = pd.DataFrame()

			# If unique values are more than 50% of the data do not 
			# convert to category, it will not optimize memory usage.
			for col in gl_obj.columns:
				num_unique_values = len(gl_obj[col].unique())
				num_total_values = len(gl_obj[col])
				if num_unique_values / num_total_values < 0.5:
					converted_obj.loc[:,col] = gl_obj[col].astype('category')
				else:
					converted_obj.loc[:,col] = gl_obj[col]

			# Apply on the csv file.     
			optimized_gl[converted_obj.columns] = converted_obj
			dtypes_col = optimized_gl.dtypes.index
			dtypes_type = [i.name for i in optimized_gl.dtypes.values]
			column_types = dict(zip(dtypes_col, dtypes_type))
			read_and_optimized = pd.read_csv('Fst_stats/pop1_pop2_flt_results_sorted.csv', \
							 dtype=column_types)

			# Rename the read and optimized csv file 
			# from the Fst analysis to "df".
			df = read_and_optimized
			df['code'] = chrom_cat.cat.codes

			# Make plot of data.
			df['ind'] = range(len(df))
			df_grouped = df.groupby(('code'))
			fig = plt.figure(figsize=(80, 20))
			ax = fig.add_subplot(111)
			colors = ['green', 'turquoise', \
				'blue', 'purple', \
				'red', 'orange', \
				'yellow']
			x_labels = []
			x_labels_pos = []
			for num, (name, group) in enumerate(df_grouped):
				group.plot(kind='scatter', x='ind', y='WEIR_AND_COCKERHAM_FST', \
				color=colors[num % len(colors)], ax=ax)
				x_labels.append(name)
				x_labels_pos.append((group['ind'].iloc[-1] \
				- (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
				ax.set_xticks(x_labels_pos)
				ax.set_xticklabels(x_labels, rotation='vertical', fontsize=10)
				ax.set_xlim([0, len(df)])
				ax.set_ylim([0, 1])
				ax.set_xlabel('contigs', fontsize=24)
				ax.set_ylabel('Fst value', fontsize=24)
				ax.set_title('Weir and Cockerham Fst', fontsize=40)
				plt.tick_params(axis='x', length=0.01)

			# Save plot as pdf. 
			plt.savefig("Fst_stats/Fst_plot.pdf")
	

	# Removing tmp-files.
	if args.clean:
		for textfile in os.listdir('.'):
			if fnmatch.fnmatch(textfile, 'name_*_list.txt'):
				os.remove(textfile)
			elif fnmatch.fnmatch(textfile, 'pop_list.txt'):
				os.remove(textfile)
		
		for fstfile in os.listdir('Fst_stats'):
			if fnmatch.fnmatch(fstfile, 'tmp.*'):
				os.remove('Fst_stats/' + fstfile)
		

if __name__ == "__main__":
	main()
