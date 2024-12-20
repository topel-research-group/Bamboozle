#!/usr/bin/env python3

#	Input from pipeline -> this script 
#	-> output Fst statistics table and graph
#	Requires relative good coverage, if the
#	coverage is low use angsd_fst.py instead.
#	Version: VCFtools/v0.1.13

import sys
import subprocess
import argparse
import fnmatch
import os
import glob

import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('agg')

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
parser.add_argument("--feature", \
		help="Input gff feature")
parser.add_argument("-w", "--window", \
		help="Size of window")
parser.add_argument("-s", "--step", \
		help="Size of step")
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
fst_out_window_in = '../Populations/pop1_pop2.windowed.weir.fst' 
fst_out_flt = 'tmp.pop1_pop2_flt.table'
fst_out_flt_results = 'tmp.pop1_pop2_flt_results.table'
fst_out_flt2_results = 'tmp.pop1_pop2_flt2_results.table'
fst_results_sorted = 'pop1_pop2_flt_results_sorted.table'
fst_results_sorted_csv = 'pop1_pop2_flt_results_sorted.csv'
fst_calculated = 'tmp.calculated.table'
rm_headers = 'tmp.remove_headers.table'
fst_headers = 'tmp.new_headers.table'
path_for_plot = 'Fst_stats/'
add = '../'

# Specify input file for the analysis by using '--feature'.
if args.feature:
	filename = '*' + args.feature +'_hdr_snpeff_annotated.vcf.gz'
else:
	filename = '*.snpeff_annotated.vcf.gz'

#######################################################################

# Perform Fst-statistics on snpEff results (gzipped vcf-files).
def vcftools():
	directories = args.pop1 + '/*/Bcftools/' + filename
	file_list = glob.glob(directories)
	for f in file_list:
		cmd1 = ['bcftools', 'index', '-c', '-f', f]
		process1 = subprocess.Popen(cmd1, \
			stdout=subprocess.PIPE)
		while process1.wait() is None:
			pass
		process1.stdout.close()

	directories1 = args.pop2 + '/*/Bcftools/' + filename
	file_list1 = glob.glob(directories1)
	for f1 in file_list1:
		cmd_1 = ['bcftools', 'index', '-c', '-f', f1]
		process_1 = subprocess.Popen(cmd_1, \
			stdout=subprocess.PIPE)
		while process_1.wait() is None:
			pass
		process_1.stdout.close()

	# Make directory for the merged vcf-files 
	# for population1 and population2.
	if not os.path.exists('Populations'):
		os.makedirs('Populations')

	# Making a list of vcf-files that will be input to bcftools merge 
	# and then merge population1.
	directories2 = args.pop1 + '/*/Bcftools/' + filename
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
	directories3 = args.pop2 + '/*/Bcftools/' + filename
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

def fst():
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

def sliding_window():
	# If using sliding window.
	for file in os.listdir('Populations'):
		if fnmatch.fnmatch(file, 'all_pop_merged.vcf.gz'):
			cmda = ['vcftools', \
				'--gzvcf', all_pop_merged, \
				'--weir-fst-pop', indv_txt_pop1, \
				'--weir-fst-pop', indv_txt_pop2, \
				'--fst-window-size', args.window, \
				'--fst-window-step', args.step, \
				'--out', fst_out]
			processa = subprocess.Popen(cmda, \
				stdout=subprocess.PIPE, \
				cwd='Populations')
			while processa.wait() is None:
				pass
			processa.stdout.close()

	for file in os.listdir('Populations'):
		if fnmatch.fnmatch(file, '*.weir.fst'):
			cmdb = ('cat %s | grep -v "nan" > %s') \
				% (fst_out_window_in, fst_out_flt)
			processb = subprocess.Popen(cmdb, \
				stdout=subprocess.PIPE, \
				shell=True, \
				cwd='Fst_stats')
			while processb.wait() is None:
				pass
			processb.stdout.close()

	# Removing the results below zero.
	for file in os.listdir('Fst_stats'):
		if fnmatch.fnmatch(file, '*flt.table'):
			cmdc = ("awk '{if ($6 >0) print}' %s > %s") \
				% (fst_out_flt, fst_out_flt_results)
			processc = subprocess.Popen(cmdc, \
				stdout=subprocess.PIPE, \
				shell=True, \
				cwd='Fst_stats')
			while processc.wait() is None:
				pass
			processc.stdout.close()
			# Calculating the midpoint position when using sliding window,
			# column (2+3)/2 and the output in column 2.
			cmde = ('''awk '{a=int(($2+$3)/2); $2=a; print}' %s > %s''') \
				% (fst_out_flt_results, fst_out_flt2_results)
			processe = subprocess.Popen(cmde, \
				stdout=subprocess.PIPE, \
				shell=True, \
				cwd='Fst_stats')
			while processe.wait() is None:
				pass
			processe.stdout.close()

			# Rearrange columns (if needed) and keep midpoint value in column 2.
			cmdd = ('''awk '{print $1 "\\t" $2 "\\t" $6}' %s > %s''') \
				% (fst_out_flt2_results, fst_calculated)
			processd = subprocess.Popen(cmdd, \
				stdout=subprocess.PIPE, \
				shell=True, \
				cwd='Fst_stats')
			while processd.wait() is None:
				pass
			processd.stdout.close()

			# Remove the header provided when using sliding window.
			cmdf = ('echo "$(tail -n +2 %s)" > %s') \
				% (fst_calculated, rm_headers) 
			processf = subprocess.Popen(cmdf, \
				stdout=subprocess.PIPE, \
				shell=True, \
				cwd='Fst_stats')
			while processf.wait() is None:
				pass
			processf.stdout.close()

			# Add headers. 
			cmdg = ('echo -e "CHROM\\tPOS\\tWEIR_AND_COCKERHAM_FST" \
				| cat - %s > %s') \
				% (rm_headers, fst_headers)
			processg = subprocess.Popen(cmdg, \
				stdout=subprocess.PIPE, \
				shell=True, \
				cwd='Fst_stats')
			while processg.wait() is None:
				pass
			processg.stdout.close()

			# Sorting the POS column (needed for x-axis in highcharts).
			cmdi = ("cat %s | sort -n > %s") \
				% (fst_headers, fst_results_sorted)
			processi = subprocess.Popen(cmdi, \
				stdout=subprocess.PIPE, \
				shell=True, \
				cwd='Fst_stats')
			while processi.wait() is None:
				pass
			processi.stdout.close()

			# Making a csv file.
			cmdj = ('cat %s | tr "\\t" ","  > %s') \
				% (fst_results_sorted, fst_results_sorted_csv)
			processj = subprocess.Popen(cmdj, \
				stdout=subprocess.PIPE, \
				shell=True, \
				cwd='Fst_stats')
			while processj.wait() is None:
				pass
			processj.stdout.close()

# Making a plot of the Fst results using pandas and matplotlib, 
# input is the csv file and the output is a pdf file with the plot.
def plot():
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

			df['ind'] = range(len(df))
			df_grouped = df.groupby(('code'))

			# Dict for the contig names and index number.
			names = dict( enumerate(df['CHROM'].cat.categories ))

			# Make plot of data.
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

			# Add legend with key values paired with the name of the contig.
			legend_list=[]
			for key, value in names.items():
				temp = [key,value]
				legend_list.append(temp)

			plt.legend(legend_list,bbox_to_anchor=(1.01, 1), \
						ncol=5, \
						borderaxespad=0)
			plt.tight_layout(pad=7)

			# Save plot as pdf. 
			plt.savefig("Fst_stats/Fst_plot_vcftools.pdf")
def main():
	# Making directory for Fst-results, 
	# input-files to pandas and matplotlib.
	if not os.path.exists('Fst_stats'):
		os.makedirs('Fst_stats')

	if args.window and args.step:
		try:
			sliding_window()
			plot()
		except:
			vcftools()
			sliding_window()
			plot()
	else:
		vcftools()
		fst()
		plot()	

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
