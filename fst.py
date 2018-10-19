#!/usr/bin/env python

import sys
import subprocess
import argparse
import fnmatch
import os
import glob

##################################################################################
parser = argparse.ArgumentParser(prog="ADD-SCRIPT-NAME-HERE")
parser.add_argument("-v", "--verbose", action="store_true", help="Be more verbose")
parser.add_argument("-c", "--clean", action="store_true", help="Remove some files")
args = parser.parse_args()
##################################################################################
current_directory = os.getcwd()
name = os.path.basename(current_directory)
merged_vcf_pop1 = name + '_merged_pop1.vcf.gz'
merged_vcf_pop2 = name + '_merged_pop2.vcf.gz'
names1 = current_directory + '/name_1_list.txt'
names2 = current_directory + '/name_2_list.txt'
indv_txt_pop1 = name + '_indv_names_pop1.txt'
indv_txt_pop2 = name + '_indv_names_pop2.txt'
population_list = current_directory + '/pop_list.txt'
all_pop_merged = 'all_pop_merged.vcf.gz'
fst_out = 'pop1_pop2'
fst_out_in = current_directory + '/Populations/pop1_pop2.weir.fst' 
fst_out_flt = 'tmp.pop1_pop2_flt.table'
fst_out_flt_results = 'tmp.pop1_pop2_flt_results.table'
fst_out_flt2_results = 'tmp.pop1_pop2_flt2_results.table'
fst_results_sorted = 'pop1_pop2_flt_results_sorted.table'
fst_results_sorted_csv = 'pop1_pop2_flt_results_sorted.csv'
##################################################################################

# Perform Fst-statistics on gziped vcf-files
def main():
	directories = current_directory + '/*/*/Bcftools/*.bcftools_filtered.vcf.gz'
	file_list = glob.glob(directories)
	for f in file_list:
		cmd1 = ['bcftools', 'index', '-c', f]
		process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
		while process1.wait() is None:
                	pass

	# Make directory for the merged vcf-files for population1 and population2
	population_directory = os.path.join(current_directory, r'Populations')
        if not os.path.exists(population_directory):
        	os.makedirs(population_directory)

	# Making a list of vcf-files that will be input to bcftools merge and then merge population1
	directories2 = current_directory + '/*_1/*/Bcftools/*.bcftools_filtered.vcf.gz'
        name_list1 = glob.glob(directories2)
	myfile = open("name_1_list.txt","w")
	for n1 in name_list1:
        	myfile.write("%s\n" % n1)

	myfile.close()
       	cmd2 = ['bcftools', 'merge', '-l', names1, '-Oz', '-o', merged_vcf_pop1]   
	process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, cwd='Populations')
	while process2.wait() is None:
            	pass

	# Making a list of vcf-files that will be input to bcftools merge and then merge population2
   	directories3 = current_directory + '/*_2/*/Bcftools/*.bcftools_filtered.vcf.gz'
        name_list2 = glob.glob(directories3)
        myfile2 = open("name_2_list.txt","w")
        for n2 in name_list2:
                myfile2.write("%s\n" % n2)

        myfile2.close()
        cmd3 = ['bcftools', 'merge', '-l', names2, '-Oz', '-o', merged_vcf_pop2]   
        process3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE, cwd='Populations')
        while process3.wait() is None:
                pass

	# Making a txt file of the names of the individuals in the populations that is needed for vcftools --wei-fst-pop
	# and indexing the merged files for population1 and population2
	for file in os.listdir('Populations'):
        	if fnmatch.fnmatch(file, '*_merged_pop1.vcf.gz'):
        		cmd5 = ('bcftools query -l %s > %s') % (merged_vcf_pop1, indv_txt_pop1) 
        		process5 = subprocess.Popen(cmd5, stdout=subprocess.PIPE, shell=True, cwd='Populations')
			while process5.wait() is None:
        			pass

			cmd6 = ['bcftools', 'index', '-c', merged_vcf_pop1]
                        process6 = subprocess.Popen(cmd6, stdout=subprocess.PIPE, cwd='Populations')
                        while process6.wait() is None:
                                pass

		elif fnmatch.fnmatch(file, '*_merged_pop2.vcf.gz'):
                        cmd7 = ('bcftools query -l %s > %s') % (merged_vcf_pop2, indv_txt_pop2)
                        process7 = subprocess.Popen(cmd7, stdout=subprocess.PIPE, shell=True, cwd='Populations')
                        while process7.wait() is None:
                                pass

			cmd8 = ['bcftools', 'index', '-c', merged_vcf_pop2]
                	process8 = subprocess.Popen(cmd8, stdout=subprocess.PIPE, cwd='Populations')
                	while process8.wait() is None:
                		pass
	
	# Making a list of vcf-files that will be input to bcftools merge and then merge population1 and population2
	# to a "all_merged" vcf file, which will be the input file to vcftools --weir-fst-pop 
	directories4 = current_directory + '/Populations/*_merged_*.vcf.gz'
        pop_list = glob.glob(directories4)
	myfile3 = open("pop_list.txt","w")
       	for p in pop_list:
            	myfile3.write("%s\n" % p)

        myfile3.close()
	cmd9 = ['bcftools', 'merge', '-l', population_list, '-Oz', '-o', all_pop_merged] 
       	process9 = subprocess.Popen(cmd9, stdout=subprocess.PIPE, cwd='Populations')
	while process9.wait() is None:
       		pass

	# Making directory for Fst-results, input-files to highcharts
	fst_directory = os.path.join(current_directory, r'Fst_stats')
        if not os.path.exists(fst_directory):
                os.makedirs(fst_directory)

	# Fst_statistics 
	for file in os.listdir('Populations'):
                if fnmatch.fnmatch(file, 'all_pop_merged.vcf.gz'):
			cmd10 = ['vcftools', '--gzvcf', all_pop_merged, '--weir-fst-pop', indv_txt_pop1, '--weir-fst-pop', indv_txt_pop2, '--out', fst_out]
			process10 = subprocess.Popen(cmd10, stdout=subprocess.PIPE, cwd='Populations')
			while process10.wait() is None:
                                pass

	# Filtering the resulting files from vcftools and making a new directory called Fst_stats with the resulting files,  
	# the csv-file will be the input file to high charts 
	for file in os.listdir('Populations'):
        	if fnmatch.fnmatch(file, '*.weir.fst'):
			cmd11 = ('cat %s | grep -v "nan" > %s') % (fst_out_in, fst_out_flt)
			process11 = subprocess.Popen(cmd11, stdout=subprocess.PIPE, shell=True, cwd='Fst_stats')
			while process11.wait() is None:
        			pass

	# Removing the results below zero
        for file in os.listdir('Fst_stats'):
        	if fnmatch.fnmatch(file, '*flt.table'):
        		cmd12 = ("awk '{if ($3 >=0) print}' %s > %s") % (fst_out_flt, fst_out_flt_results)
        		process12 = subprocess.Popen(cmd12, stdout=subprocess.PIPE, shell=True, cwd='Fst_stats')
			while process12.wait() is None:
        			pass

			# Rearrange columns
                        cmd13 = ('''awk '{print $2 "\t" $3}' %s > %s''') % (fst_out_flt_results, fst_out_flt2_results)
                        process13 = subprocess.Popen(cmd13, stdout=subprocess.PIPE, shell=True, executable='/bin/bash', cwd='Fst_stats')
                        while process13.wait() is None:
                                pass

			# Sorting the POS column (needed for x-axis in highcharts)
                        cmd14 = ("cat %s | sort -n > %s") % (fst_out_flt2_results, fst_results_sorted)
                        process14 = subprocess.Popen(cmd14, stdout=subprocess.PIPE, shell=True, cwd='Fst_stats')
                        while process14.wait() is None:
                                pass

			# Making a csv-file
                        cmd15 = ('cat %s | tr "\t" ","  > %s') % (fst_results_sorted, fst_results_sorted_csv)
                        process15 = subprocess.Popen(cmd15, stdout=subprocess.PIPE, shell=True, executable='/bin/bash', cwd='Fst_stats')
                        while process15.wait() is None:
                                pass

	# Removing tmp-files
	if args.clean:
		for textfile in os.listdir('.'):
			if fnmatch.fnmatch(textfile, 'name_*_list.txt'):
				os.remove(textfile)
			elif fnmatch.fnmatch(textfile, 'pop_list.txt'):
                                os.remove(textfile)
		
		for fstfile in os.listdir('Fst_stats'):
                        if fnmatch.fnmatch(fstfile, 'tmp.*'):
                                os.remove(current_directory + '/Fst_stats/' + fstfile)
		

if __name__ == "__main__":
        main()
