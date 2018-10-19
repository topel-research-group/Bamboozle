#!/usr/bin/env python

import sys
import subprocess
import fnmatch
import os
import glob

current_directory = os.getcwd()
name = os.path.basename(current_directory)

#directory = current_directory + '/P8352_103/Bcftools'
#directory2 = current_directory + '/P8352_108/Bcftools'
#bcftools_out = current_directory + '/P8352_103/Bcftools/P8352_103.bcftools_filtered.vcf.gz'
merged_vcf_pop1 = name + '_merged_pop1.vcf.gz'
merged_vcf_pop2 = name + '_merged_pop2.vcf.gz'
names1 = current_directory + '/name_1_list.txt'
names2 = current_directory + '/name_2_list.txt'
indv_txt_pop1 = name + '_indv_names_pop1.txt'
indv_txt_pop2 = name + '_indv_names_pop2.txt'
population_list = current_directory + '/pop_list.txt'
all_pop_merged = 'all_pop_merged.vcf.gz'
fst_out = 'pop1_pop2' 
fst_out_flt = 'pop1_pop2_flt.table'
fst_out_flt_results = 'pop1_pop2_flt_results.table'

def main():
	directories = current_directory + '/*/*/Bcftools/*.bcftools_filtered.vcf.gz'
	file_list = glob.glob(directories)
	for f in file_list:
		cmd1 = ['bcftools', 'index', '-c', f]
		process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
		while process1.wait() is None:
                	pass

	# Make directory for the merge vcf-files for one population
	population_directory = os.path.join(current_directory, r'Population')
        if not os.path.exists(population_directory):
        	os.makedirs(population_directory)

	directories2 = current_directory + '/*_1/*/Bcftools/*.bcftools_filtered.vcf.gz'
        name_list1 = glob.glob(directories2)
	myfile = open("name_1_list.txt","w")
	for n1 in name_list1:
        	myfile.write("%s\n" % n1)

	myfile.close()
       	cmd2 = ('bcftools merge -l %s -O z -o %s') % (names1, merged_vcf_pop1)   
	process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell=True, cwd='Population')
	while process2.wait() is None:
            	pass

   	directories3 = current_directory + '/*_2/*/Bcftools/*.bcftools_filtered.vcf.gz'
        name_list2 = glob.glob(directories3)
        myfile2 = open("name_2_list.txt","w")
        for n2 in name_list2:
                myfile2.write("%s\n" % n2)

        myfile2.close()
        cmd3 = ('bcftools merge -l %s -O z -o %s') % (names2, merged_vcf_pop2)   
        process3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE, shell=True, cwd='Population')
        while process3.wait() is None:
                pass

	for file in os.listdir('Population'):
        	if fnmatch.fnmatch(file, '*_merged_pop1.vcf.gz'):
        		cmd5 = ('bcftools query -l %s > %s') % (merged_vcf_pop1, indv_txt_pop1) 
        		process5 = subprocess.Popen(cmd5, stdout=subprocess.PIPE, shell=True, cwd='Population')
			while process5.wait() is None:
        			pass

			cmd6 = ['bcftools', 'index', '-c', merged_vcf_pop1]
                        process6 = subprocess.Popen(cmd6, stdout=subprocess.PIPE, cwd='Population')
                        while process6.wait() is None:
                                pass

		elif fnmatch.fnmatch(file, '*_merged_pop2.vcf.gz'):
                        cmd7 = ('bcftools query -l %s > %s') % (merged_vcf_pop2, indv_txt_pop2)
                        process7 = subprocess.Popen(cmd7, stdout=subprocess.PIPE, shell=True, cwd='Population')
                        while process7.wait() is None:
                                pass

			cmd8 = ['bcftools', 'index', '-c', merged_vcf_pop2]
                	process8 = subprocess.Popen(cmd8, stdout=subprocess.PIPE, cwd='Population')
                	while process8.wait() is None:
                		pass
	
	directories4 = current_directory + '/Population/*_merged_*.vcf.gz'
        pop_list = glob.glob(directories4)
	myfile3 = open("pop_list.txt","w")
       	for p in pop_list:
            	myfile3.write("%s\n" % p)

        myfile3.close()
	cmd9 = ('bcftools merge -l %s -O z -o %s') % (population_list, all_pop_merged) 
       	process9 = subprocess.Popen(cmd9, stdout=subprocess.PIPE, shell=True, cwd='Population')
	while process9.wait() is None:
       		pass

	fst_directory = os.path.join(current_directory, r'Fst_stats')
        if not os.path.exists(fst_directory):
                os.makedirs(fst_directory)

	for file in os.listdir('Population'):
                if fnmatch.fnmatch(file, 'all_pop_merged.vcf.gz'):
			cmd10 = ['vcftools', '--gzvcf', all_pop_merged, '--weir-fst-pop', indv_txt_pop1, '--weir-fst-pop', indv_txt_pop2, '--out', fst_out]
			process10 = subprocess.Popen(cmd10, stdout=subprocess.PIPE, cwd='Fst_stats')
			while process10.wait() is None:
                                pass

	for file in os.listdir('Fst_stat'):
        	if fnmatch.fnmatch(file, '*.weir.fst'):
			cmd11 = ('cat %s | grep -v "nan" > %s') % (fst_out, fst_out_flt)
			process11 = subprocess.Popen(cmd11, stdout=subprocess.PIPE, shell=True, cwd='Fst_stat')
			while process11.wait() is None:
        			pass

        for file in os.listdir('Fst_stat'):
        	if fnmatch.fnmatch(file, '*flt.table'):
        		cmd12 = ("awk '{if ($3 >=0) print}' %s > %s") % (fst_out_flt, fst_out_flt_results)
        		process12 = subprocess.Popen(cmd12, stdout=subprocess.PIPE, shell=True, cwd='Fst_stat')
			while process12.wait() is None:
        			pass


if __name__ == "__main__":
        main()
