#!/usr/bin/env python3

# Uses ANGSD for Fst analysis. 
# Version: angsd/v0.918

import sys
import subprocess
import argparse
import fnmatch
import os
import glob

#######################################################################

parser = argparse.ArgumentParser(prog="angsd.py")
parser.add_argument("-c", "--clean", \
                action="store_true", \
                help="Remove tmp files")
parser.add_argument("-1", "--pop1", \
                required=True, \
                help="Population 1 input directory")
parser.add_argument("-2", "--pop2", \
                required=True, \
                help="Population 2 input directory")
parser.add_argument("-f", "--ref", \
                required=True, \
                help="Reference")
args = parser.parse_args()

#######################################################################

ref = args.ref
fst_print = 'tmp.angsd_results.txt' 
fst_results = 'tmp.angsd_fst_results_flt.txt'
fst_col = 'tmp.angsd_fst_results_col5.txt'
fst_flt = 'tmp.angsd_fst.table'
fst_headers = 'angsd_fst_headers.table'
fst_csv = 'angsd_fst_headers.csv'
add = '../'

#######################################################################

def main():
	# Make directory for angsd output. 
	if not os.path.exists('ANGSD'):
		os.makedirs('ANGSD')

	directories1 = args.pop1 + '/*/Bowtie2/*.bam'
	bam_list1 = glob.glob(directories1)
	myfile = open("bam_list1.txt","w")
	for n1 in bam_list1:
		myfile.write(add+"%s\n" % n1)

	myfile.close()
	
	directories2 = args.pop2 + '/*/Bowtie2/*.bam'
	bam_list2 = glob.glob(directories2)
	myfile2 = open("bam_list2.txt","w")
	for n2 in bam_list2:
		myfile2.write(add+"%s\n" % n2)

	myfile2.close()

	# Calculate per pop site allele freq.
	cmd1 = ['angsd', '-P', '$NSLOTS', '-b', '../bam_list1.txt', \
		'-anc', add+ref, \
		'-out', 'pop1', \
		'-dosaf', '1', '-gl', \
		'1']	
	process1 = subprocess.Popen(cmd1, \
		stdout=subprocess.PIPE, \
		cwd='ANGSD')
	while process1.wait() is None:
		pass
	process1.stdout.close()

	cmd2 = ['angsd', '-P', '$NSLOTS', '-b', '../bam_list2.txt', \
		'-anc', add+ref, \
		'-out', 'pop2', \
		'-dosaf', '1', \
		'-gl', '1']	
	process2 = subprocess.Popen(cmd2, \
		stdout=subprocess.PIPE, \
		cwd='ANGSD')
	while process2.wait() is None:
		pass	
	process2.stdout.close()

	# Calculate 2dsfs prior.
	cmd3 = ('/usr/local/packages/angsd0.918/angsd/misc/realSFS \
		pop1.saf.idx pop2.saf.idx -P %s > pop1.pop2.ml') % ('$NSLOTS')
	process3 = subprocess.Popen(cmd3, \
		stdout=subprocess.PIPE, \
		shell=True, \
		cwd='ANGSD')
	while process3.wait() is None:
		pass
	process3.stdout.close()

	# Index and prepare for fst analysis and easy sliding window analysis.
	cmd4 = ['/usr/local/packages/angsd0.918/angsd/misc/realSFS', \
		'fst', 'index', 'pop1.saf.idx', 'pop2.saf.idx', \
		'-sfs', 'pop1.pop2.ml', \
		'-fstout', 'pop1.pop2']
	process4 = subprocess.Popen(cmd4, \
		stdout=subprocess.PIPE, \
		cwd='ANGSD')
	while process4.wait() is None:
		pass
	process4.stdout.close()

	# Global Fst estimate,log = Fst.Unweight Fst.Weight.
	cmd5 = ('/usr/local/packages/angsd0.918/angsd/misc/realSFS \
		fst stats pop1.pop2.fst.idx > angsd_fst_estimate.log')
	process5 = subprocess.Popen(cmd5, \
		stdout=subprocess.PIPE, \
		shell=True, \
		cwd='ANGSD')
	while process5.wait() is None:
		pass
	process5.stdout.close()

	# Print stdout of Fst analysis to file, tab-separated.
	# Columns: CHROM, POS, (a), (a+b).
	cmd6 = ('/usr/local/packages/angsd0.918/angsd/misc/realSFS \
		fst print %s > %s') \
		% ('pop1.pop2.fst.idx', fst_print)
	process6 = subprocess.Popen(cmd6, \
		stdout=subprocess.PIPE, \
		shell=True, \
		cwd='ANGSD')
	while process6.wait() is None:
		pass
	process6.stdout.close()

	# Divide col 3 and 4 (a/(a+b)) and print in col 5 (=Fst value).
	cmd7 = ('''awk -v OFS='\\t' '{$5 = ($4 != 0) ? sprintf("%.6f", $3 / $4) : "nan"}1' \
		tmp.angsd_results.txt > tmp.angsd_fst_results_flt.txt''') 
	process7 = subprocess.Popen(cmd7, \
		stdout=subprocess.PIPE, \
		shell=True, \
		cwd='ANGSD')
	while process7.wait() is None:
		pass
	process7.stdout.close()
	
	# Filtering and preparation of Fst table.
	cmd8 = ('''awk '{print $1 "\\t" $2 "\\t" $5}' \
		%s > %s''') \
		% (fst_results, fst_col) 
	process8 = subprocess.Popen(cmd8, \
		stdout=subprocess.PIPE, \
		shell=True, \
		cwd='ANGSD')
	while process8.wait() is None:
		pass
	process8.stdout.close()

	cmd9 = ('''cat %s \
		| grep -v "nan" \
		| awk '{if ($3 >0) print}' \
		> %s''') \
		% (fst_col, fst_flt)
	process9 = subprocess.Popen(cmd9, \
		stdout=subprocess.PIPE, \
		shell=True, \
		cwd='ANGSD')
	while process9.wait() is None:
		pass
	process9.stdout.close()

	cmd11 = ('echo -e "CHROM\\tPOS\\tFST" | cat - %s > %s') \
		% (fst_flt, fst_headers)
	process11 = subprocess.Popen(cmd11, \
		stdout=subprocess.PIPE, \
		shell=True, \
		cwd='ANGSD')
	while process11.wait() is None:
		pass
	process11.stdout.close()

	# Converting to csv file.
	cmd10 = ('cat %s | tr "\\t" "," > %s') \
		% (fst_headers, fst_csv)
	process10 = subprocess.Popen(cmd10, \
		stdout=subprocess.PIPE, \
		shell=True, \
		cwd='ANGSD')
	while process10.wait() is None:
		pass
	process10.stdout.close()

	# Making a plot of the Fst results using pandas and matplotlib, 
	# input is the csv file and the output is a pdf file with the plot.
	for file in os.listdir('ANGSD'):
		if fnmatch.fnmatch(file, 'angsd_fst_headers.csv'):
			# Import csv file with Fst results.
			gl = pd.read_csv('ANGSD/angsd_fst_headers.csv')

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
			read_and_optimized = pd.read_csv('ANGSD/angsd_fst_headers.csv', \
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
			colors = ['green', 'turquoise', 'blue', 'purple', 'red', 'orange', 'yellow']
			x_labels = []
			x_labels_pos = []
			for num, (name, group) in enumerate(df_grouped):
				group.plot(kind='scatter', x='ind', y='FST', \
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
			plt.savefig("ANGSD/Fst_plot.png")

	if args.clean:
		for textfile in os.listdir('.'):
			if fnmatch.fnmatch(textfile, 'bam_list*.txt'):
				os.remove(textfile)

		for fstfile in os.listdir('ANGSD'):
			if fnmatch.fnmatch(fstfile, 'tmp.*'):
				os.remove('ANGSD/' + fstfile)

if __name__ == "__main__":
	main()

