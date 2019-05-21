#!/usr/bin/env python3

# Makes a plot using pandas and matplotlib.pyplot.
# Input file should be csv with headers:
# CHROM,POS,VALUE

import argparse
import pandas as pd
import matplotlib.pyplot as plt

from functools import wraps
from time import time

#######################################################################

parser = argparse.ArgumentParser(prog="pandas_fst_plot.py")
parser.add_argument("-i", "--inputfile", \
                required=True, \
                help="csv input file, 3 col, CHROM,POS,VALUE")
parser.add_argument("-y", "--yaxis", \
                required=True, \
                help="Name of header for y axis")
args = parser.parse_args()

#######################################################################

inputfile = args.inputfile
yaxis = args.yaxis

#######################################################################

def timing(function):
    @wraps(function)
    def wrapper(*args, **kwargs):
        start = time()
        result = function(*args, **kwargs)
        end = time()
        print('Elapsed time: {}'.format(end-start))
        return result
    return wrapper

@timing
def plot():
	# Import csv from Fst statstics with vcftools.
	gl = pd.read_csv(inputfile)

	# Optimize memory usage.
	gl_int = gl.select_dtypes(include=['int'])
	converted_int = gl_int.apply(pd.to_numeric,downcast='unsigned')
	gl_float = gl.select_dtypes(include=['float'])
	converted_float = gl_float.apply(pd.to_numeric,downcast='float')
	optimized_gl = gl.copy()
	optimized_gl[converted_int.columns] = converted_int
	optimized_gl[converted_float.columns] = converted_float
	gl_obj = gl.select_dtypes(include=['object']).copy()

	# Convert CHROM column from object to category.
	chrom = gl_obj.CHROM
	chrom_cat = chrom.astype('category')
	converted_obj = pd.DataFrame()

	# If unique values are more than 50% of the 
	# data don't convert to category, 
	# it will not optimize memory usage.
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
	read_and_optimized = pd.read_csv(inputfile, \
					dtype=column_types)

	# Rename read and optimized csv file 
	# from the Fst analysis to "df".
	df = read_and_optimized 
	df['code'] = chrom_cat.cat.codes

	df['ind'] = range(len(df))
	df_grouped = df.groupby(('code'))

	# Dict for the contig names and index number.
	names = dict( enumerate(df['CHROM'].cat.categories ))

	# Make plot of data.
	fig = plt.figure(figsize=(80,20))
	ax = fig.add_subplot(111)
	colors = ['green', 'turquoise', \
		'blue', 'purple', \
		'red', 'orange', \
		'yellow']
	x_labels = []
	x_labels_pos = []
	for num, (name, group) in enumerate(df_grouped):
		group.plot(kind='scatter', x='ind', y=yaxis, \
		color=colors[num % len(colors)], ax=ax)
		x_labels.append(name)
		x_labels_pos.append((group['ind'].iloc[-1] \
		- (group['ind'].iloc[-1] \
		- group['ind'].iloc[0])/2))
		ax.set_xticks(x_labels_pos)
		ax.set_xticklabels(x_labels, \
				rotation='vertical', \
				fontsize=10)
		ax.set_xlim([0, len(df)])
		ax.set_ylim([0, 1])
		ax.set_xlabel('contigs', fontsize=24)
		ax.set_ylabel('Fst value', fontsize=24)
		ax.set_title('Fst', fontsize=40)
		plt.tick_params(axis='x', length=0.01)
	
	# Add legend with key values paired 
	# with the name of the contig.
	legend_list=[]
	for key, value in names.items():
		temp = [key,value]
		legend_list.append(temp)
	
	plt.legend(legend_list,bbox_to_anchor=(1.01, 1), \
					ncol=5, \
					borderaxespad=0)
	plt.tight_layout(pad=7)	

	# Save plot as image.
	plt.savefig("Fst_plot.png")

if __name__ == "__main__":
	plot()
