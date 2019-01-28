#!/usr/bin/env python3

import sys
import subprocess
import argparse
import fnmatch
import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#	Temperature + Time = value1 (File1) 
#	Time + Location = value2 (File2)
#	Location + Temperature = value3 (File3)
#	
#	-> Temperature = value1 - Time
#	-> Time = value2 - Location
#	-> Location = value3 - Temperature

#######################################################################

parser = argparse.ArgumentParser(prog="angsd.py")
parser.add_argument("-1", "--file1", \
                required=False, \
                help="Input file warm_cold")
parser.add_argument("-2", "--file2", \
                required=False, \
                help="Input file control_cold")
parser.add_argument("-3", "--file3", \
                required=False, \
                help="Input file control_warm")
parser.add_argument("-m", "--mean", \
		action="store_true", \
		help="Output mean")
parser.add_argument("-v", "--vcftools", \
		action="store_true", \
		help="Input from vcftools_fsy.py")
parser.add_argument("-a", "--angsd", \
		action="store_true", \
		help="Input from angsd_fst.py")
args = parser.parse_args()

#######################################################################

fst_name = ""
if args.angsd:
	fst_name = 'FST'
if args.vcftools:
	fst_name = 'WEIR_AND_COCKERHAM_FST'

def main(args):
	# Import csv files.
	fst1 = pd.read_csv(args.file1)	
	
	fst2 = pd.read_csv(args.file2)	
	
	fst3 = pd.read_csv(args.file3)	

	if args.mean:
		mean1 = fst1[fst_name].mean()
		mean2 = fst2[fst_name].mean()
		mean3 = fst3[fst_name].mean()	
	
		# Do the calculations (using equation systems).
		temp = (mean1-mean2+mean3)/2
		loc = mean3 - temp
		time = mean2 - loc

		mean_file = open('Mean_' + fst_name + '.txt', 'w')
		mean_file.write('Temperature:' + str(temp) + '\n' \
				'Location:' + str(loc) + '\n' \
				'Time:' + str(time) + '\n')
		mean_file.close()

	
	# Select WEIR_AND_COCKERHAM_FST column from each file and put in variable.
	new = pd.DataFrame({'value1': fst1[fst_name], 'value2' : fst2[fst_name], 'value3': fst3[fst_name]})

	# If NaN print 0.
	new = new.fillna(0)

	# Do the calculations (using equation systems).
	temp = (new['value1']-new['value2']+new['value3'])/2
	loc = new['value3'] - temp
	time = new['value2'] - loc

	# If minus value print 0.
	temp[temp < 0] = 0
	loc[loc < 0] = 0
	time[time < 0] = 0	

	# Put results in new csv files.
	fst1['Temperature'] = temp
	del fst1[fst_name]
	fst1.to_csv('temp.csv', index=False)

	fst2['Time'] = time
	del fst2[fst_name]
	fst2.to_csv('time.csv', index=False)

	fst3['Location'] = loc 
	del fst3[fst_name]
	fst3.to_csv('local.csv', index=False)

	# Plot the results.
	df = pd.DataFrame({'Temperature': temp, 'Location': loc, 'Time': time})
	ax = df.plot(kind='box',
		color=dict(boxes='r', whiskers='black', medians='r', caps='black'),
		flierprops=dict(linestyle='-', linewidth=1.5),
		medianprops=dict(linestyle='-', linewidth=1.5),
		whiskerprops=dict(linestyle='-', linewidth=1.5),
		capprops=dict(linestyle='-', linewidth=1.5),
		showfliers=False, grid=False, rot=0)
	ax.set_ylabel('Fst')
	plt.savefig('environmental_' + fst_name + '_plot.png')

if __name__ == "__main__":
	main(args)	
