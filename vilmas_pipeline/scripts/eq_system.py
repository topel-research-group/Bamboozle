#!/usr/bin/env python3

import argparse
import time

import pandas as pd
import matplotlib.pyplot as plt
from functools import reduce

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
                help="Input csv file warm_cold")
parser.add_argument("-2", "--file2", \
                required=False, \
                help="Input csv file control_cold")
parser.add_argument("-3", "--file3", \
                required=False, \
                help="Input csv file control_warm")
parser.add_argument("-m", "--mean", \
		action="store_true", \
		help="Output mean")
parser.add_argument("-v", "--vcftools", \
		action="store_true", \
		help="Input from vcftools_fst.py")
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
	
	# Merge csv files on CHROM and POS.
	frames = [fst1, fst2, fst3]
	new = reduce(lambda left,right: pd.merge(left,right,\
			on=['CHROM','POS'],\
			suffixes=('_1', '_2')), \
			frames)

	# Do the calculations (using equation systems).
	# FST_1 = value1, FST_2 = value2, FST = value3.
	temp = (new['FST_1']-new['FST_2']+new['FST'])/2
	loc = new['FST'] - temp
	time = new['FST_2'] - loc

	# If minus value print 0.
	temp[temp < 0] = 0
	loc[loc < 0] = 0
	time[time < 0] = 0	

	new.to_csv('results.csv', index=False)

	# Put results in new csv files.
	new['Temperature'] = temp
	del new['FST_1']
	del new['FST_2']
	del new['FST']
	new.to_csv('temp.csv', index=False)

	new['Location'] = loc 
	del new['Temperature']
	new.to_csv('local.csv', index=False)

	new['Time'] = time
	del new['Location']
	new.to_csv('time.csv', index=False)

	# Plot the results.
	df = pd.DataFrame({'Temperature': temp, 'Location': loc, 'Time': time})
	ax = df.plot(kind='box',
		color=dict(boxes='black', whiskers='black', medians='r', caps='black'),
		flierprops=dict(linestyle='-', linewidth=1.5),
		medianprops=dict(linestyle='-', linewidth=1.5),
		whiskerprops=dict(linestyle='-', linewidth=1.5),
		capprops=dict(linestyle='-', linewidth=1.5),
		showfliers=False, meanline=True, showmeans=True, grid=False, rot=0)
	ax.set_ylabel('Value')
	plt.savefig('environmental_' + fst_name + '_plot.png')

if __name__ == "__main__":
	start_time = time.time()
	main(args)	
	print("--- %s seconds ---" % (time.time() - start_time))	

