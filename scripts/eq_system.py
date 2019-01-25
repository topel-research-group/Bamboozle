#!/usr/bin/env python3

import sys
import subprocess
import argparse
import fnmatch
import os

import pandas as pd


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
args = parser.parse_args()

#######################################################################


def main(args):
	# Import csv files.
	fst1 = pd.read_csv(args.file1)	
	
	fst2 = pd.read_csv(args.file2)	
	
	fst3 = pd.read_csv(args.file3)	
	
	# Select FST column from each file and put in variable.
	new = pd.DataFrame({'value1': fst1['FST'], 'value2' : fst2['FST'], 'value3': fst3['FST']})

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
	del fst1['FST']
	fst1.to_csv('temp.csv', index=False)

	fst2['Time'] = time
	del fst2['FST']
	fst2.to_csv('time.csv', index=False)

	fst3['Location'] = loc 
	del fst3['FST']
	fst3.to_csv('local.csv', index=False)


if __name__ == "__main__":
	main(args)	
