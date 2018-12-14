#!/usr/bin/env python3


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
                myfile.write("%s\n" % n1)

        myfile.close()
	
	directories2 = args.pop2 + '/*/Bowtie2/*.bam'
	bam_list2 = glob.glob(directories2)
	myfile2 = open("bam_list2.txt","w")
	for n2 in bam_list2:
		myfile.write("%s\n" % n2)

	myfile2.close()
	
	cmd1 = ['angsd', '-b', 'bam_list1.txt', '-anc', add+ref, '-out', 'pop1', '-dosaf', '1', '-gl', '1']	
	process1 = subprocess.Popen(cmd, \
		stdout=subprocess.PIPE, \
		cwd='ANGSD')

	cmd1 = ['angsd', '-b', 'bam_list2.txt', '-anc', add+ref, '-out', 'pop2', '-dosaf', '1', '-gl', '1']	
	process1 = subprocess.Popen(cmd, \
		stdout=subprocess.PIPE, \
		cwd='ANGSD')

if __name__ == "__main__":
	main()

