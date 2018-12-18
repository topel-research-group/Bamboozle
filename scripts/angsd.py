#!/usr/bin/env python3

# Uses ANGSD for Fst analysis.

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
		myfile.write(add+"%s\n" % n1)

	myfile.close()
	
	directories2 = args.pop2 + '/*/Bowtie2/*.bam'
	bam_list2 = glob.glob(directories2)
	myfile2 = open("bam_list2.txt","w")
	for n2 in bam_list2:
		myfile2.write(add+"%s\n" % n2)

	myfile2.close()
	
	cmd1 = ['angsd', '-b', '../bam_list1.txt', \
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

	cmd2 = ['angsd', '-b', '../bam_list2.txt', \
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

	cmd3 = ('/usr/local/packages/angsd0.918/angsd/misc/realSFS \
		pop1.saf.idx pop2.saf.idx > pop1.pop2.ml') 
	process3 = subprocess.Popen(cmd3, \
		stdout=subprocess.PIPE, \
		shell=True, \
		cwd='ANGSD')
	while process3.wait() is None:
		pass
	process3.stdout.close()

	cmd4 = ['/usr/local/packages/angsd0.918/angsd/misc/realSFS', \
		'pop1.saf.idx', 'pop2.saf.idx', \
		'-sfs', 'pop1.pop2.ml', \
		'-fstout', 'here']
	process4 = subprocess.Popen(cmd4, \
		stdout=subprocess.PIPE, \
		cwd='ANGSD')
	while process4.wait() is None:
		pass
	process4.stdout.close()

	cmd5 = ['/usr/local/packages/angsd0.918/angsd/misc/realSFS', \
		'fst', 'stats', 'here.fst.idx']
	process5 = subprocess.Popen(cmd5, \
		stdout=subprocess.PIPE, \
		cwd='ANGSD')
	while process5.wait() is None:
		pass
	process5.stdout.close()

	cmd6 = ('/usr/local/packages/angsd0.918/angsd/misc/realSFS \
		fst print here.fst.idx > angsd_results.txt')
	process6 = subprocess.Popen(cmd6, \
		stdout=subprocess.PIPE, \
		shell=True, \
		cwd='ANGSD')
	while process6.wait() is None:
		pass
	process6.stdout.close()

	cmd7 = ('''awk -v OFS='\\t' '{$5 = ($4 != 0) ? sprintf("%.6f", $3 / $4) : "UND"}1' \
		angsd_result.txt > angsd_fst_results.txt''')
	process7 = subprocess.Popen(cmd7, \
		stdout=subprocess.PIPE, \
		shell=True, \
		cwd='ANGSD')
	while process7.wait() is None:
		pass
	process7.stdout.close()
	
	cmd8 = ('''awk '{print $1 "\\t" $2 "\\t" $5}' \
		angsd_fst_results.txt > angsd_fst_results_flt.txt''') 
	process8 = subprocess.Popen(cmd8, \
		stdout=subprocess.PIPE, \
		shell=True, \
		cwd='ANGSD')
	while process8.wait() is None:
		pass
	process8.stdout.close()

	cmd9 = ('''cat %s \
		| grep -v "nan" \
		| "awk '{if ($3 >0) print}' \
		| sort -n > %s''') \
		% ('angsd_fst_results_flt.txt', 'angsd_fst.table')
	process9 = subprocess.Popen(cmd9, \
		stdout=subprocess.PIPE, \
		shell=True, \
		cwd='ANGSD')
	while process9.wait() is None:
		pass
	process9.stdout.close()

	cmd10 = ('cat %s | tr "\\t" "," > %s') \
		% ('angsd_fst.table', 'angsd_fst.csv')
	process10 = subprocess.Popen(cmd10, \
		stdout=subprocess.PIPE, \
		shell=True, \
		cwd='ANGSD')
	while process10.wait() is None:
		pass
	process10.stdout.close()


if __name__ == "__main__":
	main()

