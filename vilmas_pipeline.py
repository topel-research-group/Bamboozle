#!/usr/bin/env python

import sys
import subprocess
import argparse
import fnmatch
import os 
													        
parser = argparse.ArgumentParser(prog="ADD-SCRIPT-NAME-HERE")
parser.add_argument("-v", "--verbose", action="store_true", help="Be more verbose")
parser.add_argument("-b", "--bowtie2", action="store_true", help="Run Bowtie2")
parser.add_argument("-s", "--samtools", action="store_true", help="Run Samtools")

args = parser.parse_args()

ref = 'Skeletonema_marinoi_Ref_v1.1.1.fst'
base = 'P8352_102.contigs'
file1 = 'P8352_102_S1_L001_R1_001.fastq.gz'
file2 = 'P8352_102_S1_L001_R2_001.fastq.gz'
sam = 'P8352_102.sam'
bam = 'P8352_102.bam'

def bowtie2():
	cmd1 = ['bowtie2-build', ref, base]
	process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
	for line in process1.stdout:
		print line
	
	for file in os.listdir('.'):
		if fnmatch.fnmatch(file, '*.rev.2.bt2'):
			cmd2 = ['bowtie2', '--no-unal', '--very-sensitive', '-x', base, '-1', file1, '-2', file2, '-S', sam]	
			process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE)

def samtools():
	for file in os.listdir('.'):
		if fnmatch.fnmatch(file, '*.sam'):
			cmd3 = ['samtools', 'view', '-b', sam, '>', bam]
			process3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE)
			
def main():
	if args.bowtie2:
		bowtie2()
	elif args.samtools:
		samtools()

if __name__ == "__main__":
    main()

