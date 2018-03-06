#!/usr/bin/env python

import sys
import subprocess
import argparse
import time

###################################################################################

parser = argparse.ArgumentParser(description='This script is a work in progress.')
#parser.add_argument('InputBam', metavar='input', type=str, help='Bam file to take stats from.')
#parser.add_argument('Threshold', metavar='threshold', type=int, help='Threshold for calculating coverage percentage.')

args = parser.parse_args()

###################################################################################

start_time = time.time()

samtools = str('/usr/local/packages/samtools-1.3.1/samtools')
#input_bam = sys.argv[1]
input_bam = str('/proj/data17/Skeletonema_marinoi_adaptation_to_warming_project/01_mapping/P8352_101/P8352_101_sorted.bam')
#threshold = sys.argv[2]
threshold = 40


# Generate a per-base coverage for each contig (samtools can apparently achieve this,
# so use this rather than bedtools [reduces no. of dependencies by one])
# Output = Contig-Position-Coverage, tab-separated
# Save output to new dictionary coverage_stats

coverage_stats = {}

cmd = [samtools, 'depth', '-aa', input_bam]
process = subprocess.Popen(cmd, stdout=subprocess.PIPE)

# coverage_stats should have a key for each coverage value of threshold and above
# For each line, if coverage is >= threshold, add 1 to the value of that coverage's key

with process.stdout as result:
	rows = (line.split('\t') for line in result)
	num_lines = 0
	for row in rows:
		coverage = int(row[2])
		if coverage >= threshold:
			if coverage in coverage_stats:
				coverage_stats[coverage] += 1
			else:
				coverage_stats[coverage] = 1
		num_lines += 1
print "Length of assembly: " + str(num_lines)

# Sum all values in the dictionary and print something informative

value = 100.0 / num_lines * sum(coverage_stats.values())

print str(value) + "% of the assembly has >=" + str(threshold) + "x coverage."

print "Time taken = " + str(time.time() - start_time) + " seconds."
