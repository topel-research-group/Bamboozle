#!/usr/bin/env python

import sys
import subprocess
import argparse
import time

###################################################################################

parser = argparse.ArgumentParser(description='This script is a work in progress.')
parser.add_argument('-a', help='Gives whole-assembly coverage stats.', action="store_true")
parser.add_argument('-c', help='Gives per-contig coverage stats.', action="store_true")
parser.add_argument('InputBam', metavar='input', type=str, help='Bam file to take stats from.')
parser.add_argument('Threshold', metavar='threshold', type=int, help='Threshold for calculating coverage %.')

args = parser.parse_args()

###################################################################################

start_time = time.time()

samtools = str('/usr/local/packages/samtools-1.3.1/samtools')
input_bam = sys.argv[2]
threshold = int(sys.argv[3])

# Generate a per-base coverage for each contig
# Output = Contig-Position-Coverage, tab-separated; save output to new dictionary coverage_stats

cmd = [samtools, 'depth', '-aa', input_bam]
process = subprocess.Popen(cmd, stdout=subprocess.PIPE)

# coverage_stats should have a key for each coverage value of threshold and above
# For each line, if coverage is >= threshold, add 1 to the value of that coverage's key
# Then sum all values in the dictionary and print something informative

def whole_assembly():
	print "Obtaining whole-genome stats for " + input_bam + "; coverage >+" + str(threshold) + "%."
	coverage_stats = {}
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
	value = 100.0 / num_lines * sum(coverage_stats.values())
	print str(value) + "% of the assembly has >=" + str(threshold) + "x coverage."

# Should be as above, but instead save each contig as a dictionary key at its first appearance
# Then add one each time that key appears with a value >= the given threshold
# Finally, run the statistics for each contig and print as above

def per_contig():
	print "Obtaining per-contig stats for " + input_bam + "; coverage >+" + str(threshold) + "%."
	max_length = {}
	GTOET_threshold = {}
	with process.stdout as result:
		rows = (line.split('\t') for line in result)
		for row in rows:
			coverage = int(row[2])
			if str(row[0]) not in max_length:
				max_length[str(row[0])] = 0
			if str(row[0]) not in GTOET_threshold:
				GTOET_threshold[str(row[0])] = 0
			if coverage >= threshold:
				GTOET_threshold[str(row[0])] += 1
			max_length[str(row[0])] += 1
	for key in max_length:
		value = 100.0 / max_length[key] * GTOET_threshold[key]
		print "Length of contig " + key + " = " + str(max_length[key])
		print str(value) + "% of contig " + key + " has >=" + str(threshold) + "x coverage."
		print "---"

if args.a == True:
	whole_assembly()
elif args.c == True:
	per_contig()

print "Time taken = " + str(time.time() - start_time) + " seconds."
