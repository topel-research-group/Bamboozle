#!/usr/bin/env python3

import sys
import subprocess
import argparse
import time

###################################################################################

parser = argparse.ArgumentParser(description='Obtain statistics regarding percentage coverage from bam files. \
                                              The script gives percentage of positions in an assembly/contig \
                                              with coverage greater than or equal to a given threshold')

#parser.add_argument('-c', '--contig', type=str, nargs='+', help='Gives cov. stats for the specified contigs')
parser.add_argument('-c', '--contig', type=str, nargs=1, help='Gives cov. stats for the specified contig')
parser.add_argument('-b', '--bam', type=str, nargs=1, help='Bam file to take stats from')
parser.add_argument('-t' ,'--threshold', type=int, nargs='?', const=1, default='20', \
                    help='Threshold for calculating cov. percentage; default 20')
parser.add_argument('--dev', help=argparse.SUPPRESS, action="store_true")
parser.add_argument('-v', '--verbose', help='Gives verbose output')

args = parser.parse_args()

###################################################################################


if args.dev == True:
	start_time = time.time()

# Devel. Add try except in separate module
samtools = str('/usr/local/packages/samtools-1.3.1/samtools')

input_bam = args.bam[0]


# Generate a per-base coverage for each contig
# Output = Contig-Position-Coverage, tab-separated; save output to new dictionary cov_stats

cmd = [samtools, 'depth', '-aa', input_bam]
process = subprocess.Popen(cmd, stdout=subprocess.PIPE)

# cov_stats should have a key for each coverage value of threshold and above
# For each line, if coverage is >= threshold, add 1 to the value of that coverage's key
# Then sum all values in the dictionary and print something informative

# Later, this should be scaled for multiple contigs, with a dictionary for each

def coverage_stats():
	if args.verbose == True:
		if args.contig is not None:
			print("Obtaining stats for " + args.contig[0] + " in " + input_bam + "; coverage >+" + str(args.threshold) + "%.")
		else:
			print("Obtaining whole-genome stats for " + input_bam + "; coverage >+" + str(args.threshold) + "%.")
	cov_stats = {}
	num_lines = 0
	with process.stdout as result:
		rows = (line.decode().split('\t') for line in result)
		for row in rows:
			coverage = int(row[2])
			if args.contig is not None:
				if str(row[0]) == args.contig[0]:
					num_lines += 1
					if coverage >= args.threshold:
						if coverage in cov_stats:
							cov_stats[coverage] += 1
						else:
							cov_stats[coverage] = 1
				
			else:
				num_lines += 1
				if coverage >= args.threshold:
					if coverage in cov_stats:
						cov_stats[coverage] += 1
					else:
						cov_stats[coverage] = 1


	print("Length of assembly/contig: " + str(num_lines))
	value = 100.0 / num_lines * sum(cov_stats.values())
	print(str(value) + "% of the assembly has >=" + str(args.threshold) + "x coverage.")

def main():
	coverage_stats()


if __name__ == "__main__":
	main()

if args.dev == True:
	print("Time taken = " + str(time.time() - start_time) + " seconds.")
