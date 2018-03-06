#!/usr/bin/env python3

import sys
import subprocess
import argparse
import time

###################################################################################

parser = argparse.ArgumentParser(description='Obtain statistics regarding percentage coverage from bam files. \
                                              The script gives percentage of positions in an assembly/contig \
                                              with coverage greater than or equal to a given threshold')
parser.add_argument('-w', '--whole', help='Gives whole-assembly coverage stats', action="store_true")
parser.add_argument('-c', '--contig', help='Gives per-contig coverage stats', action="store_true")
#parser.add_argument('-e', help='Contigs to ignore during analysis', nargs='*', type=str)
parser.add_argument('Input', type=str, nargs=1, help='Bam file to take stats from')
parser.add_argument('-t' ,'--threshold', type=int, nargs='?', const=1, default='20', help='Threshold for calculating coverage percentage; default 20')
parser.add_argument('--dev', help=argparse.SUPPRESS, action="store_true")

args = parser.parse_args()

###################################################################################

# Devel. This part should be handled by the arg parser
if args.whole == True and args.contig == True:
	print("Please select only one protocol!")
	exit()

if args.whole == False and args.contig == False:
	print("Please select a protocol!")
	print("For details, see help screen.")
	exit()

##########

if args.dev == True:
	start_time = time.time()

# Devel. Add try except in separate mdule
samtools = str('/usr/local/packages/samtools-1.3.1/samtools')
input_bam = args.Input[0]
threshold = args.threshold


#if args.e[0] == True:
#	print(args.e[0])
#	ignore_contigs = [args.e[0:]]
#	print(ignore_contigs)


# Generate a per-base coverage for each contig
# Output = Contig-Position-Coverage, tab-separated; save output to new dictionary coverage_stats

cmd = [samtools, 'depth', '-aa', input_bam]
process = subprocess.Popen(cmd, stdout=subprocess.PIPE)

# coverage_stats should have a key for each coverage value of threshold and above
# For each line, if coverage is >= threshold, add 1 to the value of that coverage's key
# Then sum all values in the dictionary and print something informative

def whole_assembly():
	print("Obtaining whole-genome stats for " + input_bam + "; coverage >+" + str(threshold) + "%.")
	coverage_stats = {}
	with process.stdout as result:
		rows = (line.decode().split('\t') for line in result)
		num_lines = 0
		for row in rows:
			coverage = int(row[2])
			if coverage >= threshold:
				if coverage in coverage_stats:
					coverage_stats[coverage] += 1
				else:
					coverage_stats[coverage] = 1
			num_lines += 1
	print("Length of assembly: " + str(num_lines))
	value = 100.0 / num_lines * sum(coverage_stats.values())
	print(str(value) + "% of the assembly has >=" + str(threshold) + "x coverage.")

# Should be as above, but instead save each contig as a dictionary key at its first appearance
# Then add one each time that key appears with a value >= the given threshold
# Finally, run the statistics for each contig and print as above

def per_contig():
	print("Obtaining per-contig stats for " + input_bam + "; coverage >+" + str(threshold) + "%.")
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
		print("Length of contig " + key + " = " + str(max_length[key]))
		print(str(value) + "% of contig " + key + " has >=" + str(threshold) + "x coverage.")
		print("---")

def extract_sequence():
	# This function extracts the sequence of the mapped reads 
	# from a part of the reference sequence specified by args.range
	command = ("samtools mpileup -uf %s %s -r %s:%s | bcftools view -cg -") \
	% (args.refference, args.bam, args.chromosome, args.range)
	bam = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
	header = ">" + args.chromosome + ":" + args.range
	seq = ""

	for line in bam.stdout:
		row = line.decode(encoding="utf-8", errors="strict")
		if row[0] == "#":
			pass
		else:
			if row.split("\t")[4] == ".":
				seq += row.split("\t")[3]
			else:
				seq += row.split("\t")[4]
		print(header)
		print(seq)

def main():

	if args.whole == True:
		whole_assembly()
	elif args.contig == True:
		per_contig()


if __name__ == "__main__":
	main()

if args.dev == True:
	print("Time taken = " + str(time.time() - start_time) + " seconds.")
