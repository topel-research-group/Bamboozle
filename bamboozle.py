#!/usr/bin/env python3

import sys
import subprocess
import argparse
import time

parser = argparse.ArgumentParser(description='Obtain statistics regarding percentage coverage from bam files. \
                                              The script gives percentage of positions in an assembly/contig \
                                              with coverage greater than or equal to a given threshold')
parser.add_argument('-c', '--contig', help='Gives per-contig coverage stats')
#parser.add_argument('-c', '--contig', type=str, nargs='+', help='Gives cov. stats for the specified contigs')
parser.add_argument('-t' ,'--threshold', type=int, nargs='?', const=1, default='20', help='Threshold for calculating coverage percentage; default 20')
parser.add_argument("-r", "--refference", help="Reference sequence file")
parser.add_argument("-b", "--bam", help="Bam file")
parser.add_argument("--range", help="somethingsomsing")
parser.add_argument("-z", "--zero", action="store_true", help="Find regions of 0x coverage")
parser.add_argument("-v", "--verbose", action="store_true", help="Be more verbose")
parser.add_argument('--dev', help=argparse.SUPPRESS, action="store_true")
args = parser.parse_args()


if args.dev == True:
	start_time = time.time()


def coverage_stats():
	# This function calculates the percentage of positions in an assembly/contig
	# with read coverage >= a given threshold (default: 20x)

	try:
		subprocess.check_output('samtools depth 2>&1 | grep -- "-aa"', stderr=subprocess.PIPE, shell=True)
	except subprocess.CalledProcessError:
		print("This version of samtools does not support the `depth -aa` option; please update samtools.")
		exit()

	cmd = ["samtools depth -aa %s" % args.bam]
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

	if args.verbose == True:
		if args.contig:
			print("Obtaining stats for " + args.contig + " in " + args.bam + "; coverage >+" + str(args.threshold) + "%.")
		else:
			print("Obtaining whole-genome stats for " + args.bam + "; coverage >+" + str(args.threshold) + "%.")
	cov_stats = {}
	num_lines = 0
	with process.stdout as result:
		rows = (line.decode().split('\t') for line in result)
		for row in rows:
			coverage = int(row[2])
			if args.contig:
				if str(row[0]) == args.contig:
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
	print(str(value) + "% of the assembly/contig has >=" + str(args.threshold) + "x coverage.")



def zero_regions():
	# This function identifies regions of 0x coverage, then prints the reference sequence
	# and GC content at these coordinates

	print("Blah")



def extract_sequence():
	# This function extracts the sequence of the mapped reads 
	# from a part of the reference sequence specified by args.range
	command = ("samtools mpileup -uf %s %s -r %s:%s | bcftools view -cg -") \
	% (args.refference, args.bam, args.contig, args.range)
	bam = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
	header = ">" + args.contig + ":" + args.range
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
	if args.zero:
		zero_regions()
	elif args.range:
		extract_sequence()
	else:
		coverage_stats()


if __name__ == "__main__":
	main()

if args.dev == True:
	print("Time taken = " + str(time.time() - start_time) + " seconds.")
