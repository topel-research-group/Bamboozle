#!/usr/bin/env python3

# With a gff and the output of coverage_stats, identify which genes get mapped to
# Intended for use with the S. subsalsum data
#
# Will require altering coverage_stats

#######################################################################
# IMPORTS
#######################################################################

import argparse
import subprocess
import re
import os
from time import time

#######################################################################
# ARGUMENTS
#######################################################################

parser = argparse.ArgumentParser(prog="Bamboozle")

parser.add_argument("--sortbam", \
			nargs='*', \
			help="Sorted BAM infile (N.B. only the BarcodeSearch function accepts multiple inputs)")
parser.add_argument('-c', '--contig', \
			help='Gives per-contig coverage stats')
parser.add_argument('-d', '--threshold', \
			type=int, \
			nargs='?', \
			const='1', \
			default='20', \
			help='Threshold for calculating the coverage percentage; default 20')
parser.add_argument("--gff", \
			help="gff infile")
parser.add_argument("-o", "--outfile", \
			help="Output filename")
parser.add_argument("-v", "--verbose", \
			action="store_true", \
			help="Be more verbose")

parser.add_argument('--dev', \
			help=argparse.SUPPRESS, action="store_true")

args = parser.parse_args()

#######################################################################
# HANDLING OF MULTIPLE SORTED BAM INPUTS
#	Ensure that if multiple sorted BAM inputs are specified,
#	BarcodeSearch is the function being run
#	Else warn the user and exit
# CAN BE REMOVED WHEN INTEGRATED INTO BAMBOOZLE
#######################################################################

if args.sortbam:
	if len(args.sortbam) == 1:
		args.sortbam = args.sortbam[0]

	elif len(args.sortbam) > 1 and not args.barcode:
		print("Please note that only BarcodeSearch currently accepts multiple BAM inputs.")
		exit()

#######################################################################
# TESTING RUNTIME (--dev)
# CAN BE REMOVED WHEN INTEGRATED INTO BAMBOOZLE
#######################################################################

if args.dev == True:
	start_time = time()

#######################################################################

def coverage_stats_2():

# Generate a dictionary of dictionaries of lists; contig -> gene name -> coordinates

	if args.gff:
		in_gff_loci = {}

		with open(args.gff, 'r') as input:
			lines = (entry.split("\t") for entry in input)
			for line in lines:

				# Is this standard in gff files?
				if "FASTA" in line[0]:
					break

				elif not line[0].startswith("#"):
					contig_name = line[0]
					feature_type = line[2]
					feature_start = line[3]
					feature_stop = line[4]
					feature_name = re.split("=|;",line[8])[1]

					if not contig_name in in_gff_loci.keys():
						in_gff_loci[contig_name] = {}

					if feature_type == "gene":
						in_gff_loci[contig_name][feature_name] = (feature_start,feature_stop)

#######################################################################

# Record all contig lengths (may be overkill for single-contig analyses)
# Taken from `barcodesearch.py`

	contig_lengths = {}

	cmd2 = ["samtools idxstats %s" % (args.sortbam)]
	process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell=True)

	with process2.stdout as result2:
		rows2 = (line2.decode() for line2 in result2)
		for row2 in rows2:
			if row2.split("\t")[0] != "*":
				contig_lengths[row2.split("\t")[0]] = int(row2.split("\t")[1])

	assembly_length = sum(contig_lengths.values())

#######################################################################

	if args.contig:
		cmd = ["samtools depth -aa %s -r %s" % (args.sortbam, args.contig)]
		sequence = "contig"
	else:
		cmd = ["samtools depth -aa %s" % args.sortbam]
		sequence = "assembly"

	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

	if args.verbose == True:
		if args.contig:
			print("Obtaining stats for ",args.contig," in ",os.path.basename(args.sortbam),"; coverage >=",args.threshold,"x.",sep="")
		else:
			print("Obtaining whole-genome stats for ",os.path.basename(args.sortbam),"; coverage >=",args.threshold,"x.",sep="")

	with open(args.outfile, "a") as output_file:
		output_file.write("track name=Coverage description=\"Coverage above " + str(args.threshold) + "x\"\n")

	# Define starting variables

		# Start and stop of region of interest
		start_coord = 0
		stop_coord = 0
		current_contig = ""
		recording = False

		# For calculating percentages (can be refined)
		cov_stats = {}
		for contig in contig_lengths.keys():
			cov_stats[contig] = 0

# USE CODE FROM BAMPARSER'S COVERAGE_LIMITS AS GUIDELINE FOR SAVING LOCI

		with process.stdout as result:
			rows = (line.decode().split('\t') for line in result)
			for row in rows:
				contig = row[0]
				position = int(row[1])
				coverage = int(row[2])

				if coverage >= args.threshold:
					cov_stats[contig] += 1
					if not recording:
						start_coord = position
						recording = True
				elif coverage < args.threshold and recording == True:
					stop_coord = position - 1
					output_line = contig + "\t" + str(start_coord - 1) + "\t" + str(stop_coord) + "\n"
					output_file.write(output_line)
					recording = False

# Start saving start_coord and stop_coord, cross-referencing with gff file if given, and writing to output_file

				# ...

				# if args.gff:
					# output_file.write(contig + "\t" + (start_coord - 1) + stop_coord + ...)
				# else:
					# output_file.write(contig + "\t" + (start_coord - 1) + stop_coord)

			if args.contig:
				print("Length of " + args.contig + ":",contig_lengths[args.contig])
				value = 100.0 / contig_lengths[args.contig] * cov_stats[args.contig]
				print(round(value, 3),"% of contig ", args.contig, " with >=",args.threshold,"x coverage.",sep="")
			else:
				print("Length of assembly:",str(assembly_length))
				value = 100.0 / assembly_length * sum(cov_stats.values())
				print(round(value, 3),"% of assembly with >=",args.threshold,"x coverage.",sep="")

		output_file.close()

coverage_stats_2()


# CONTIG	START	STOP	LABEL

#######################################################################
# CAN BE REMOVED WHEN INTEGRATED INTO BAMBOOZLE
#######################################################################

if args.dev == True:
	print("Time taken =",(time() - start_time),"seconds.")
