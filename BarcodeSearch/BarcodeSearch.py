#!/usr/bin/env python3

# Look for ~5,000 bp regions with conserved ends and a variable centre versus the reference
# Then check whether this region is also variable in subsequent samples
# This will allow identification of barcodes flanked by conserved primer sites

## Imports

import argparse
import subprocess
import os
import sys
import io

#######################################################################

## Arguments - define reference, input bams, window size, and primer site size

parser = argparse.ArgumentParser(prog="BarcodeSearch")

parser.add_argument("-f", "--ref", \
			help="Reference")
parser.add_argument("-B", "--BAMs", \
			nargs="+", \
			help="BAM files of samples")
parser.add_argument("--window_size", \
			type=int, \
			default="5000", \
			help="Window size for barcode search")
parser.add_argument("--primer_size", \
			type=int, \
			default="21", \
			help="Desired size of conserved regions at beginning and end of barcode")
parser.add_argument("-t", "--threads", \
			default=1, \
			help="Threads")
parser.add_argument("-c", "--contig", \
			help="Specify contig to investigate")
parser.add_argument("-o", "--out_prefix", \
			help="Prefix for output files")

args = parser.parse_args()

#######################################################################

outFile_windows = str(args.out_prefix + "_windows.txt")

if os.path.isfile(outFile_windows) == True:
	print("Warning: Output file",outFile_windows,"already exists. Please choose another output prefix.")
	sys.exit(0)

print("Searching for potential barcodes in",len(args.BAMs),"file(s).\n")

## ENSURE CORRECT MODULES ARE LOADED - SAMTOOLS + BCFTOOLS

# Record length of contig of interest
cmd2 = ["samtools idxstats %s" % (args.BAMs[0])]
process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell=True)

with process2.stdout as result2:
	rows2 = (line2.decode() for line2 in result2)
	for row2 in rows2:
		if row2.split("\t")[0] == args.contig:
			contig_length = int(row2.split("\t")[1])

#######################################################################

# This ensures that window starts and stops align as intended
primer_size = args.primer_size - 1

all_SNPs = []
all_indels = []

def BarFind(infile,outdict):
	
	# Parse BAM files
	SNP_loci = []
	sample_SNPs = []
	sample_indels = []
	cmd = ["bcftools mpileup --threads %s --fasta-ref %s -r %s %s | bcftools call --threads %s -mv" % \
		(args.threads, args.ref, args.contig, infile, args.threads)]

	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

	# Generate a list of loci where variants occur
	with process.stdout as result:
		rows = (line.decode() for line in result)
		for row in rows:
			if not row.startswith("#"):
				SNP_loci.append(row.split("\t")[1])
				if "INDEL" in row.split("\t")[7]:
					all_indels.append(row.split("\t")[1])
				else:
					all_SNPs.append(row.split("\t")[1])

	# Start stepping through the file

	outdict = {}
	window_coords = []

	# Assign start and stop locations for window and primers
	for window in range(0,(contig_length - args.window_size)):
		window_start = int(window + 1)
		window_stop = int(window_start + args.window_size)
		primer1_stop = int(window_start + args.primer_size)
		primer2_start = int(window_stop - args.primer_size)
		window_coords = (window_start,primer1_stop,primer2_start,window_stop)
		unsuitable = 0

		# If any variants fall within primer sites, skip the window
		for SNP in SNP_loci:
			if window_start <= int(SNP) <= primer1_stop:
				unsuitable += 1
				break
			elif primer2_start <= int(SNP) <= window_stop:
				unsuitable += 1
				break

		# If no variants in primer sites, save the coordinates
		if unsuitable == 0:
			outdict[window_coords[0]] = window_coords

	return(outdict)

#######################################################################

all_files = {}

# Generate dictionary of dictionaries of lists, representing all suitable loci in all samples
for n in range(len(args.BAMs)):
	print(os.path.splitext(os.path.basename(args.BAMs[n]))[0])
	all_files[os.path.splitext(os.path.basename(args.BAMs[n]))[0]] = \
		BarFind(args.BAMs[n],os.path.splitext(os.path.basename(args.BAMs[n]))[0])

master_dict = {}

# If key appears in all sub-dictionaries, save list to master dictionary
	# For every key in subdictionary 1...
for key1 in all_files[os.path.splitext(os.path.basename(args.BAMs[0]))[0]]:
	occurences = 1
	# CHECK NUMBERING BELOW!
	# ... for each sub-dictionary in the rest of the main dictionary...
	for dict1 in range(1,len(args.BAMs)):
	# ... check whether key1 appears...
		if key1 in all_files[os.path.splitext(os.path.basename(args.BAMs[dict1]))[0]]:
			occurences += 1
	# ... and if it appears in all sub-dictionaries...
	if occurences == len(args.BAMs):
	# ... add it to the master dictionary
		master_dict[key1] = all_files[os.path.splitext(os.path.basename(args.BAMs[0]))[0]][key1]

#######################################################################

# Merge overlapping windows

final_dict = {}

saved_window = (0,0,0,0)

for key2 in master_dict:
	if saved_window == (0,0,0,0):
		saved_window = master_dict[key2]
	elif saved_window[0] <= master_dict[key2][0] <= saved_window[1]:
		saved_window = (saved_window[0],master_dict[key2][1],saved_window[2],master_dict[key2][3])
	else:
		final_dict[saved_window[0]] = saved_window
		saved_window = master_dict[key2]

#######################################################################

# Count variants and report results

with open(outFile_windows, "a") as output_file:
	output_file.write("P1_1\tP1_2\tP2_1\tP2_2\tSNPs\tIndels\n")

	for key3 in final_dict:

		window_SNPs = 0
		window_indels = 0

		for SNP in all_SNPs:
			if final_dict[key3][1] < int(SNP) < final_dict[key3][2]:
				window_SNPs += 1
		for indel in all_indels:
			if final_dict[key3][1] < int(indel) < final_dict[key3][2]:
				window_indels += 1

		result = str(final_dict[key3][0]) + "\t" + str(final_dict[key3][1]) + "\t" + \
			str(final_dict[key3][2]) + "\t" + str(final_dict[key3][3]) + "\t" + \
			str(window_SNPs) + "\t" + str(window_indels) + "\n"
		output_file.write(result)

	output_file.close()

#######################################################################
