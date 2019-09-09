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

from time import time
import datetime

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
parser.add_argument("-o", "--outfile", \
			help="Output filename")

parser.add_argument("--dev", \
			help=argparse.SUPPRESS, action="store_true")

args = parser.parse_args()

#######################################################################

# For determining how long functions take to run
if args.dev == True:
	start_time = time()

#######################################################################
# HOUSEKEEPING
#######################################################################

# Ensure the output BED file doesn't already exist

if os.path.isfile(args.outfile) == True:
	print("Warning: Output file",args.outfile,"already exists. Please choose another output prefix.")
	sys.exit(0)

#######################################################################

# Record length of contig of interest

contig_lengths = {}

cmd2 = ["samtools idxstats %s" % (args.BAMs[0])]
process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell=True)

with process2.stdout as result2:
	rows2 = (line2.decode() for line2 in result2)
	for row2 in rows2:
		if row2.split("\t")[0] != "*":
			contig_lengths[row2.split("\t")[0]] = row2.split("\t")[1]

#######################################################################

# This ensures that window starts and stops align as intended

primer_size = args.primer_size - 1

#######################################################################

# Get a reasonable name for each sample

def FileName(long_name):
	return(os.path.splitext(os.path.basename(long_name))[0])

#######################################################################

## Define a class to use for SNPs and indels
#
#class Variant:
#	def __init__(self, contig, position):
#		self.contig = contig
#		self.position = position
#	def getVariant(self):
#		return self.contig, self.position
#	def setContig(self, contig):
#		self.contig = contig
#	def setPosition(self, position):
#		self.position = position
#	def __eq__(self, other):
#		return (self.contig, self.position) == (other.contig, other.position)

#######################################################################

# Define a class to use for windows

class Window:
	def __init__(self, contig, win_start, p1_stop, p2_start, win_stop, validity):
		self.contig = contig
		self.win_start = win_start
		self.p1_stop = p1_stop
		self.p2_start = p2_start
		self.win_stop = win_stop
		self.validity = validity
	def getWindow(self):
		return self.contig, self.win_start, self.p1_stop, self.p2_start, self.win_stop, self.validity
	def setWin_start(self, win_start):
		self.win_start = win_start
	def setP1_stop(self, p1_stop):
		self.p1_stop = p1_stop
	def setP2_start(self, p2_start):
		self.p2_start = p2_start
	def setWin_stop(self, win_stop):
		self.win_stop = win_stop
	def setValidity(self, validity):
		self.validity = validity
	def __eq__(self, other):
		return (self.contig, self.win_start, self.p1_stop, self.p2_start, self.win_stop, self.validity) \
			== (other.contig, other.win_start, other.p1_stop, other.p2_start, other.win_stop, other.validity)

#######################################################################

# Set initial global lists/dictionaries

all_SNPs = {}
for contig in contig_lengths:
	all_SNPs[contig] = []

all_indels = {}
for contig in contig_lengths:
	all_indels[contig] = []

all_files = {}

master_list = []
final_list = []

#######################################################################
# MAIN CODE
#######################################################################

print("Searching for potential barcodes in",len(args.BAMs),"file(s).\n")


for bam in args.BAMs:

	print(FileName(bam))

	# Parse BAM files
	variant_loci = {}
	sample_SNPs = {}
	sample_indels = {}

	cmd = ["bcftools mpileup --threads %s --fasta-ref %s %s | bcftools call --threads %s -mv" % \
		(args.threads, args.ref, bam, args.threads)]

	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

	# Generate a list of loci where variants occur

	print("Finding variants...")

	with process.stdout as result:
		rows = (line.decode() for line in result)
		for row in rows:
			if not row.startswith("#"):
				contig_name = row.split("\t")[0]
				variant_position = row.split("\t")[1]

				if not contig_name in variant_loci.keys():
					print(contig_name)
					variant_loci[contig_name] = []
				variant_loci[contig_name].append(variant_position)

				if "INDEL" in row.split("\t")[7]:
					if not contig_name in sample_indels.keys():
						sample_indels[contig_name] = []
					sample_indels[contig_name].append(variant_position)

				else:
					if not contig_name in sample_SNPs.keys():
						sample_SNPs[contig_name] = []
					sample_SNPs[contig_name].append(variant_position)

	print("Adding variants to list...")
	for contig in sample_SNPs:
		for SNP in sample_SNPs[contig]:
			if not SNP in all_SNPs[contig]:
				all_SNPs[contig].append(SNP)

	for contig in sample_indels:
		for indel in sample_indels[contig]:
			if not indel in all_indels[contig]:
				all_indels[contig].append(indel)

	# Start stepping through the file

	outlist = []

	# Assign start and stop locations for window and primers

	print("Checking windows...")

	for contig in contig_lengths:

		print(contig)
		for window in range(0,(int(contig_lengths[contig]) - args.window_size)):
			window_start = int(window + 1)
			window_stop = int(window_start + args.window_size)
			primer1_stop = int(window_start + args.primer_size)
			primer2_start = int(window_stop - args.primer_size)
			window_coords = Window(contig,window_start,primer1_stop,primer2_start,window_stop,"true")

			# If any variants fall within primer sites, skip the window
			for variant in variant_loci[contig]:
				if (window_start <= int(variant) <= primer1_stop) or \
				(primer2_start <= int(variant) <= window_stop):
					window_coords.setValidity("false")
					break

			# If no variants in primer sites, save the coordinates
			if window_coords.validity == "true":
				outlist.append(window_coords)

	all_files[FileName(bam)] = outlist

#######################################################################
# 	CHECK THIS BLOCK AND DOWN!
#######################################################################

# If key appears in all sub-dictionaries, save list to master dictionary

# For every valid window identified in the first BAM file parsed...

print("Checking valid windows for presence in every sample...")
for window in all_files[FileName(args.BAMs[0])]:

	occurences = 1

	# CHECK NUMBERING BELOW!

		# ... for each sub-dictionary in the rest of the main dictionary...
	for bam in range(1,len(args.BAMs)):

		# ... check whether the contig/window pair appears...
		if window in all_files[FileName(args.BAMs[bam])]:
			occurences += 1

	# ... and if it appears in all sub-dictionaries...
	if occurences == len(args.BAMs):

		# ... add it to the master dictionary
		print(window)
		master_list.append(window)

for window in master_list:
	print(window.getWindow())

#######################################################################
#
## Merge overlapping windows
#
#saved_window = []
#
#for key2 in master_list:
#
#	if not saved_window:
#		saved_window = master_list[key2]
#
#	elif saved_window[0] <= master_list[key2][0] <= saved_window[1]:
#		saved_window = [saved_window[0],master_list[key2][1],saved_window[2],master_list[key2][3]]
#
#	else:
#		final_list[saved_window[0]] = saved_window
#		saved_window = master_list[key2]
#
#######################################################################
#
# Report results in BED format
#
#with open(args.outfile, "a") as output_file:
#	output_file.write("track name=PotentialBarcodes description=\"Potential barcodes\"\n")
#
#	for contig in final_list:
#
#		window_number = 0
#
#		for key3 in final_list[contig]:
#
#			window_number += 1
#			window_SNPs = 0
#			window_indels = 0
#
#			for SNP in all_SNPs[contig]:
#				if final_list[key3][1] < int(SNP) < final_list[key3][2]:
#					window_SNPs += 1
#
#			for indel in all_indels:
#				if final_list[key3][1] < int(indel) < final_list[key3][2]:
#					window_indels += 1
#
#			# Fields 7 and 8 (thickStart and thickEnd) represent the start and stop positions of the non-primer part of the window
#
#			window = str(contig) + "\t" + str(final_list[key3][0] - 1) + "\t" + str(final_list[key3][3]) + "\t" + \
#					"window_" + str(window_number) + "_SNPs_" + str(window_SNPs) + "_indels_" + str(window_indels) + \
#					"\t0\t.\t" + str(final_list[key3][1]) + "\t" + str(final_list[key3][2] - 1) + "\n"
#
#			output_file.write(window)
#
#	output_file.close()
#
#######################################################################

if args.dev == True:
	print("Time taken =",(time() - start_time),"seconds.")
