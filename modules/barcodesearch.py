#!/usr/bin/env python3

#	Look for genomic regions with conserved ends and a variable centre versus the reference,
#	then check whether these regions are also variable in subsequent samples.
#	This will allow identification of barcode regions flanked by conserved primer sites
#
#	Copyright (C) 2019 Matthew Pinder. matt_pinder13@hotmail.com
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <https://www.gnu.org/licenses/>.

#######################################################################
# IMPORTS
#######################################################################

import subprocess
import os
import sys
import io

#######################################################################
# HOUSEKEEPING
#######################################################################

def barcode(args):

# Ensure the output BED file doesn't already exist

	if os.path.isfile(args.outfile) == True:
		print("Warning: Output file",args.outfile,"already exists. Please choose another output prefix.")
		sys.exit(0)

#######################################################################

# Record all contig lengths

	contig_lengths = {}

	cmd2 = ["samtools idxstats %s" % (args.sortbam[0])]
	process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell=True)

	with process2.stdout as result2:
		rows2 = (line2.decode() for line2 in result2)
		for row2 in rows2:
			if row2.split("\t")[0] != "*":
				contig_lengths[row2.split("\t")[0]] = int(row2.split("\t")[1])

#######################################################################

# This ensures that window starts and stops align as intended

	primer_size = args.primer_size - 1

#######################################################################

# Get a reasonable name for each sample

	def FileName(long_name):
		return(os.path.splitext(os.path.basename(long_name))[0])

#######################################################################

# Set initial global lists/dictionaries

	all_SNPs = {}
	for contig in contig_lengths:
		all_SNPs[contig] = []

	all_indels = {}
	for contig in contig_lengths:
		all_indels[contig] = []

	all_files = {}
	for contig in contig_lengths:
		all_files[contig] = []

	master_dict = {}
	for contig in contig_lengths:
		master_dict[contig] = {}

	final_dict = {}
	for contig in contig_lengths:
		final_dict[contig] = {}

#######################################################################
# MAIN CODE
#######################################################################

	print("Searching for potential barcodes in",len(args.sortbam),"file(s).")

	file_number = 0

	for bam in args.sortbam:

		print("\n" + FileName(bam))

# Parse BAM files
		variant_loci = {}
		for contig in contig_lengths:
			variant_loci[contig] = []

		sample_SNPs = {}
		sample_indels = {}

		cmd = ["bcftools mpileup --threads %s --fasta-ref %s %s | bcftools call --threads %s -mv" % \
			(args.threads, args.ref, bam, args.threads)]

		process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

# Generate a list of loci where variants occur

		print("\nFinding variants...")

		with process.stdout as result:
			rows = (line.decode() for line in result)

			current_contig = ""

			for row in rows:
				if not row.startswith("#"):
					contig_name = row.split("\t")[0]
					variant_position = row.split("\t")[1]

					if not current_contig == contig_name:
						current_contig = contig_name
						print(current_contig)

					variant_loci[contig_name].append(variant_position)

					if "INDEL" in row.split("\t")[7]:
						if not contig_name in sample_indels.keys():
							sample_indels[contig_name] = []
						sample_indels[contig_name].append(variant_position)

					else:
						if not contig_name in sample_SNPs.keys():
							sample_SNPs[contig_name] = []
						sample_SNPs[contig_name].append(variant_position)

		print("\nAdding variants to list...")
		for contig in sample_SNPs:
			for SNP in sample_SNPs[contig]:
				if not SNP in all_SNPs[contig]:
					all_SNPs[contig].append(SNP)

		for contig in sample_indels:
			for indel in sample_indels[contig]:
				if not indel in all_indels[contig]:
					all_indels[contig].append(indel)

# Start stepping through the file

# Assign start and stop locations for window and primers

		print("\nChecking windows...")

		for contig in contig_lengths:

			print(contig)
			for window in range(0,(contig_lengths[contig] - args.window_size)):
				window_start = int(window + 1)
				window_stop = int(window_start + args.window_size)
				primer1_stop = int(window_start + args.primer_size)
				primer2_start = int(window_stop - args.primer_size)
				window_coords = [window_start,primer1_stop,primer2_start,window_stop,1]

				if bam == args.sortbam[0]:
				# If any variants fall within primer sites, skip the window
					for variant in variant_loci[contig]:
						validity = "true"
						if int(variant) > window_stop:
							break
						elif ((window_start <= int(variant) <= primer1_stop) or (primer2_start <= int(variant) <= window_stop)):
							validity = "false"
							break
				# If no variants in primer sites, save the coordinates
					if validity == "true":
						master_dict[contig][window_start] = window_coords

				elif window_start in master_dict[contig].keys() and master_dict[contig][window_start][4] == file_number:
				# If any variants fall within primer sites, skip the window
					for variant in variant_loci[contig]:
						validity = "true"
						if int(variant) > window_stop:
							break
						elif ((window_start <= int(variant) <= primer1_stop) or (primer2_start <= int(variant) <= window_stop)):
							validity = "false"
							break
				# If no variants in primer sites, increase instances by 1
					if validity == "true":
						master_dict[contig][window_start] = \
						[window_start,primer1_stop,primer2_start,window_stop,(master_dict[contig][window_start][4] + 1)]

		file_number += 1

#######################################################################

# Merge overlapping windows

	print("\nMerging overlapping windows...")

	for contig in master_dict:
		if master_dict[contig]:

			saved_window = []

			print(contig)
			for window_start in master_dict[contig]:
				if master_dict[contig][window_start][4] == len(args.sortbam):
					if not saved_window:
						saved_window = master_dict[contig][window_start]

					elif saved_window[0] <= master_dict[contig][window_start][0] <= saved_window[1]:
						saved_window = [saved_window[0],master_dict[contig][window_start][1],\
							saved_window[2],master_dict[contig][window_start][3]]

					else:
						final_dict[contig][saved_window[0]] = \
						[saved_window[0],saved_window[1],saved_window[2],saved_window[3]]

						saved_window = master_dict[contig][window_start]

			final_dict[contig][saved_window[0]] = [saved_window[0],saved_window[1],saved_window[2],saved_window[3]]

#######################################################################

# Report results in BED format

	print("\nGenerating BED file...")

	with open(args.outfile, "a") as output_file:
		output_file.write("track name=PotentialBarcodes description=\"Potential barcodes\"\n")

		for contig in final_dict:

			window_number = 0

			for window in final_dict[contig]:

				window_number += 1
				window_SNPs = 0
				window_indels = 0

				for SNP in all_SNPs[contig]:
					if final_dict[contig][window][1] < int(SNP) < final_dict[contig][window][2]:
						window_SNPs += 1

				for indel in all_indels[contig]:
					if final_dict[contig][window][1] < int(indel) < final_dict[contig][window][2]:
						window_indels += 1

			# Fields 7 and 8 (thickStart and thickEnd) represent the start and stop positions of the non-primer part of the window

				window_out = str(contig) + "\t" + str(final_dict[contig][window][0] - 1) + "\t" + str(final_dict[contig][window][3]) + "\t" + \
						"window_" + str(window_number) + "_SNPs_" + str(window_SNPs) + "_indels_" + str(window_indels) + \
						"\t0\t.\t" + str(final_dict[contig][window][1]) + "\t" + str(final_dict[contig][window][2] - 1) + "\n"

				output_file.write(window_out)

		output_file.close()
