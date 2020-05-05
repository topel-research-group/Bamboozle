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

import numpy as np
from Levenshtein import distance

#######################################################################
# HOUSEKEEPING
#######################################################################

def barcode(args):

# Ensure the output files doesn't already exist

	out_bed = args.outprefix + ".bed"
	out_txt = args.outprefix + ".txt"

	if os.path.isfile(out_bed) == True:
		print("Warning: Output file",out_bed,"already exists. Please choose another output prefix.")
		sys.exit(0)
	if os.path.isfile(out_txt) == True:
		print("Warning: Output file",out_txt,"already exists. Please choose another output prefix.")
		sys.exit(0)

#######################################################################

# Record all contig lengths

	contig_lengths = {}

	cmd = ["samtools idxstats %s" % (args.sortbam[0])]
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

	with process.stdout as result:
		rows = (line.decode() for line in result)
		for row in rows:
			if row.split("\t")[0] != "*":
				contig_lengths[row.split("\t")[0]] = int(row.split("\t")[1])

#######################################################################

# This ensures that window starts and stops align as intended

#	primer_size = args.primer_size - 1

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

	pre_final_dict = {}
	for contig in contig_lengths:
		pre_final_dict[contig] = {}

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

		filter_qual = "%QUAL>" + str(args.quality)

		cmd2 = ["bcftools mpileup --threads %s --fasta-ref %s %s | bcftools call --threads %s -mv | \
			bcftools filter --threads %s -i '%s'" % \
			(args.threads, args.ref, bam, args.threads, args.threads, filter_qual)]

		process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell=True)

# Generate a list of loci where variants occur
# Also generate (and compress and index) the resultant VCF file, for consensus generation later

		print("\nFinding variants...")

		vcf_file = FileName(bam) + ".vcf"
		with process2.stdout as result2, open(vcf_file, "a") as output_file:
			rows2 = (line.decode() for line in result2)

			current_contig = ""

			for row2 in rows2:
				output_file.write(row2)
				if not row2.startswith("#"):
					contig_name = row2.split("\t")[0]
					variant_position = row2.split("\t")[1]

					if not current_contig == contig_name:
						current_contig = contig_name
						print(current_contig)

					variant_loci[contig_name].append(variant_position)

					if "INDEL" in row2.split("\t")[7]:
						if not contig_name in sample_indels.keys():
							sample_indels[contig_name] = []
						sample_indels[contig_name].append(variant_position)

					else:
						if not contig_name in sample_SNPs.keys():
							sample_SNPs[contig_name] = []
						sample_SNPs[contig_name].append(variant_position)

			output_file.close()
			os.system("bgzip -@ " + str(args.threads) + " " + vcf_file)
			os.system("bcftools index --threads " + str(args.threads) + " " + vcf_file + ".gz")

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
#				primer1_stop = int(window_start + primer_size)
#				primer2_start = int(window_stop - primer_size)
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
						pre_final_dict[contig][saved_window[0]] = \
						[saved_window[0],saved_window[1],saved_window[2],saved_window[3]]

						saved_window = master_dict[contig][window_start]

			pre_final_dict[contig][saved_window[0]] = [saved_window[0],saved_window[1],saved_window[2],saved_window[3]]

#######################################################################
#
# Find a way to skip loci where the following type of error occurs:
# `Warning: ignoring overlapping variant starting at 000215F:2953`
#
#######################################################################

# Compare consensus sequences for all BAMs, to ensure they are truly unique,
# not merely differing from the reference in all the same positions

	print("\nChecking consensuses...")

	for contig in pre_final_dict:
		for window in pre_final_dict[contig]:
			compare_seqs = {}
			for bam in args.sortbam:
				compare_seqs[FileName(bam)] = ""
				vcf_zipped = FileName(bam) + ".vcf.gz"
				cmd3 = ["samtools faidx %s %s:%s-%s | bcftools consensus %s" % \
					(args.ref, contig, pre_final_dict[contig][window][0], pre_final_dict[contig][window][3], vcf_zipped)]
				process3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE, shell=True)
				with process3.stdout as result3:
					rows3 = (line.decode() for line in result3)
					for row3 in rows3:
						compare_seqs[FileName(bam)] += (row3.upper().strip("\n"))

			if not (len(compare_seqs) != len(set(compare_seqs.values()))):	# All values unique
				final_dict[contig][window] = pre_final_dict[contig][window]

# TO DO: calculate the maximum and minimum difference between two sequences

				list1 = list(compare_seqs.values())
				list2 = list(compare_seqs.values())
				matrix = np.zeros((len(list1), len(list2)), dtype=np.int)

				for i in range(0,len(list1)):
					for j in range(0,len(list2)):
						matrix[i,j] = distance(list1[i],list2[j])

# Add minimum and maximum differences to each window
				final_dict[contig][window].append(np.min(matrix[np.nonzero(matrix)]))
				final_dict[contig][window].append(np.max(matrix[np.nonzero(matrix)]))
			

# Define a function to generate the required sequence for conserved and variable regions

	def return_region(start, stop):
		cmd4 = ["samtools faidx %s %s:%s-%s" % (args.ref, contig, start, stop)]
		process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, shell=True)

		return_me = ""

		with process4.stdout as result4:
			rows4 = (line.decode() for line in result4)
			for row4 in rows4:
				if not row4.startswith(">"):
					return_me += row4.strip("\n")

		return(return_me)


# Report results in TXT and BED format

	print("\nGenerating output files...")

	with open(out_bed, "a") as output_bed, open(out_txt, "a") as output_txt:
		output_bed.write("track name=PotentialBarcodes description=\"Potential barcodes\"\n")
		output_txt.write("contig\tconserved_1_start\tconserved_1_end\tconserved_1_seq\tvariable_start\tvariable_end\tvariable_seq\tconserved_2_start\tconserved_2_end\tconserved_2_seq\tvariable_length\tminimum_diffs\tmaximum_diffs\n")

# TO DO: print length of variable region

		for contig in final_dict:

			window_number = 0

			for window in final_dict[contig]:

				window_number += 1
				window_SNPs = 0
				window_indels = 0

				conserved_1_start = final_dict[contig][window][0]
				conserved_1_stop = final_dict[contig][window][1] - 1
				variable_start = final_dict[contig][window][1]
				variable_stop = final_dict[contig][window][2]
				conserved_2_start = final_dict[contig][window][2] + 1
				conserved_2_stop = final_dict[contig][window][3]
				variable_len = conserved_2_start - variable_start
				min_diffs = final_dict[contig][window][4] 
				max_diffs = final_dict[contig][window][5]

				for SNP in all_SNPs[contig]:
					if variable_start < int(SNP) < variable_stop:
						window_SNPs += 1

				for indel in all_indels[contig]:
					if variable_start < int(indel) < variable_stop:
						window_indels += 1

		# Fields 7 and 8 (thickStart and thickEnd) represent the start and stop positions of the non-primer part of the window

				window_out = str(contig) + "\t" + str(conserved_1_start - 1) + "\t" + str(conserved_2_stop) + "\t" + \
						"window_" + str(window_number) + "_SNPs_" + str(window_SNPs) + "_indels_" + str(window_indels) + \
						"\t0\t.\t" + str(conserved_1_stop) + "\t" + str(variable_stop) + "\n"

				line_out = str(contig) + "\t" + str(conserved_1_start) + "\t" + str(conserved_1_stop) + "\t" + return_region(conserved_1_start, conserved_1_stop) + "\t" + \
						str(variable_start) + "\t" + str(variable_stop) + "\t" + return_region(variable_start, variable_stop) + "\t" + \
						str(conserved_2_start) + "\t" + str(conserved_2_stop) + "\t" + return_region(conserved_2_start, conserved_2_stop) + "\t" + \
						str(variable_len) + "\t" + str(min_diffs) + "\t" + str(max_diffs) + "\n"

				output_bed.write(window_out)
				output_txt.write(line_out)

		output_bed.close()
		output_txt.close()
