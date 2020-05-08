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

from time import time
from multiprocessing import Pool

#######################################################################
# GET REASONABLE SAMPLE NAMES
#######################################################################

def FileName(long_name):
	return os.path.splitext(os.path.basename(long_name))[0]

#######################################################################
# RECORD ALL CONTIG LENGTHS
#######################################################################

def get_contig_lengths(firstbam):

	contig_lengths = {}

	cmd = ["samtools idxstats %s" % (firstbam)]
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

	with process.stdout as result:
		rows = (line.decode() for line in result)
		for row in rows:
			if row.split("\t")[0] != "*":
				contig_lengths[row.split("\t")[0]] = int(row.split("\t")[1])

	return(contig_lengths)

#######################################################################
# PARSE BAM FILES USING BCF TOOLS
#######################################################################

def bcf(infile, contig_list, quality, threads, reference):

	print("\n" + FileName(infile))

	if os.path.isfile(FileName(infile) + ".vcf.gz") == True:
		import gzip
		print("VCF file already exists for",FileName(infile),"- reading file...")
		process2 = gzip.open(FileName(infile) + ".vcf.gz", 'rt')

	else:
		cmd2 = ["bcftools mpileup --threads %s --fasta-ref %s %s | bcftools call --threads %s -mv | \
			bcftools filter --threads %s -i '%s'" % \
			(threads, reference, infile, threads, threads, quality)]

		process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell=True)

	return(process2)

#######################################################################
# GENERATE A LIST OF LOCI WHERE VARIANTS OCCUR
#######################################################################

def get_variants(vcf_row, variant_dict, indel_dict, SNP_dict, contig):
	contig_name = vcf_row.split("\t")[0]
	variant_position = vcf_row.split("\t")[1]

	if not contig == contig_name:
		contig = contig_name
		print(contig)

	variant_dict[contig_name].append(variant_position)

	if "INDEL" in vcf_row.split("\t")[7]:
		indel_dict[contig_name].append(variant_position)
	else:
		SNP_dict[contig_name].append(variant_position)
	return(contig)

#######################################################################
# IDENTIFY VARIABLE WINDOWS AND CONSERVED PRIMER SITES
#######################################################################

def find_windows(contig, contig_list, window_len, primer_len, variant_list):
	windows = {}

	for window in range(0,(contig_list[contig] - window_len)):
		window_start = int(window + 1)
		window_stop = int(window_start + window_len)
		primer1_stop = int(window_start + primer_len)
		primer2_start = int(window_stop - primer_len)
		window_coords = [window_start,primer1_stop,primer2_start,window_stop]
	# If any variants fall within primer sites, skip the window
		for variant in variant_list[contig]:
			validity = "True"
			if int(variant) > window_stop:
				break
			elif ((window_start <= int(variant) <= primer1_stop) or (primer2_start <= int(variant) <= window_stop)):
				validity = "False"
				break
	# If no variants in primer sites, save the coordinates
		if validity == "True":
			windows[window_start] = window_coords

	return(windows)

#######################################################################
# MERGE OVERLAPPING WINDOWS
# THIS TAKES A LONG TIME AND SHOULD BE MULTITHREADED
#######################################################################

def merge_windows(contig, window_dict):

	merged = {}	
	saved_window = []

	for window_start in window_dict[contig]:
		if not saved_window:
			saved_window = window_dict[contig][window_start]

		elif saved_window[0] <= window_dict[contig][window_start][0] <= saved_window[1]:
			saved_window = [saved_window[0],window_dict[contig][window_start][1],\
				saved_window[2],window_dict[contig][window_start][3]]

		else:
			merged[saved_window[0]] = [saved_window[0],saved_window[1],saved_window[2],saved_window[3]]
			saved_window = window_dict[contig][window_start]

		merged[saved_window[0]] = [saved_window[0],saved_window[1],saved_window[2],saved_window[3]]

	return(merged)

#######################################################################
# ENSURE UNIQUENESS OF VARIABLE REGIONS BETWEEN STRAINS
# THIS TAKES A VERY LONG TIME AND SHOULD BE MULTITHREADED
#######################################################################

def check_unique_windows(windows, contig, reference, infiles):
	final = {}
	for window in windows[contig]:
		compare_seqs = {}

		for bam in infiles:
			compare_seqs[FileName(bam)] = ""
			vcf_zipped = FileName(bam) + ".vcf.gz"
			cmd3 = ["samtools faidx %s %s:%s-%s | bcftools consensus %s" % \
				(reference, contig, windows[contig][window][0], windows[contig][window][3], vcf_zipped)]
			process3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE, shell=True)
			with process3.stdout as result3:
				rows3 = (line.decode() for line in result3)
				for row3 in rows3:
					compare_seqs[FileName(bam)] += (row3.upper().strip("\n"))

# Ensure all variable regions are unique between samples

		if not (len(compare_seqs) != len(set(compare_seqs.values()))):
			final[window] = windows[contig][window]

# Calculate the number of differences between each pair of sequences

			list1 = list(compare_seqs.values())
			list2 = list(compare_seqs.values())
			matrix = np.zeros((len(list1), len(list2)), dtype=np.int)

			for i in range(0,len(list1)):
				for j in range(0,len(list2)):
					matrix[i,j] = distance(list1[i],list2[j])

# Add minimum and maximum differences to each window
			final[window].append(np.min(matrix[np.nonzero(matrix)]))
			final[window].append(np.max(matrix[np.nonzero(matrix)]))

	return(final)


#######################################################################
# GENERATE THE REQUIRED SEQUENCE FOR CONSERVED AND VARIABLE REGIONS
#######################################################################

def return_region(start, stop, reference, contig):

	return_me = ""

	cmd4 = ["samtools faidx %s %s:%s-%s" % (reference, contig, start, stop)]
	process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, shell=True)

	with process4.stdout as result4:
		rows4 = (line.decode() for line in result4)
		for row4 in rows4:
			if not row4.startswith(">"):
				return_me += row4.strip("\n")

	return(return_me)

#######################################################################
# HOUSEKEEPING
#######################################################################

def barcode(args):

	if args.dev == True:
		start_time = time()

# Set number of threads

	pool = Pool(processes = int(args.threads))

# Ensure the output files doesn't already exist

	out_bed = args.outprefix + ".bed"
	out_txt = args.outprefix + ".txt"

	if os.path.isfile(out_bed) == True:
		print("Warning: Output file",out_bed,"already exists. Please choose another output prefix.")
		sys.exit(0)
	if os.path.isfile(out_txt) == True:
		print("Warning: Output file",out_txt,"already exists. Please choose another output prefix.")
		sys.exit(0)

# Record all contig lengths

	contig_lengths = get_contig_lengths(args.sortbam[0])

# Set initial global lists/dictionaries

	all_SNPs = {}
	for contig in contig_lengths:
		all_SNPs[contig] = []

	all_indels = {}
	for contig in contig_lengths:
		all_indels[contig] = []

	all_variants = {}
	for contig in contig_lengths:
		all_variants[contig] = []

	all_files = {}
	for contig in contig_lengths:
		all_files[contig] = []

	master_dict = {}
	for contig in contig_lengths:
		master_dict[contig] = {}

	merged_dict = {}
	for contig in contig_lengths:
		merged_dict[contig] = {}

	final_dict = {}
	for contig in contig_lengths:
		final_dict[contig] = {}

# Set other values
	
	filter_qual = "%QUAL>" + str(args.quality)

	if args.dev == True:
		print("Housekeeping =",(time() - start_time),"seconds.")
		start_time = time()

#######################################################################
# MAIN CODE
#######################################################################

	print("Searching for potential barcodes in",len(args.sortbam),"file(s).")

	for bam in args.sortbam:

# Generate or read in a VCF file for the current BAM

		process2 = bcf(bam, contig_lengths, filter_qual, args.threads, args.ref)
		if args.dev == True:
			print("Parsing",bam,"=",(time() - start_time),"seconds.")
			start_time = time()

# Generate a list of loci where variants occur
# ADJUST TO MAKE THIS PARALLELISABLE!

		print("\nFinding variants...")

		current_contig = ""

		if os.path.isfile(FileName(bam) + ".vcf.gz") == True:
			for row2 in process2:
				if not row2.startswith("#"):
					current_contig = get_variants(row2, all_variants, all_indels, all_SNPs, current_contig)

		else:

			vcf_file = FileName(bam) + ".vcf"
			with process2.stdout as result2, open(vcf_file, "a") as output_file:
				rows2 = (line.decode() for line in result2)

				for row2 in rows2:
					output_file.write(row2)
					if not row2.startswith("#"):
						current_contig = get_variants(row2, all_variants, all_indels, all_SNPs, current_contig)

			output_file.close()

# If the VCF file didn't exist previously, then bgzip the newly-written one, and index it

			os.system("bgzip -@ " + str(args.threads) + " " + vcf_file)
		if os.path.isfile(FileName(bam) + ".vcf.gz.csi") == False:
			os.system("bcftools index --threads " + str(args.threads) + " " + vcf_file + ".gz")

		if args.dev == True:
			print("Finding variants in",bam,"=",(time() - start_time),"seconds.")
			start_time = time()

# Get sorted lists of variant/SNP/indel positions, with duplicates removed

	for contig in all_variants:
		all_variants[contig] = sorted(list(set(all_variants[contig])), key=int)

	for contig in all_SNPs:
		all_SNPs[contig] = sorted(list(set(all_SNPs[contig])), key=int)

	for contig in all_indels:
		all_indels[contig] = sorted(list(set(all_indels[contig])), key=int)

# PARALLELISE THIS STEP! Step through each contig, assigning start and stop locations for window and primers

	print("\nChecking windows...")

	to_master = pool.starmap(find_windows, \
	[(contig, contig_lengths, args.window_size, args.primer_size, all_variants) for contig in contig_lengths])

	for entry in range(0,len(contig_lengths)):
		master_dict[list(contig_lengths.keys())[entry]] = to_master[entry]

	if args.dev == True:
		print("Finding windows =",(time() - start_time),"seconds.")
		start_time = time()

#######################################################################

# Merge overlapping windows
# THIS TAKES A LONG TIME AND SHOULD BE MULTITHREADED

	print("\nMerging overlapping windows...")

	for contig in master_dict:
		if master_dict[contig]:
			print(contig)
			merged_dict[contig] = merge_windows(contig, master_dict)

	if args.dev == True:
		print("Merging windows =",(time() - start_time),"seconds.")
		start_time = time()

#######################################################################
#
# Find a way to skip loci where the following type of error occurs:
# `Warning: ignoring overlapping variant starting at 000215F:2953`
#
#######################################################################

# Compare consensus sequences for all BAMs, to ensure they are truly unique,
# not merely differing from the reference in all the same positions

# THIS TAKES A VERY LONG TIME AND SHOULD BE MULTITHREADED

	print("\nChecking consensuses...")

	for contig in merged_dict:
		final_dict[contig] = check_unique_windows(merged_dict, contig, args.ref, args.sortbam)

	if args.dev == True:
		print("Comparing consensus sequences =",(time() - start_time),"seconds.")
		start_time = time()

# Report results in TXT and BED format

	print("\nGenerating output files...")

	with open(out_bed, "a") as output_bed, open(out_txt, "a") as output_txt:
		output_bed.write("track name=PotentialBarcodes description=\"Potential barcodes\"\n")
		output_txt.write("contig\tconserved_1_start\tconserved_1_end\tconserved_1_seq\tvariable_start\tvariable_end\tvariable_seq\tconserved_2_start\tconserved_2_end\tconserved_2_seq\tvariable_length\tminimum_diffs\tmaximum_diffs\n")

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

				line_out = str(contig) + "\t" + str(conserved_1_start) + "\t" + str(conserved_1_stop) + "\t" + return_region(conserved_1_start, conserved_1_stop, args.ref, contig) + "\t" + \
						str(variable_start) + "\t" + str(variable_stop) + "\t" + return_region(variable_start, variable_stop, args.ref, contig) + "\t" + \
						str(conserved_2_start) + "\t" + str(conserved_2_stop) + "\t" + return_region(conserved_2_start, conserved_2_stop, args.ref, contig) + "\t" + \
						str(variable_len) + "\t" + str(min_diffs) + "\t" + str(max_diffs) + "\n"

				output_bed.write(window_out)
				output_txt.write(line_out)

		output_bed.close()
		output_txt.close()

		if args.dev == True:
			print("Writing output =",(time() - start_time),"seconds.")
			start_time = time()
