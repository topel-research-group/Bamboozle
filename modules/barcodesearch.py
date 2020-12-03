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
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program. If not, see <https://www.gnu.org/licenses/>.

#######################################################################
# IMPORTS
#######################################################################

import subprocess
import os
import sys
import io
from itertools import combinations
from Levenshtein import distance
from multiprocessing import Pool
from time import time

from collections import Counter
from statistics import median
import vcfpy

#######################################################################
# GET REASONABLE SAMPLE NAMES
#######################################################################

def FileName(long_name):
	return os.path.splitext(os.path.basename(long_name))[0]

def NoExt(long_name):
	return os.path.splitext(long_name)[0]

#######################################################################
# RECORD ALL CONTIG LENGTHS BY CHECKING FIRST BAM FILE
#######################################################################

def get_contig_lengths(firstbam):
	contig_lengths = {}

	cmd = ["samtools", "idxstats", firstbam]
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=False)

	with process.stdout as result:
		rows = (line.decode() for line in result)
		for row in rows:
			if row.split("\t")[0] != "*":
				contig_lengths[row.split("\t")[0]] = int(row.split("\t")[1])
	return(contig_lengths)

#######################################################################
# RECORD ALL CONTIGS' AVERAGE COVERAGE PER BAM FILE
#######################################################################

def get_median(bamfile, contigs):
	
	coverage_stats = {}
	for contig in contigs.keys():
		coverage_stats[contig] = 0

	cmd = ["samtools", "depth", "-aa", bamfile]
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=False)
	current_contig = "None"
	this_contig = {}
	with process.stdout as result:
		rows = (line.decode().split('\t') for line in result)
		for row in rows:
			ctg = str(row[0])
			position = int(row[1])
			coverage = int(row[2])
			if current_contig == "None":
				current_contig = ctg
			this_contig[position] = coverage
			if position == contigs[ctg]:
				coverage_stats[current_contig] = median(this_contig.values())
				this_contig = {}
				current_contig = "None"
	return(coverage_stats)

#######################################################################
# IDENTIFY AREAS OF LOW COVERAGE FOR EACH SAMPLE
#######################################################################

def get_lowcov(bam, contigs, coverage_stats):
	low_cov = {}
	for contig in contigs:
		low_cov[contig] = {}

	cmdX = ["bedtools", "genomecov", "-bga", "-ibam", bam]
	procX = subprocess.Popen(cmdX, stdout=subprocess.PIPE, shell=False)

	with procX.stdout as result:
		rows = (line.decode().split("\t") for line in result)
		for row in rows:
			ctg = str(row[0])
			start_pos = int(row[1]) + 1
			end_pos = int(row[2])
			coverage = int(row[3])
			if (coverage > coverage_stats[bam][ctg]*2) or (coverage < coverage_stats[bam][ctg]*0.5):
				low_cov[ctg][start_pos] = [start_pos, end_pos, coverage]
	return(low_cov)

#######################################################################
# PARSE BAM FILES USING BCFTOOLS
#	DevNote - adjust so that the VCF saves in the BAM's directory
#######################################################################

def bcf(infile, contig_list, quality, threads, reference):
	print("\n" + FileName(infile))

	if os.path.isfile(NoExt(infile) + ".vcf.gz") == True:
		import gzip
		print("VCF file already exists for",FileName(infile),"- reading file...")
		process2 = gzip.open(NoExt(infile) + ".vcf.gz", 'rt')
	else:
		cmdA = ["bcftools", "mpileup", "--threads", threads, "--fasta-ref", reference, infile]
		procA = subprocess.Popen(cmdA, stdout=subprocess.PIPE, shell=False)

		cmdB = ["bcftools", "call", "--threads", threads, "-mv"]
		procB = subprocess.Popen(cmdB, stdin=procA.stdout, stdout=subprocess.PIPE, shell=False)

		cmdC = ["bcftools", "filter", "--threads", threads, "-i", quality]
		process2 = subprocess.Popen(cmdC, stdin=procB.stdout, stdout=subprocess.PIPE, shell=False)

	return(process2)

#######################################################################
# GENERATE A LIST OF LOCI WHERE VARIANTS OCCUR
#######################################################################

def get_variants(vcf_row, variant_dict, indel_dict, SNP_dict, contig):
	contig_name = vcf_row.split("\t")[0]
	variant_position = int(vcf_row.split("\t")[1])

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
# DevNote - Any way to further speed up if/elif section?
#######################################################################

def find_windows(contig, contig_list, window_len, primer_len, variant_list):
	windows = {}

	for window in range(0,(contig_list[contig] - window_len)):
		window_start = int(window + 1)
		window_stop = int(window_start + window_len)
		primer1_stop = int(window_start + primer_len)
		primer2_start = int(window_stop - primer_len)
		conserved = list(range(window_start,primer1_stop+1)) + list(range(primer2_start,window_stop+1))

	# If no variants in primer sites, save the coordinates

		if not list(set(conserved) & set(variant_list[contig])):
			windows[window_start] = [window_start,primer1_stop,primer2_start,window_stop]

	return(windows)

#######################################################################
# MERGE OVERLAPPING WINDOWS
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

		# Include the final merged window
		merged[saved_window[0]] = [saved_window[0],saved_window[1],saved_window[2],saved_window[3]]
	return(merged)

#######################################################################
# SEQUENCE COMPARISON
#######################################################################

def compare(s, t, query, i):
	if Counter(s) == Counter(t):
		return f"{query} == {i}"
	else:
		pass

#######################################################################
# GET CONSENSUS SEQUENCE FOR BOTH ALLELES AT A LOCUS
# FROM A GIVEN SAMPLE
#	SUBFUNCTION OF VERIFY_WINDOWS()
#######################################################################

def get_allele_consensus(phase, bam, ref, my_range):
	allele = ""
	inbam = bam.replace(".bam", "." + str(phase) + ".bam")
	invcf = bam.replace(".bam", "." + str(phase) + ".vcf.gz")

	cmdD = ["samtools", "faidx", ref, my_range]
	procD = subprocess.Popen(cmdD, stdout=subprocess.PIPE, shell=False)

	cmdE = ["bcftools", "consensus", "-i", "GT=\"hom\"", invcf]
	process3 = subprocess.Popen(cmdE, stdin=procD.stdout, stdout=subprocess.PIPE, shell=False)

	with process3.stdout as result3:
		rows3 = (line.decode() for line in result3)
		for row3 in rows3:
			if not row3.startswith(">"):
				allele += row3.upper().strip("\n")
	return(allele)

#######################################################################
# CHECK COVERAGE OF VARIANTS IN A REGION IS HIGH ENOUGH
#	SUBFUNCTION OF VERIFY_WINDOWS()
# DevNote - currently doesn't work as intended!
#######################################################################

def check_coverage(median_list, bamfile, this_contig, window_start, window_end):
	# Ensure the absence of low-support variants which may be problematic
	## DevNote - adjust to check for suitably-sized overlap between 
	## the window of interest and the low coverage regions

	suitability = True
#	my_median = median_list[bamfile][this_contig]
#	fullvcf = os.path.splitext(bamfile)[0] + ".vcf.gz"
#	reader = vcfpy.Reader.from_path(fullvcf)
#
#	for record in reader.fetch(this_contig, window_start, window_end):
#		if record.INFO['DP'] < (my_median * 0.4):
#			suitability = False
#			break
	return(suitability)

#######################################################################
# ENSURE UNIQUENESS AND COVERAGE OF VARIABLE REGIONS BETWEEN STRAINS
# DevNote - needs speeding up!
#######################################################################

def verify_windows(windows, contig, reference, infiles, medians):
	final = {}

	for window in windows[contig]:
		# Temporary - print checked windows
		with open("CheckedWindows." + contig + ".log", 'a') as outfile:
			outfile.write(str(windows[contig][window]) + "\n")
		outfile.close()
		window_start = windows[contig][window][0]
		window_end = windows[contig][window][3]
		con_range = contig + ":" + str(window_start) + "-" + str(window_end)

		alleles = {}
		results = []
		for bam in infiles:
			# DevNote - this currently generates allele consensuses for all BAM files,
			# even if an earlier one shows poor coverage; use a while condition?
			if check_coverage(medians, bam, contig, window_start, window_end) == True:
				# Add each allele to a list in the relevant nested dictionary
				alleles[bam] = ["", ""]
				for phase in [0, 1]:
					alleles[bam][phase] = get_allele_consensus(phase, bam, reference, con_range)
			else:
				continue

#		# Searches for unique COMBINATIONS of alleles
#		for seq in alleles:
#			query = seq
#			for i in alleles:
#				if query == i:
#					continue
#				results.append(compare(alleles[query], alleles[i], query, i))
#		if all(result is None for result in results):

		if len(alleles) == len(infiles):

			# Searches for at least one unique allele per sample
			all_alleles = [item for sublist in list(alleles.values()) for item in sublist]
			unique_values = [x for x in all_alleles if all_alleles.count(x)==1]

			alleles_are_unique = True
			for sample in alleles.keys():
				if not list(set(alleles[sample]) & set(unique_values)):
					alleles_are_unique = False
					break

			if alleles_are_unique:
				final[window] = windows[contig][window]

				# Re-add the calculator for between-allele differences
				diff_counts = []
				list1 = list(set(all_alleles))
				for i, j in list(combinations(range(0,len(list1)), 2)):
					diff_counts.append(distance(list1[i],list1[j]))
				final[window].append(min(diff_counts))
				final[window].append(max(diff_counts))

			# Get the required sequences for conserved and variable regions
				window_start = windows[contig][window][0]
				primer1_stop = windows[contig][window][1]
				primer2_start = windows[contig][window][2] + 1
				window_stop = windows[contig][window][3]

				full_window = ""

				cmd4 = ["samtools", "faidx", reference, con_range]
				process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, shell=False)
				with process4.stdout as result1:
					rows1 = (line.decode() for line in result1)
					for row1 in rows1:
						if not row1.startswith(">"):
							full_window += row1.strip("\n")
				final[window].append(full_window[(window_start - window_start):(primer1_stop - window_start)])
				final[window].append(full_window[(primer1_stop - window_start):(primer2_start - window_start)])
				final[window].append(full_window[(primer2_start - window_start):((window_stop - window_start)+1)])
			else:
				with open("CheckedWindows." + contig + ".log", 'a') as outfile:
					outfile.write("Not all samples at this locus contain a unique allele.")
				outfile.close()
		else:
			with open("CheckedWindows." + contig + ".log", 'a') as outfile:
				outfile.write("This window was very invalid!")
			outfile.close()
	return(final)

#######################################################################
# PRINT DURATION OF THE STEP
#######################################################################

def print_time(step_name, start_time):
	running_time = step_name + ": " + str(round(time() - start_time, 1)) + " seconds.\n"
	with open("time.log", 'a') as outlog:
		outlog.write(running_time)
	outlog.close()

#######################################################################
# MAIN
#######################################################################

def main(args):
	# Timing - get start time
	full_time = time()
	start_time = time()

	# Set number of threads
	pool = Pool(processes = int(args.threads))

	# Ensure the output files doesn't already exist
	out_bed = args.outprefix + ".bed"
	out_txt = args.outprefix + ".txt"

	if os.path.isfile(out_bed) == True:
		sys.exit("[Error] Output BED file already exists. Please choose another output prefix.")
	if os.path.isfile(out_txt) == True:
		sys.exit("[Error] Output TXT file already exists. Please choose another output prefix.")

	# Record all contig lengths
	contig_lengths = get_contig_lengths(args.sortbam[0])

	# Timing - time taken to get contig lengths
	print_time("Get contig lengths", start_time)
	start_time = time()

	# Record all files' median contig coverages
	cov_stats = {}
	for bam in args.sortbam:
		cov_stats[bam] = {}

	to_coverage = pool.starmap(get_median, \
		[(bam, contig_lengths) for bam in args.sortbam])
	for entry in range(0,len(args.sortbam)):
		cov_stats[args.sortbam[entry]] = to_coverage[entry]

	# Timing - time taken to get contig medians
	print_time("Get median contig coverage stats", start_time)
	start_time = time()

	# Record all low-coverage areas per input BAM
	low_cov = {}
	for bam in args.sortbam:
		low_cov[bam] = {}

	to_lowcov = pool.starmap(get_lowcov, \
		[(bam, contig_lengths, cov_stats) for bam in args.sortbam])
	for entry in range(0,len(args.sortbam)):
		low_cov[args.sortbam[entry]] = to_lowcov[entry]

	# Timing - time taken to get low coverage stats
	print_time("Get low coverage stats", start_time)
	start_time = time()

	with open("low_cov.txt", 'a') as outfile:
		outfile.write(str(low_cov))
	outfile.close()

def main2(args):

	# Timing - get start time
	full_time = time()
	start_time = time()

	# Set number of threads
	pool = Pool(processes = int(args.threads))

	# Ensure the output files doesn't already exist
	out_bed = args.outprefix + ".bed"
	out_txt = args.outprefix + ".txt"

	if os.path.isfile(out_bed) == True:
		sys.exit("[Error] Output BED file already exists. Please choose another output prefix.")
	if os.path.isfile(out_txt) == True:
		sys.exit("[Error] Output TXT file already exists. Please choose another output prefix.")

	# Record all contig lengths
	contig_lengths = get_contig_lengths(args.sortbam[0])

	# Timing - time taken to get contig lengths
	print_time("Get contig lengths", start_time)
	start_time = time()

	# Record all files' median contig coverages
	cov_stats = {}
	for bam in args.sortbam:
		cov_stats[bam] = {}

	to_coverage = pool.starmap(get_median, \
		[(bam, contig_lengths) for bam in args.sortbam])
	for entry in range(0,len(args.sortbam)):
		cov_stats[args.sortbam[entry]] = to_coverage[entry]

	# Timing - time taken to get contig medians
	print_time("Get median contig coverage stats", start_time)
	start_time = time()

	# Record all low-coverage areas per input BAM
	low_cov = {}
	for bam in args.sortbam:
		low_cov[bam] = {}

	to_lowcov = pool.starmap(get_lowcov, \
		[(bam, contig_lengths, cov_stats) for bam in args.sortbam])
	for entry in range(0,len(args.sortbam)):
		low_cov[args.sortbam[entry]] = to_lowcov[entry]

	# Timing - time taken to get low coverage stats
	print_time("Get low coverage stats", start_time)
	start_time = time()

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


	print("Searching for potential barcodes in",len(args.sortbam),"file(s).")

	for bam in args.sortbam:
		vcf_file = NoExt(bam) + ".vcf"

		# Generate or read in a VCF file for the current BAM
		process2 = bcf(bam, contig_lengths, filter_qual, args.threads, args.ref)

		# Generate a list of loci where variants occur
		print("\nFinding variants...")
		current_contig = ""
		if os.path.isfile(NoExt(bam) + ".vcf.gz") == True:
			if os.path.isfile(NoExt(bam) + ".vcf.gz.csi") == False:
				os.system("bcftools index --threads " + str(args.threads) + " " + vcf_file + ".gz")
			for row2 in process2:
				if not row2.startswith("#"):
					current_contig = get_variants(row2, all_variants, all_indels, all_SNPs, current_contig)

		# If a VCF file doesn't already exist, generate one, generate variant the list, then bgzip and index the new VCF
		else:
			with process2.stdout as result2, open(vcf_file, "a") as output_file:
				rows2 = (line.decode() for line in result2)
				for row2 in rows2:
					output_file.write(row2)
					if not row2.startswith("#"):
						current_contig = get_variants(row2, all_variants, all_indels, all_SNPs, current_contig)
			output_file.close()
			os.system("bgzip -@ " + str(args.threads) + " " + vcf_file)
			os.system("bcftools index -f --threads " + str(args.threads) + " " + vcf_file + ".gz")

	# Get sorted lists of variant/SNP/indel positions, with duplicates removed
	for contig in all_variants:
		all_variants[contig] = sorted(list(set(all_variants[contig])), key=int)
	for contig in all_SNPs:
		all_SNPs[contig] = sorted(list(set(all_SNPs[contig])), key=int)
	for contig in all_indels:
		all_indels[contig] = sorted(list(set(all_indels[contig])), key=int)

	# Timing - time taken to get lists of variants
	print_time("Get lists of variants", start_time)
	start_time = time()

	# Step through each contig, assigning start and stop locations for window and primers
	# DevNote - this needs speeding up!
	print("\nChecking windows...")

	to_master = pool.starmap(find_windows, \
	[(contig, contig_lengths, args.window_size, args.primer_size, all_variants) for contig in contig_lengths])

	for entry in range(0,len(contig_lengths)):
		master_dict[list(contig_lengths.keys())[entry]] = to_master[entry]

	# Timing - time taken to find valid windows
	print_time("Find valid windows", start_time)
	start_time = time()

	# Merge overlapping windows
	print("\nMerging overlapping windows...")

	for contig in master_dict:
		if master_dict[contig]:
			print(contig)
			merged_dict[contig] = merge_windows(contig, master_dict)

	# Timing - time taken to merge windows
	print_time("Merge windows", start_time)
	start_time = time()

# Dev Note: Find a way to skip loci where the following type of error occurs:
# `Warning: ignoring overlapping variant starting at 000215F:2953`

	# Compare consensus sequences for all BAMs, to ensure they are truly unique,
	# not merely differing from the reference in all the same positions
	# DevNote - this needs speeding up!

	# DevNote - Update in progress to interrogate phased information

	print("\nChecking consensuses...")

	to_final = pool.starmap(verify_windows, \
		[(merged_dict, contig, args.ref, args.sortbam, cov_stats) for contig in contig_lengths])

	for entry in range(0,len(contig_lengths)):
		final_dict[list(contig_lengths.keys())[entry]] = to_final[entry]

	# Timing - time taken to get unique windows
	print_time("Get unique windows", start_time)
	start_time = time()

	# Report results in TXT and BED format
	print("\nGenerating output files...")

	with open(out_bed, "a") as output_bed, open(out_txt, "a") as output_txt:
		output_bed.write("track name=PotentialBarcodes description=\"Potential barcodes\"\n")

#		output_txt.write("window_name\tcontig\tconserved_1_start\tconserved_1_end\tconserved_1_seq\tvariable_start\tvariable_end\tvariable_seq\tconserved_2_start\tconserved_2_end\tconserved_2_seq\tvariable_length\n")
		output_txt.write("window_name\tcontig\tconserved_1_start\tconserved_1_end\tconserved_1_seq\tvariable_start\tvariable_end\tvariable_seq\tconserved_2_start\tconserved_2_end\tconserved_2_seq\tvariable_length\tmin_diffs\tmax_diffs\n")

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
				conserved_1_seq = final_dict[contig][window][6]
				variable_seq = final_dict[contig][window][7]
				conserved_2_seq = final_dict[contig][window][8]

#				conserved_1_seq = final_dict[contig][window][4]
#				variable_seq = final_dict[contig][window][5]
#				conserved_2_seq = final_dict[contig][window][6]

				for SNP in all_SNPs[contig]:
					if variable_start < int(SNP) < variable_stop:
						window_SNPs += 1
				for indel in all_indels[contig]:
					if variable_start < int(indel) < variable_stop:
						window_indels += 1

				window_name = str(contig) + "_" + str(window_number) + "_SNPs_" + str(window_SNPs) + "_indels_" + str(window_indels)

				# Fields 7 and 8 (thickStart and thickEnd) represent the start and stop positions of the non-primer part of the window
				window_out = str(contig) + "\t" + str(conserved_1_start - 1) + "\t" + str(conserved_2_stop) + "\t" + \
						str(window_name) + "\t0\t.\t" + str(conserved_1_stop) + "\t" + str(variable_stop) + "\n"
				line_out = str(window_name) + "\t" + str(contig) + "\t" + str(conserved_1_start) + "\t" + str(conserved_1_stop) + "\t" + conserved_1_seq + "\t" + \
						str(variable_start) + "\t" + str(variable_stop) + "\t" + variable_seq + "\t" + \
						str(conserved_2_start) + "\t" + str(conserved_2_stop) + "\t" + conserved_2_seq + "\t" + \
						str(variable_len) + "\t" + str(min_diffs) + "\t" + str(max_diffs) + "\n"
#				line_out = str(window_name) + "\t" + str(contig) + "\t" + str(conserved_1_start) + "\t" + str(conserved_1_stop) + "\t" + conserved_1_seq + "\t" + \
#						str(variable_start) + "\t" + str(variable_stop) + "\t" + variable_seq + "\t" + \
#						str(conserved_2_start) + "\t" + str(conserved_2_stop) + "\t" + conserved_2_seq + "\t" + \
#						str(variable_len) + "\n"

				output_bed.write(window_out)
				output_txt.write(line_out)
		output_bed.close()
		output_txt.close()

	# Timing - time taken to write output
	print_time("Write output", start_time)

	# Timing - total pipeline time
	total_time = "Total time: " + str(round(time() - full_time, 1)) + " seconds.\n"
	with open("time.log", 'a') as outlog:
		outlog.write("\n")
		outlog.write(total_time)
	outlog.close()

#######################################################################

if __name__ == "__main__":
	main()
