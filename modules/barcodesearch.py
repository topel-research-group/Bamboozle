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
import vcfpy
import gzip
import datetime
from itertools import combinations
from Levenshtein import distance
from multiprocessing import Pool
from time import time
from collections import Counter
from statistics import median

#######################################################################
# DEFINE WINDOW CLASS
#######################################################################

class Window:
	def __init__(self, contig, winstart, p1stop, p2start, winstop):
		self.contig = contig
		self.winstart = winstart
		self.p1stop = p1stop
		self.p2start = p2start
		self.winstop = winstop
	def __str__(self):
		return f"[{self.contig}, {self.winstart}, {self.p1stop}, {self.p2start}, {self.winstop}]"

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
# IDENTIFY AREAS OF IRREGULAR COVERAGE FOR EACH SAMPLE
#	DevNote - Currently, bad regions are defined as regions
#		either less than half or more than twice the contig
#		median; this should be adjusted
#	DevNote - Any efficient way to merge adjacent windows?
#######################################################################

def get_badcov(bam, contigs, coverage_stats):
	bad_cov = {}
	for contig in contigs:
		bad_cov[contig] = {}

	cmdX = ["bedtools", "genomecov", "-bga", "-ibam", bam]
	procX = subprocess.Popen(cmdX, stdout=subprocess.PIPE, shell=False)

	with procX.stdout as result:
		rows = (line.decode().split("\t") for line in result)
		for row in rows:
			ctg = str(row[0])
			start_pos = int(row[1]) + 1
			end_pos = int(row[2])
			coverage = int(row[3])
			if (coverage > coverage_stats[ctg]*2) or (coverage < coverage_stats[ctg]*0.5):
				bad_cov[ctg][start_pos] = [start_pos, end_pos]
	return(bad_cov)

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
# OUT: [Window1, Window2, Window3, ...]
# DevNote - Any way to further speed up if section?
#######################################################################

def find_windows(contig, con_len, window_len, primer_len, variant_list, logfile):
	windows = []

	for window in range(0,(con_len - window_len)):
		window_start = int(window + 1)
		window_stop = int(window_start + window_len)
		primer1_stop = int(window_start + primer_len)
		primer2_start = int(window_stop - primer_len)
		conserved = list(range(window_start,primer1_stop+1)) + list(range(primer2_start,window_stop+1))

	# If no variants in primer sites, save the coordinates

		if not list(set(conserved) & set(variant_list)):
			windows.append(Window(contig,window_start,primer1_stop,primer2_start,window_stop))

	with open(logfile, 'a') as outfile:
		outfile.write(str(len(windows)) + " preliminary windows found in contig " + contig + ".\n")
	outfile.close()

	return(windows)

#######################################################################
# MERGE OVERLAPPING WINDOWS
# OUT: [Window1, Window2, Window3, ...]
#######################################################################

def merge_windows(window_list, contig):
	merged = []
	saved_window = Window(contig,0,0,0,0)

	for window in window_list:
		if saved_window == Window(contig,0,0,0,0):
			saved_window = Window(contig,window.winstart, window.p1stop, window.p2start, window.winstop)
		elif saved_window.winstart <= window.winstart <= saved_window.p1stop:
			saved_window = Window(contig,saved_window.winstart, window.p1stop, saved_window.p2start, window.winstop)
		else:
			merged.append(saved_window)
			saved_window = Window(contig,window.winstart, window.p1stop, window.p2start, window.winstop)

	# Include the final merged window
	merged.append(saved_window)

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
# CHECK COVERAGE OF VARIANTS IN A REGION IS GOOD ENOUGH
#	SUBFUNCTION OF VERIFY_WINDOWS()
#######################################################################

def check_coverage(median_list, bad_regions, bamfile, this_contig, window_start, window_end):
	# Ensure the window doesn't heavily overlap a bad-coverage region
	## In this case, render unsuitable if such a region overlaps the
	## proposed window by 10 bases (saves < 3aa indels)

	suitability = True

	window_range = set(range(window_start, window_end + 1))

	bad_overlap = 10

	for badplace in bad_regions:
		bad_start = bad_regions[badplace][0]
		bad_end = bad_regions[badplace][1] + 1
		bad_range = range(bad_start, bad_end)

		if bad_start > window_end:
			break
		elif bad_end < window_start:
			continue
		elif len(window_range.intersection(bad_range)) >= bad_overlap:
			suitability = False
			break

	return(suitability)

#######################################################################
# ENSURE GOOD COVERAGE OF VARIABLE REGIONS BETWEEN STRAINS
# OUT: [Window1, Window2, Window3, ...]
# DevNote - formerly part of verify_windows()
# DevNote - this currently generates allele consensuses for all BAM files,
# even if an earlier one shows poor coverage; use a while condition?
#######################################################################

def good_coverage(windows, median_list, bad_regions, infiles, contig, logfile):
	good_windows = []

	for window in windows:
		suitability = True
		for bam in infiles:
			window_range = set(range(window.winstart, window.winstop + 1))
			bad_overlap = 10

			for badplace in bad_regions[bam][window.contig]:
				bad_start = bad_regions[bam][window.contig][badplace][0]
				bad_end = bad_regions[bam][window.contig][badplace][1] + 1
				bad_range = range(bad_start, bad_end)

				if bad_start > window.winstop:
					break
				elif bad_end < window.winstart:
					continue
				elif len(window_range.intersection(bad_range)) >= bad_overlap:
					suitability = False
					break
		if suitability == True:
			good_windows.append(window)

	with open(logfile, 'a') as outfile:
		outfile.write(str(len(good_windows)) + " good coverage windows found in contig " + contig + ".\n")
	outfile.close()

	return(good_windows)

#######################################################################
# CHUNK UP good_cov_list TO SHARE EQUALLY BETWEEN PROCESSES IN THE
#  FINAL STEP OF THE PIPELINE
#######################################################################

def chunks(lst, n):
	for i in range(0, len(lst), n):
		yield lst[i:i + n]

#######################################################################
# ENSURE UNIQUENESS AND COVERAGE OF VARIABLE REGIONS BETWEEN STRAINS
# OUT: [[Window, min_diffs, max_diffs, p1_seq, var_seq, p2_seq], [...]]
# DevNote - needs speeding up!
#######################################################################

def verify_windows(windows, reference, infiles, medians, badcov):
	final = []

	good_windows = 0

	for window in windows:
		con_range = window.contig + ":" + str(window.winstart) + "-" + str(window.winstop)
		alleles = {}
		results = []
		for bam in infiles:
			# Add each allele to a list in the relevant nested dictionary
			alleles[bam] = ["", ""]
			for phase in [0, 1]:
				alleles[bam][phase] = get_allele_consensus(phase, bam, reference, con_range)

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
				final.append([])
				final[good_windows].append(window)

			# Calculate between-allele differences
				diff_counts = []
				list1 = list(set(all_alleles))
				for i, j in list(combinations(range(0,len(list1)), 2)):
					diff_counts.append(distance(list1[i],list1[j]))
				final[good_windows].append(min(diff_counts))
				final[good_windows].append(max(diff_counts))

			# Get the required sequences for conserved and variable regions

				full_window = ""

				cmd4 = ["samtools", "faidx", reference, con_range]
				process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, shell=False)
				with process4.stdout as result1:
					rows1 = (line.decode() for line in result1)
					for row1 in rows1:
						if not row1.startswith(">"):
							full_window += row1.strip("\n")
				final[good_windows].append(full_window[0:(window.p1stop - window.winstart)])
				final[good_windows].append(full_window[(window.p1stop - window.winstart):((window.p2start + 1) - window.winstart)])
				final[good_windows].append(full_window[((window.p2start + 1) - window.winstart):((window.winstop - window.winstart) + 1)])
				good_windows += 1

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

# For debugging input_files.py
def main2(args):
	for infile in args.sortbam:
		for FileExt in [".0.bam", ".0.bai", ".0.vcf.gz", ".0.vcf.gz.csi", ".1.bam", ".1.bai", ".1.vcf.gz", ".1.vcf.gz.csi"]:
			if os.path.isfile(NoExt(infile) + FileExt):
				print(NoExt(infile) + FileExt + " is present!")
			else:
				print("Error: " + NoExt(infile) + FileExt + " is missing!")

def main(args):

	#######################################################################
	# STEP 1 - LAYING THE GROUNDWORK
	#######################################################################

	# Timing - get start time
	full_time = time()

	# Set number of threads
	pool = Pool(processes = int(args.threads))

	# Ensure the output files doesn't already exist
	out_bed = args.outprefix + ".bed"
	out_txt = args.outprefix + ".txt"

	barcode_log = barcode_log = "barcoding." + datetime.datetime.now().strftime("%d-%b-%Y") + ".log"

	if os.path.isfile(out_bed):
		sys.exit("[Error] Output BED file already exists. Please choose another output prefix.")
	if os.path.isfile(out_txt):
		sys.exit("[Error] Output TXT file already exists. Please choose another output prefix.")
	if os.path.isfile(barcode_log):
		sys.exit("[Error] Another log file exists from today. Please rename or delete it and retry.")

	#######################################################################
	# STEP 2 - GET CONTIG LENGTH STATISTICS
	#	OUT: contig_lengths = {contig1: length}
	#######################################################################
	# Included custom functions: get_contig_lengths
	#######################################################################
	start_time = time()

	contig_lengths = get_contig_lengths(args.sortbam[0])

	# Timing - time taken to get contig lengths
	print_time("Get contig lengths", start_time)

	########################################################################
	# STEP 3 - SET SOME ADDITIONAL VARIABLES
	########################################################################

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
		master_dict[contig] = []
	merged_dict = {}
	for contig in contig_lengths:
		merged_dict[contig] = []

	good_cov_dict = {}
	for contig in contig_lengths:
		good_cov_dict[contig] = []

#	final_dict = {}
#	for contig in contig_lengths:
#		final_dict[contig] = []

	# Set other values
	filter_qual = "%QUAL>" + str(args.quality)

	#######################################################################
	# STEP 4 - GET CONTIG MEDIANS PER BAM FILE
	#	OUT: cov_stats = {bam1: {contig1: median}}
	#######################################################################
	# Included custom functions: get_median
	#######################################################################
	start_time = time()

	cov_stats = {}
	for bam in args.sortbam:
		cov_stats[bam] = {}

	to_coverage = pool.starmap(get_median, \
		[(bam, contig_lengths) for bam in args.sortbam])
	for entry in range(0,len(args.sortbam)):
		cov_stats[args.sortbam[entry]] = to_coverage[entry]

	# Timing - time taken to get contig medians
	print_time("Get median contig coverage stats", start_time)

	#######################################################################
	# STEP 5 - GET IRREGULAR COVERAGE REGIONS PER BAM FILE
	#	OUT: bad_cov = {bam1: {contig1: {window1: [start, end]}}}
	#######################################################################
	# Included custom functions: get_badcov
	#######################################################################
	start_time = time()

	bad_cov = {}
	for bam in args.sortbam:
		bad_cov[bam] = {}

	to_badcov = pool.starmap(get_badcov, \
		[(bam, contig_lengths, cov_stats[bam]) for bam in args.sortbam])
	for entry in range(0,len(args.sortbam)):
		bad_cov[args.sortbam[entry]] = to_badcov[entry]

	# Timing - time taken to get bad coverage stats
	print_time("Get irregular coverage stats", start_time)

	#######################################################################
	# STEP 6 - GENERATE A LIST OF POSITIONS WHERE VARIANTS OCCUR
	#	OUT: all_variants = {contig1: [var1, var2, var3]}
	#		all_SNPs = {contig1: [SNP1, SNP2, SNP3]}
	#		all_indels = {contig1: [ind1, ind2, ind3]}
	#######################################################################
	# Included custom functions: get_variants
	#######################################################################
	start_time = time()

	print("Searching for potential barcodes in",len(args.sortbam),"file(s).")

	# HOW TO FIX THIS SO IT READS GZIP FILES STRAIGHTAWAY?

	for bam in args.sortbam:
		print("\nFinding variants...")
		current_contig = ""

		for row2 in gzip.open(NoExt(bam) + ".vcf.gz", 'rt'):
			if not row2.startswith("#"):
				current_contig = get_variants(row2, all_variants, all_indels, all_SNPs, current_contig)

	# Get sorted lists of variant/SNP/indel positions, with duplicates removed
	for contig in all_variants:
		all_variants[contig] = sorted(list(set(all_variants[contig])), key=int)
	for contig in all_SNPs:
		all_SNPs[contig] = sorted(list(set(all_SNPs[contig])), key=int)
	for contig in all_indels:
		all_indels[contig] = sorted(list(set(all_indels[contig])), key=int)

	# Timing - time taken to get lists of variants
	print_time("Get lists of variants", start_time)

	#######################################################################
	# STEP 7 - STEP THROUGH EACH CONTIG LOOKING FOR WINDOWS WITH
	#		VARIABLE CENTRES AND CONSERVED PRIMER SITES
	# OUT: master_dict = {contig1: [Window1, Window2, Window3]}
	#######################################################################
	# Included custom functions: find_windows
	#######################################################################
	# DevNote - this needs speeding up!
	#######################################################################
	start_time = time()

	print("\nChecking windows...")

	to_master = pool.starmap(find_windows, \
	[(contig, contig_lengths[contig], args.window_size, args.primer_size, all_variants[contig], barcode_log) for contig in contig_lengths])

	for entry in range(0,len(contig_lengths)):
		master_dict[list(contig_lengths.keys())[entry]] = to_master[entry]

	# Timing - time taken to find valid windows
	print_time("Find valid windows", start_time)

	#######################################################################
	# STEP 8 - MERGE OVERLAPPING WINDOWS
	# OUT: merged_dict = {contig1: [Window1, Window2, Window3]}
	#######################################################################
	# Included custom functions: merge_windows
	#######################################################################
	start_time = time()

	print("\nMerging overlapping windows...")

	for contig in master_dict:
		if master_dict[contig]:
			print(contig)
			merged_dict[contig] = merge_windows(master_dict[contig], contig)

#	merged_list = [item for sublist in list(merged_dict.values()) for item in sublist]

	# Timing - time taken to merge windows
	print_time("Merge windows", start_time)

	#######################################################################
	# STEP 8.5 - EXCLUDE WINDOWS OF LOW COVERAGE
	# OUT: good_cov_dict = {contig1: [Window1, Window2, Window3]}
	#######################################################################
	start_time = time()

	print("\nChecking window coverage...")

	to_good_cov = pool.starmap(good_coverage, \
		[(merged_dict[contig], cov_stats, bad_cov, args.sortbam, contig, barcode_log) for contig in contig_lengths])

	for entry in range(0,len(contig_lengths)):
		good_cov_dict[list(contig_lengths.keys())[entry]] = to_good_cov[entry]

	good_cov_list = [item for sublist in list(good_cov_dict.values()) for item in sublist]

	# Timing - time taken to get good-coverage windows
	print_time("Get good-coverage windows", start_time)

	#######################################################################
	# STEP 9 - ENSURE THAT EACH SAMPLE CONTAINS AT LEAST ONE UNIQUE ALLELE
	#		AT EACH LOCUS, AND THAT THE COVERAGE IS ACCEPTABLE
	# OUT: final_dict = {contig1: [[Window1, min_diffs, max_diffs, p1_seq, var_seq, p2_seq],[...]]}
	#######################################################################
	# Included custom functions: verify_windows
	#	Subfunctions: compare, get_allele_consensus, check_coverage
	#######################################################################
	# DevNote - find a way to skip loci where the following type of error occurs:
	# `Warning: ignoring overlapping variant starting at 000215F:2953`
	#######################################################################
	# DevNote - this needs speeding up!
	#######################################################################
	start_time = time()

	print("\nChecking consensuses...")

	final_list = []

	chunk_len = int(len(good_cov_list) / int(args.threads))

	to_final = pool.starmap(verify_windows, \
		[(chunky, args.ref, args.sortbam, cov_stats, bad_cov) for chunky in chunks(good_cov_list, chunk_len)])

	chunk_no = int(args.threads)

	for entry in range(0,chunk_no):
		final_list.append(to_final[entry])

	really_final_list = [item for sublist in final_list for item in sublist]

	# Timing - time taken to get unique windows
	print_time("Get unique windows", start_time)

	#######################################################################
	# STEP 10 - REPORT RESULTS IN TXT AND BED FORMAT
	#######################################################################
	start_time = time()

	print("\nGenerating output files...")

	with open(out_bed, "a") as output_bed, open(out_txt, "a") as output_txt:
		output_bed.write("track name=PotentialBarcodes description=\"Potential barcodes\"\n")

		output_txt.write("window_name\tcontig\tconserved_1_start\tconserved_1_end\tconserved_1_seq\tvariable_start\tvariable_end\tvariable_seq\tconserved_2_start\tconserved_2_end\tconserved_2_seq\tvariable_length\tmin_diffs\tmax_diffs\n")

		window_number = 0
		for window in really_final_list:
			window_number += 1
			window_SNPs = 0
			window_indels = 0

			conserved_1_start = window[0].winstart
			conserved_1_stop = window[0].p1stop - 1
			variable_start = window[0].p1stop
			variable_stop = window[0].p2start
			conserved_2_start = window[0].p2start + 1
			conserved_2_stop = window[0].winstop
			variable_len = conserved_2_start - variable_start

			min_diffs = window[1]
			max_diffs = window[2]
			conserved_1_seq = window[3]
			variable_seq = window[4]
			conserved_2_seq = window[5]

			for SNP in all_SNPs[window[0].contig]:
				if variable_start < int(SNP) < variable_stop:
					window_SNPs += 1
			for indel in all_indels[window[0].contig]:
				if variable_start < int(indel) < variable_stop:
					window_indels += 1

			window_name = str(window[0].contig) + "_" + str(window_number) + "_SNPs_" + str(window_SNPs) + "_indels_" + str(window_indels)

			# Fields 7 and 8 (thickStart and thickEnd) represent the start and stop positions of the non-primer part of the window
			window_out = str(window[0].contig) + "\t" + str(conserved_1_start - 1) + "\t" + str(conserved_2_stop) + "\t" + \
					str(window_name) + "\t0\t.\t" + str(conserved_1_stop) + "\t" + str(variable_stop) + "\n"
			line_out = str(window_name) + "\t" + str(window[0].contig) + "\t" + str(conserved_1_start) + "\t" + str(conserved_1_stop) + "\t" + conserved_1_seq + "\t" + \
					str(variable_start) + "\t" + str(variable_stop) + "\t" + variable_seq + "\t" + \
					str(conserved_2_start) + "\t" + str(conserved_2_stop) + "\t" + conserved_2_seq + "\t" + \
					str(variable_len) + "\t" + str(min_diffs) + "\t" + str(max_diffs) + "\n"

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
