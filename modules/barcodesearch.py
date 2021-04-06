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
import pickle
import shutil
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

	cmd = ["samtools", "depth", bamfile]
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=False)
	current_contig = "None"
	this_contig = []
	with process.stdout as result:
		rows = (line.decode().split('\t') for line in result)
		for row in rows:
			ctg = str(row[0])
			position = int(row[1])
			coverage = int(row[2])
			if current_contig == "None":
				current_contig = ctg
			if ctg != current_contig:
				coverage_stats[current_contig] = median(this_contig)
				this_contig = []
				current_contig = ctg
			this_contig.append(coverage)
	coverage_stats[current_contig] = median(this_contig)
	
	return(coverage_stats)

#######################################################################
# IDENTIFY AREAS OF IRREGULAR COVERAGE FOR EACH SAMPLE
#	DevNote - Currently, bad regions are defined as regions
#		either less than half or more than twice the contig
#		median; this should be adjusted
#######################################################################

def get_badcov(bam, contigs, coverage_stats):
	bad_cov = {}
	for contig in contigs:
		bad_cov[contig] = []

	cmdX = ["bedtools", "genomecov", "-bga", "-ibam", bam]
	procX = subprocess.Popen(cmdX, stdout=subprocess.PIPE, shell=False)

	with procX.stdout as result:
		rows = (line.decode().split("\t") for line in result)
		this_window = Window("","","","","")
		for row in rows:
			ctg = str(row[0])
			start_pos = int(row[1]) + 1
			end_pos = int(row[2])
			coverage = int(row[3])

		# If this is the first window checked, save the contig name
			if this_window.contig == "":
				this_window.contig = ctg

		# Is it the same contig as the previous window?
		# If not, save the existing window, and shift contig
			if this_window.contig != ctg and this_window != Window(ctg,"","","",""):
				bad_cov[this_window.contig].append([this_window.winstart, this_window.winstop])
				this_window = Window(ctg,"","","","")

		# Check whether the window fulfils the criteria
			if (coverage > coverage_stats[ctg]*2) or (coverage < coverage_stats[ctg]*0.5):

		# If yes, and this is the first window in the contig, save the position
				if this_window.winstart == "":
					this_window = Window(ctg,start_pos,"","",end_pos)

		# If yes, and it's adjacent to the previous saved window, extend the existing window
				elif start_pos == (this_window.winstop + 1):
					this_window.winstop = end_pos

		# Otherwise, export the existing window to bad_cov and save the new one
				else:
					bad_cov[this_window.contig].append([this_window.winstart, this_window.winstop])
					this_window = Window(ctg,start_pos,"","",end_pos)

		# Output the final window
		if this_window != Window(ctg,"","","",""):
			bad_cov[this_window.contig].append([this_window.winstart, this_window.winstop])

	return(bad_cov)

#######################################################################
# MERGE OVERLAPPING BADCOV WINDOWS FROM DIFFERENT BAM FILES
# Code adapted from https://stackoverflow.com/questions/58535825/combine-overlapping-ranges-of-numbers
#######################################################################

def any_items_overlap(l):

	# For each possible pair of sublists
	for item1, item2 in combinations(l, 2):

		min1, max1 = item1
		min2, max2 = item2

		# If no overlap, ignore this pair
		if min1 > max2 or max1 < min2:
			continue

		# If overlap exists, return pair
		else:
			return item1, item2

	return None

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

	refgen = vcf_row.split("\t")[3]
	this_var = vcf_row.split("\t")[4]

	if "," in this_var:
		var1 = this_var.split(",")[0]
		var2 = this_var.split(",")[1]
		if len(var1) != len(refgen) or len(var2) != len(refgen):
			indel_dict[contig_name].append(variant_position)
		else:
			SNP_dict[contig_name].append(variant_position)
	elif len(this_var) != len(refgen):
		indel_dict[contig_name].append(variant_position)
	else:
		SNP_dict[contig_name].append(variant_position)

	return(contig)

#######################################################################
# IDENTIFY VARIABLE WINDOWS AND CONSERVED PRIMER SITES
# OUT: [Window1, Window2, Window3, ...]
# DevNote - Any way to further speed up if section?
#######################################################################

def find_windows(contig, con_len, window_len, primer_len, variant_list, logfile, tempdirectory):
	tempfile = tempdirectory + "/find_windows." + contig + ".pickle"
	if os.path.isfile(tempfile):
		with open(tempfile, "rb") as infile:
			windows = pickle.load(infile)
		infile.close()

	else:
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

		with open(logfile, 'a') as outfile1, open(tempfile, 'wb') as outfile2:
			outfile1.write(str(len(windows)) + " preliminary windows found in contig " + contig + ".\n")
			pickle.dump(windows, outfile2)
		outfile1.close()
		outfile2.close()

	return(windows)

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
# AS get_allele_consensus(), BUT WHEN DEALING WITH HAPLOID SAMPLES
#       SUBFUNCTION OF VERIFY_WINDOWS()
# DevNote - merge with get_allele_consensus()?
#######################################################################

def get_consensus_haploid(bam, ref, my_range):
	allele = ""
	invcf = bam.replace(".bam", ".vcf.gz")

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
# ENSURE GOOD COVERAGE OF VARIABLE REGIONS BETWEEN STRAINS
# OUT: [Window1, Window2, Window3, ...]
# DevNote - formerly part of verify_windows()
# DevNote - this currently generates allele consensuses for all BAM files,
# even if an earlier one shows poor coverage; use a while condition?
#######################################################################

def good_coverage(windows, median_list, bad_regions, badcov_threshold, infiles, contig, logfile):
	good_windows = []

	for window in windows:
		suitability = True
		for bam in infiles:
			window_range = set(range(window.winstart, window.winstop + 1))
			bad_overlap = badcov_threshold

			for badplace in bad_regions[window.contig]:
				bad_start = badplace[0]
				bad_end = badplace[1] + 1
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
# CHUNK UP good_cov_list TO SHARE EQUALLY BETWEEN PROCESSES IN THE
#  FINAL STEP OF THE PIPELINE
#######################################################################

def chunks(lst, n):
	for i in range(0, len(lst), n):
		yield lst[i:i + n]

#######################################################################
# ENSURE UNIQUENESS AND COVERAGE OF VARIABLE REGIONS BETWEEN STRAINS
# OUT: [[Window, min_diffs, max_diffs, p1_seq, var_seq, p2_seq], [...]]
# DevNote - Needs adjusting to accept haploids
#######################################################################

def verify_windows(windows, reference, infiles, medians, ploidy, out_dir):
	final = []

	good_windows = 0

	for window in windows:
		con_range = window.contig + ":" + str(window.winstart) + "-" + str(window.winstop)
		alleles = {}
		results = []

		if ploidy == "haploid":
			for bam in infiles:
				alleles[bam] = get_consensus_haploid(bam, reference, con_range)

		elif ploidy == "diploid":
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

			if ploidy == "diploid":
				# Searches for at least one unique allele per sample
				all_alleles = [item for sublist in list(alleles.values()) for item in sublist]
			if ploidy == "haploid":
				all_alleles = list(alleles.values())
			unique_alleles = [x for x in all_alleles if all_alleles.count(x)==1]

			# Compare set of sample's alleles and all unique alleles
			# Each sample should contain at least one unique allele
			# For haploids, all alleles should be unique
			alleles_are_unique = True

			if ploidy == "haploid":
				if len(all_alleles) != len(unique_alleles):
					alleles_are_unique = False

			if ploidy == "diploid":
				for sample in alleles.keys():
					if not list(set(alleles[sample]) & set(unique_alleles)):
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

			# Print the sequences to a FASTA file in the FASTA output directory
				output_file = out_dir + "/" + window.contig + "_" + str(window.winstart) + "-" + str(window.winstop) + "_alleles.fasta"
				with open(output_file, "a") as outfile:
					for bam in infiles:
						if ploidy == "haploid":
							fasta_header = FileName(bam) + "_" + window.contig + "_" + str(window.winstart) + "-" + str(window.winstop)
							outfile.write(">" + fasta_header + "\n")
							outfile.write(alleles[bam] + "\n")
						elif ploidy == "diploid":
							for phase in [0, 1]:
								fasta_header = FileName(bam) + "_" + window.contig + "_" + str(window.winstart) + "-" + str(window.winstop) + "_alleles_" + str(phase)
								outfile.write(">" + fasta_header + "\n")
								outfile.write(alleles[bam][phase] + "\n")
				outfile.close()

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
	# LAYING THE GROUNDWORK
	#######################################################################

	# Timing - get start time
	full_time = time()

	# Set number of threads
	pool = Pool(processes = int(args.threads))

	# If using --resume, check which step to continue from
	if args.resume:
		step_list = ["Get contig lengths", "Get median contig coverage stats", "Get irregular coverage stats", "Get lists of variants", \
			"Find valid windows", "Merge windows", "Get good-coverage windows", "Get unique windows"]

		with open("time.log", "r") as logfile:
			lines = logfile.read().splitlines()
		logfile.close()
		gofrom = step_list.index(lines[-1].split(":")[0]) + 1
	else:
		gofrom = 0

	# Ensure the output files doesn't already exist
	out_bed = args.outprefix + ".bed"
	out_txt = args.outprefix + ".txt"
	out_fasta_dir = args.outprefix + "_alleles"
	out_pickle_tmp = args.outprefix + "_tmpdir"

	barcode_log = barcode_log = "barcoding." + datetime.datetime.now().strftime("%d-%b-%Y") + ".log"

	if os.path.isfile(out_bed):
		sys.exit("[Error] Output BED file already exists. Please choose another output prefix.")
	if os.path.isfile(out_txt):
		sys.exit("[Error] Output TXT file already exists. Please choose another output prefix.")
	if os.path.isfile(barcode_log) and not args.resume:
		sys.exit("[Error] Another log file exists from today. Please rename or delete it and retry.")

	if os.path.isdir(out_fasta_dir) and not args.resume:
		sys.exit("[Error] Output FASTA directory already exists. Please choose another output prefix.")

	if not args.resume:
		os.mkdir(out_pickle_tmp)

	#######################################################################
	# STEP 1 - GET CONTIG LENGTH STATISTICS
	#	OUT: contig_lengths = {contig1: length}
	#######################################################################
	# Included custom functions: get_contig_lengths
	#######################################################################
	if gofrom < 1:
		start_time = time()

		contig_lengths = get_contig_lengths(args.sortbam[0])

		# Export contig_lengths to pickle
		with open(out_pickle_tmp + "/contig_lengths.pickle", "wb") as outfile:
			pickle.dump(contig_lengths, outfile)
		outfile.close()

		# Timing - time taken to get contig lengths
		print_time("Get contig lengths", start_time)

	else:
		with open(out_pickle_tmp + "/contig_lengths.pickle", "rb") as infile:
			contig_lengths = pickle.load(infile)
		infile.close()

	########################################################################
	# SET SOME ADDITIONAL VARIABLES
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
	good_cov_dict = {}
	for contig in contig_lengths:
		good_cov_dict[contig] = []
	merged_dict = {}
	for contig in contig_lengths:
		merged_dict[contig] = []

	# Set other values
	filter_qual = "%QUAL>" + str(args.quality)

	#######################################################################
	# STEP 2 - GET CONTIG MEDIANS PER BAM FILE
	#	OUT: cov_stats = {bam1: {contig1: median}}
	#######################################################################
	# Included custom functions: get_median
	#######################################################################
	if gofrom < 2:
		start_time = time()

		cov_stats = {}
		for bam in args.sortbam:
			cov_stats[bam] = {}

		to_coverage = pool.starmap(get_median, \
			[(bam, contig_lengths) for bam in args.sortbam])
		for entry in range(0,len(args.sortbam)):
			cov_stats[args.sortbam[entry]] = to_coverage[entry]

		# Export cov_stats to pickle
		with open(out_pickle_tmp + "/cov_stats.pickle", "wb") as outfile:
			pickle.dump(cov_stats, outfile)
		outfile.close()

		# Timing - time taken to get contig medians
		print_time("Get median contig coverage stats", start_time)

	else:
		with open(out_pickle_tmp + "/cov_stats.pickle", "rb") as infile:
			cov_stats = pickle.load(infile)
		infile.close()

	#######################################################################
	# STEP 3 - GET IRREGULAR COVERAGE REGIONS PER BAM FILE
	#	OUT: bad_cov = {bam1: {contig1: {window1: [start, end]}}}
	#######################################################################
	# Included custom functions: get_badcov
	#######################################################################
	if gofrom < 3:
		start_time = time()

		bad_cov = {}
		for contig in contig_lengths:
			bad_cov[contig] = []

		bad_cov_temp = {}
		for bam in args.sortbam:
			bad_cov_temp[bam] = {}

		to_badcov = pool.starmap(get_badcov, \
			[(bam, contig_lengths, cov_stats[bam]) for bam in args.sortbam])
		for entry in range(0,len(args.sortbam)):
			bad_cov_temp[args.sortbam[entry]] = to_badcov[entry]

		# Merge overlapping badcov windows
		## Code adapted from https://stackoverflow.com/questions/58535825/combine-overlapping-ranges-of-numbers
		## and https://stackoverflow.com/questions/36955553/sorting-list-of-lists-by-the-first-element-of-each-sub-list

		for contig in contig_lengths:
			for bamfile in bad_cov_temp:
				for window in bad_cov_temp[bamfile][contig]:
					bad_cov[contig].append(window)

		for contig in bad_cov:
			while True:
				if not any_items_overlap(bad_cov[contig]):
					# No items overlapped - break the loop and finish
					break
				else:
					item1, item2 = any_items_overlap(bad_cov[contig])

					# Remove the items from the main list
					bad_cov[contig].remove(item1)
					bad_cov[contig].remove(item2)

					# Replace them with a merged version
					item_values = item1 + item2
					bad_cov[contig].append([min(item_values), max(item_values)])
					# Start the loop again to check for any other overlaps

			bad_cov[contig] = sorted(bad_cov[contig], key=lambda x: x[0])

		# Export bad_cov to pickle
		with open(out_pickle_tmp + "/bad_cov.pickle", "wb") as outfile:
			pickle.dump(bad_cov, outfile)
		outfile.close()

		# Timing - time taken to get bad coverage stats
		print_time("Get irregular coverage stats", start_time)

	else:
		with open(out_pickle_tmp + "/bad_cov.pickle", "rb") as infile:
			bad_cov = pickle.load(infile)
		infile.close()

	#######################################################################
	# STEP 4 - GENERATE A LIST OF POSITIONS WHERE VARIANTS OCCUR
	#	OUT: all_variants = {contig1: [var1, var2, var3]}
	#		all_SNPs = {contig1: [SNP1, SNP2, SNP3]}
	#		all_indels = {contig1: [ind1, ind2, ind3]}
	#######################################################################
	# Included custom functions: get_variants
	#######################################################################
	if gofrom < 4:
		start_time = time()

		print("Searching for potential barcodes in",len(args.sortbam),"file(s).")

		# HOW TO FIX THIS SO IT READS GZIP FILES STRAIGHTAWAY?

		for bam in args.sortbam:
			print("\nFinding variants in " + FileName(bam) + "...")
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

		# Export all_variants, all_SNPs, and all_indels to pickles

		with open(out_pickle_tmp + "/all_variants.pickle", "wb") as outfile1, open(out_pickle_tmp + "/all_SNPs.pickle", "wb") as outfile2, open(out_pickle_tmp + "/all_indels.pickle", "wb") as outfile3:
			pickle.dump(all_variants, outfile1)
			pickle.dump(all_SNPs, outfile2)
			pickle.dump(all_indels, outfile3)
		outfile1.close()
		outfile2.close()
		outfile3.close()

		# Timing - time taken to get lists of variants
		print_time("Get lists of variants", start_time)

	else:
		with open(out_pickle_tmp + "/all_variants.pickle", "rb") as infile1, open(out_pickle_tmp + "/all_SNPs.pickle", "rb") as infile2, open(out_pickle_tmp + "/all_indels.pickle", "rb") as infile3:
			all_variants = pickle.load(infile1)
			all_SNPs = pickle.load(infile2)
			all_indels = pickle.load(infile3)
		infile1.close()
		infile2.close()
		infile3.close()

	#######################################################################
	# STEP 5 - STEP THROUGH EACH CONTIG LOOKING FOR WINDOWS WITH
	#		VARIABLE CENTRES AND CONSERVED PRIMER SITES
	# OUT: master_dict = {contig1: [Window1, Window2, Window3]}
	#######################################################################
	# Included custom functions: find_windows
	#######################################################################
	# DevNote - this needs speeding up!
	#######################################################################
	if gofrom < 5:
		start_time = time()

		print("\nChecking windows...")

		to_master = pool.starmap(find_windows, \
		[(contig, contig_lengths[contig], args.window_size, args.primer_size, all_variants[contig], barcode_log, out_pickle_tmp) for contig in contig_lengths])

		for entry in range(0,len(contig_lengths)):
			master_dict[list(contig_lengths.keys())[entry]] = to_master[entry]

		# Export master_dict to pickle
		with open(out_pickle_tmp + "/master_dict.pickle", "wb") as outfile:
			pickle.dump(master_dict, outfile)
		outfile.close()

		# Delete per-contig find_windows temporary files
		for contig in contig_lengths:
			os.remove(out_pickle_tmp + "/find_windows." + contig + ".pickle")

		# DevNote - trying to save memory space
		all_variants.clear()

		# Timing - time taken to find valid windows
		print_time("Find valid windows", start_time)

	else:
		with open(out_pickle_tmp + "/master_dict.pickle", "rb") as infile:
			master_dict = pickle.load(infile)
		infile.close()

	#######################################################################
	# STEP 6 - EXCLUDE WINDOWS OF LOW COVERAGE
	# OUT: good_cov_dict = {contig1: [Window1, Window2, Window3]}
	#######################################################################
	# Included custom functions: good_coverage
	#######################################################################
	if gofrom < 6:
		start_time = time()

		print("\nChecking window coverage...")

		to_good_cov = pool.starmap(good_coverage, \
			[(master_dict[contig], cov_stats, bad_cov, args.badcov, args.sortbam, contig, barcode_log) for contig in contig_lengths])

		for entry in range(0,len(contig_lengths)):
			good_cov_dict[list(contig_lengths.keys())[entry]] = to_good_cov[entry]

		# Export good_cov_dict to pickle
		with open(out_pickle_tmp + "/good_cov_dict.pickle", "wb") as outfile:
			pickle.dump(good_cov_dict, outfile)
		outfile.close()

		# DevNote - trying to save memory space
		master_dict.clear()

		# Timing - time taken to get good-coverage windows
		print_time("Get good-coverage windows", start_time)

	else:
		with open(out_pickle_tmp + "/good_cov_dict.pickle", "rb") as infile:
			good_cov_dict = pickle.load(infile)
		infile.close()

	#######################################################################
	# STEP 7 - MERGE OVERLAPPING WINDOWS
	# OUT: merged_dict = {contig1: [Window1, Window2, Window3]}
	#######################################################################
	# Included custom functions: merge_windows
	#######################################################################
	if gofrom < 7:
		start_time = time()

		print("\nMerging overlapping windows...")

		for contig in good_cov_dict:
			if good_cov_dict[contig]:
				print(contig)
				merged_dict[contig] = merge_windows(good_cov_dict[contig], contig)

		merged_list = [item for sublist in list(merged_dict.values()) for item in sublist]

		# Export merged_list to pickle
		with open(out_pickle_tmp + "/merged_list.pickle", "wb") as outfile:
			pickle.dump(merged_list, outfile)
		outfile.close()

		# Write number of merged windows to logfile
		with open(barcode_log, 'a') as outlog:
			for contig in merged_dict:
				outlog.write(str(len(merged_dict[contig])) + " merged windows found in contig " + contig + ".\n")
		outlog.close()

		# Timing - time taken to merge windows
		print_time("Merge windows", start_time)

	else:
		with open(out_pickle_tmp + "/merged_list.pickle", "rb") as infile:
			merged_list = pickle.load(infile)
		infile.close()

	#######################################################################
	# STEP 8 - ENSURE THAT EACH SAMPLE CONTAINS AT LEAST ONE UNIQUE ALLELE
	#		AT EACH LOCUS, AND THAT THE COVERAGE IS ACCEPTABLE
	# OUT: really_final_list = [[Window1, min_diffs, max_diffs, p1_seq, var_seq, p2_seq],[Window2, ...]]
	#######################################################################
	# Included custom functions: verify_windows
	#	Subfunctions: compare, get_allele_consensus, check_coverage
	#######################################################################
	# DevNote - find a way to skip loci where the following type of error occurs:
	# `Warning: ignoring overlapping variant starting at 000215F:2953`
	#######################################################################
	# DevNote - this needs speeding up!
	#######################################################################
	if gofrom < 8:
		start_time = time()

		print("\nChecking consensuses...")

		# Make directory for saving alleles
		if not os.path.isdir(out_fasta_dir):
			os.mkdir(out_fasta_dir)

		final_list = []

		# If each chunk would be less than 2 Gb, then split evenly among threads
		# Otherwise, split so that each chunk is < 2 Gb

		chunk_len = int(len(merged_list) / int(args.threads))

		if args.dev:
			print("merged_list pickle: " + str(len(pickle.dumps(merged_list))) + " bytes.\n")
			print("args.ref pickle: " + str(len(pickle.dumps(args.ref))) + " bytes.\n")
			print("args.sortbam pickle: " + str(len(pickle.dumps(args.sortbam))) + " bytes.\n")
			print("cov_stats pickle: " + str(len(pickle.dumps(cov_stats))) + " bytes.\n")
			print("bad_cov pickle " + str(len(pickle.dumps(bad_cov))) + " bytes.\n")

		to_final = pool.starmap(verify_windows, \
			[(chunky, args.ref, args.sortbam, cov_stats, args.ploidy, out_fasta_dir) for chunky in chunks(merged_list, chunk_len)])

		chunk_no = int(args.threads)

		# DevNote - print to_final chunks to file
		for entry in range(0,chunk_no):
#			if not os.isfile(out_pickle_tmp + "/final_list.chunk" + str(entry) + ".pickle"):
#				with open(out_pickle_tmp + "/final_list.chunk" + str(entry) + ".pickle", "wb") as outfile:
#					pickle.dump(to_final[entry], outfile)
#				outfile.close()
#			else:
#				with open(out_pickle_tmp + "/final_list.chunk" + str(entry) + ".pickle", "rb") as infile:
#					to_final[entry] = 

			final_list.append(to_final[entry])

		really_final_list = [item for sublist in final_list for item in sublist]

		# Export really_final_list to pickle
		with open(out_pickle_tmp + "/really_final_list.pickle", "wb") as outfile:
			pickle.dump(really_final_list, outfile)
		outfile.close()

		# Timing - time taken to get unique windows
		print_time("Get unique windows", start_time)

	else:
		with open(out_pickle_tmp + "/really_final_list.pickle", "rb") as infile:
			really_final_list = pickle.load(infile)
		infile.close()

	#######################################################################
	# STEP 9 - REPORT RESULTS IN TXT AND BED FORMAT
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

	# Delete checkpoint files if not in --dev mode
	if not args.dev:
		shutil.rmtree(out_pickle_tmp)

	# Timing - total pipeline time
	total_time = "Total time: " + str(round(time() - full_time, 1)) + " seconds.\n"
	with open("time.log", 'a') as outlog:
		outlog.write("\n")
		outlog.write(total_time)
	outlog.close()

#######################################################################

if __name__ == "__main__":
	main()
