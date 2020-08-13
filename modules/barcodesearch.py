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
from itertools import combinations
from Levenshtein import distance
from multiprocessing import Pool

from functools import wraps
from time import time
import datetime

#######################################################################
# TIME DECORATOR (TAKEN FROM VILMA'S PIPELINE
#######################################################################

def timing(function):
        @wraps(function)
        def wrapper(*args, **kwargs):
                now = datetime.datetime.now()
                start = time()
                result = function(*args, **kwargs)
                end = time()
                fh = open("time.log", "a")
                lines_of_text = now.strftime("%Y-%m-%d %H:%M") \
                                + ' Function: ' \
                                + function.__name__ \
                                + ' Elapsed time: {}'.format(end-start) \
                                + ' seconds \n'
                fh.writelines(lines_of_text)
                fh.close()
                return result
        return wrapper

#######################################################################
# GET REASONABLE SAMPLE NAMES
#######################################################################

def FileName(long_name):
	return os.path.splitext(os.path.basename(long_name))[0]

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
# PARSE BAM FILES USING BCFTOOLS
#######################################################################

def bcf(infile, contig_list, quality, threads, reference):
	print("\n" + FileName(infile))

	if os.path.isfile(FileName(infile) + ".vcf.gz") == True:
		import gzip
		print("VCF file already exists for",FileName(infile),"- reading file...")
		process2 = gzip.open(FileName(infile) + ".vcf.gz", 'rt')
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
# ENSURE UNIQUENESS OF VARIABLE REGIONS BETWEEN STRAINS
# DevNote - this needs speeding up!
#######################################################################

def check_unique_windows(windows, contig, reference, infiles):
	final = {}

	for window in windows[contig]:

		compare_seqs = {}
		for bam in infiles:
			compare_seqs[FileName(bam)] = ""
			vcf_zipped = FileName(bam) + ".vcf.gz"

			# In the interest of identifying reliable barcodes, only homozygous sites are considered

			con_range = contig + ":" + str(windows[contig][window][0]) + "-" + str(windows[contig][window][3])
			cmdD = ["samtools", "faidx", reference, con_range]
			procD =  subprocess.Popen(cmdD, stdout=subprocess.PIPE, shell=False)

			cmdE = ["bcftools", "consensus", "-i", "GT=\"hom\"", "--sample", bam, vcf_zipped]
			process3 = subprocess.Popen(cmdE, stdin=procD.stdout, stdout=subprocess.PIPE, shell=False)

			with process3.stdout as result3:
				rows3 = (line.decode() for line in result3)
				for row3 in rows3:
					compare_seqs[FileName(bam)] += (row3.upper().strip("\n"))

		# Ensure all variable regions are unique between samples
		if (len(compare_seqs) == len(set(compare_seqs.values()))):

			final[window] = windows[contig][window]

			# Calculate the number of differences between each pair of sequences
			# DevNote - this step is taking the longest in this function

			diff_counts = []
			list1 = list(compare_seqs.values())
			for i, j in list(combinations(range(0,len(list1)), 2)):
				diff_counts.append(distance(list1[i],list1[j]))

			# Add minimum and maximum differences to each window
			final[window].append(min(diff_counts))
			final[window].append(max(diff_counts))

			# Get the required sequences for conserved and variable regions

			window_start = windows[contig][window][0]
			primer1_stop = windows[contig][window][1]
			primer2_start = windows[contig][window][2] + 1
			window_stop = windows[contig][window][3]

			window_range = contig + ":" + str(window_start) + "-" + str(window_stop)

			full_window = ""

			cmd4 = ["samtools faidx %s %s" % (reference, window_range)]
			process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, shell=True)
#			process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, shell=False)
			with process4.stdout as result1:
				rows1 = (line.decode() for line in result1)
				for row1 in rows1:
					if not row1.startswith(">"):
						full_window += row1.strip("\n")
			final[window].append(full_window[(window_start - window_start):(primer1_stop - window_start)])
			final[window].append(full_window[(primer1_stop - window_start):(primer2_start - window_start)])
			final[window].append(full_window[(primer2_start - window_start):((window_stop - window_start)+1)])

	return(final)

#######################################################################
# PRINT DURATION OF THE STEP
#######################################################################

def print_time(step_name, start_time):
	running_time = step_name + ": " + str(time() - start_time) + " seconds.\n"
	with open("time.log", 'a') as outlog:
		outlog.write(running_time)
	outlog.close()

#######################################################################
# MAIN
#######################################################################

@timing
def main(args):

	# Timing - get start time
	if args.dev:
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
	if args.dev:
		print_time("Get contig lengths", start_time)
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
		vcf_file = FileName(bam) + ".vcf"

		# Generate or read in a VCF file for the current BAM
		process2 = bcf(bam, contig_lengths, filter_qual, args.threads, args.ref)

		# Generate a list of loci where variants occur
		print("\nFinding variants...")
		current_contig = ""
		if os.path.isfile(FileName(bam) + ".vcf.gz") == True:
			if os.path.isfile(FileName(bam) + ".vcf.gz.csi") == False:
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
	if args.dev:
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
	if args.dev:
		print_time("Find valid windows", start_time)
		start_time = time()

	# Merge overlapping windows
	print("\nMerging overlapping windows...")

	for contig in master_dict:
		if master_dict[contig]:
			print(contig)
			merged_dict[contig] = merge_windows(contig, master_dict)

	# Timing - time taken to merge windows
	if args.dev:
		print_time("Merge windows", start_time)
		start_time = time()

# Dev Note: Find a way to skip loci where the following type of error occurs:
# `Warning: ignoring overlapping variant starting at 000215F:2953`

	# Compare consensus sequences for all BAMs, to ensure they are truly unique,
	# not merely differing from the reference in all the same positions
	# DevNote - this needs speeding up!
	print("\nChecking consensuses...")

	to_final = pool.starmap(check_unique_windows, \
		[(merged_dict, contig, args.ref, args.sortbam) for contig in contig_lengths])

	for entry in range(0,len(contig_lengths)):
		final_dict[list(contig_lengths.keys())[entry]] = to_final[entry]

	# Timing - time taken to get unique windows
	if args.dev:
		print_time("Get unique windows", start_time)
		start_time = time()

	# Report results in TXT and BED format
	print("\nGenerating output files...")

	with open(out_bed, "a") as output_bed, open(out_txt, "a") as output_txt:
		output_bed.write("track name=PotentialBarcodes description=\"Potential barcodes\"\n")

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

				output_bed.write(window_out)
				output_txt.write(line_out)
		output_bed.close()
		output_txt.close()

	# Timing - time taken to write output
	if args.dev:
		print_time("Write output", start_time)

#######################################################################

if __name__ == "__main__":
	main()
