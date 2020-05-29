#!/usr/bin/env python3


#	Pipeline for retrieving coverage-related statistics from BAM files.
#
#	Copyright (C) 2018 Matthew Pinder. matt_pinder13@hotmail.com
#	Copyright (C) 2018 Mats Topel. mats.topel@marine.gu.se
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

import re
import sys
import subprocess
import argparse
import os.path
from statistics import median

#######################################################################
# COVERAGE STATS
#	Calculate percentage of positions in assembly/contig with
#	read coverage >= a given threshold (default: 20x)
#	If output file (and optional GFF) specified, outputs a BED file
#	defining regions with coverage about the threshold (and any
#	gene models overlapping these regions)
#######################################################################

def coverage_stats(args):

	# If an output file of the desired name already exists, print warning and quit
	if args.outprefix:
		outfile = args.outprefix + ".bed"

	if args.outprefix and os.path.isfile(outfile) == True:
		print("The specified output file already exists; please adjust output prefix name [-o] or delete existing file.")
		sys.exit()

	# Generate a dictionary of dictionaries of lists; contig -> gene name -> coordinates
	if args.gff:
		print("Parsing GFF file...")
		in_gff_loci = {}

		with open(args.gff, 'r') as input:
			lines = (entry.split("\t") for entry in input)
			for line in lines:

				# DevNote - Is this standard in gff files?
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

	# Function for printing window to BED file
	def print_to_bed(start_coord, stop_coord, output_file):
		print_genes = []

		# Determine whether the window overlaps with a feature in the gff file
		# DevNote - Any way to make this less wordy?
		if args.gff:
			if contig in in_gff_loci.keys():
				for gene in in_gff_loci[contig]:
					if start_coord <= int(in_gff_loci[contig][gene][0]) <= stop_coord \
					or start_coord <= int(in_gff_loci[contig][gene][1]) <= stop_coord \
					or int(in_gff_loci[contig][gene][0]) <= start_coord <= int(in_gff_loci[contig][gene][1]) \
					or int(in_gff_loci[contig][gene][0]) <= stop_coord <= int(in_gff_loci[contig][gene][1]):
						print_genes.append(gene)

		# Define BED file entry depending on whether any genes were found
		if print_genes:
			output_line = current_contig + "\t" + str(start_coord - 1) + "\t" + str(stop_coord) + "\t" + ','.join(print_genes) + "\n"
		else:
			output_line = current_contig + "\t" + str(start_coord - 1) + "\t" + str(stop_coord) + "\n"

		# Print to output file
		with open(outfile, "a") as output_file:
			output_file.write(output_line)
		output_file.close()

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

	# Run samtools depth to obtain required coverage info

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

	if args.outprefix:
		with open(outfile, "a") as output_file:
			output_file.write("track name=Coverage description=\"Coverage above " + str(args.threshold) + "x\"\n")
		output_file.close()

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

	with process.stdout as result:
		rows = (line.decode().split('\t') for line in result)
		for row in rows:
			contig = row[0]
			position = int(row[1])
			coverage = int(row[2])

			if not current_contig:
				current_contig = contig

			if contig == current_contig:
				# Start recording
				if coverage >= args.threshold:
					cov_stats[contig] += 1
					if args.outprefix and not recording:
						start_coord = position
						recording = True

				# Print window and stop recording
				elif args.outprefix and (coverage < args.threshold):
					if recording:
						print_to_bed(start_coord, (position - 1), output_file)
					recording = False

				# For cases where recorded window continues to end of contig
				if args.outprefix and (position == contig_lengths[current_contig]):
					if recording:
						print_to_bed(start_coord, position, output_file)
					recording = False

			# New contig
			elif contig != current_contig:
				current_contig = contig

				if coverage >= args.threshold:
					cov_stats[contig] += 1
					if args.outprefix and not recording:
						start_coord = position
						recording = True

	# Print percentage stats to standard out
	if args.contig:
		print("Length of " + args.contig + ":",contig_lengths[args.contig])
		value = 100.0 / contig_lengths[args.contig] * cov_stats[args.contig]
		print(round(value, 3),"% of contig ", args.contig, " with >=",args.threshold,"x coverage.",sep="")
	else:
		print("Length of assembly:",str(assembly_length))
		value = 100.0 / assembly_length * sum(cov_stats.values())
		print(round(value, 3),"% of assembly with >=",args.threshold,"x coverage.",sep="")

	if args.verbose and args.outprefix:
		print("Results printed to",outfile + ".\n")

#######################################################################
# DELETION
#	Scan for potential heterozygous/deletion sites -
#	per base, discrete events, or frameshifts
#######################################################################

## DevNote - currently skips first position

def deletion(args):

	mutation_list = []

	def print_deletion(m, n):
		m.append(n)
		if args.deletion2 or (args.deletion3 and n % 3 != 0):
			print(m[0],m[1],m[2],sep="\t")
		elif args.deletionx or args.homohetero:
			mutation_list.append(m)

	if args.verbose == True:
		if not args.contig:
			args.contig = "assembly"
		if args.deletion1:
			print("Finding individual deletions in",args.contig)
		elif args.deletion2:
			print("Finding deletion events in",args.contig)
		elif args.deletion3:
			print("Finding frameshift deletion events in",args.contig)
		elif args.deletionx:
			print("Finding deletion events within exons in",args.contig)
		elif args.homohetero:
			print("Determining potentially homozygous/heterozygous deletions in",args.contig)

	if args.contig and args.contig != "assembly":
		cmd = ["samtools depth -aa %s -r %s" % (args.sortbam, args.contig)]
	else:
		cmd = ["samtools depth -aa %s" % args.sortbam]

	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
	
	old_position = 0
	deletion = []
	del_size = 1

	with process.stdout as result:
		rows = (line.decode().split('\t') for line in result)
		ctg = ""
		previous_ctg = ""
		for row in rows:
			position = int(row[1])
			coverage = int(row[2])
			if ctg != str(row[0]):
				previous_ctg = ctg
				ctg = str(row[0])
				window = {}
				reported = []

			if len(window) == 12:
				del window[position - 12]
				window[position] = coverage
				base1 = window[position - 11]
				if ((base1*0.8) <= window[position] <= (base1*1.25)) and (base1 > 0) and (window[position] >= args.threshold):
					for x, y in window.items():
						if y < (base1*0.6) and x not in reported:
							reported.append(x)
							if args.deletion1:
								print(ctg,x,sep="\t")
							else:
								if (int(x) - int(old_position)) != 1:
									if len(deletion) != 0:
										print_deletion(deletion, del_size)
										deletion = []
									deletion.extend([ctg,x])
									del_size = 1
								else:
									del_size +=1
							old_position = x
			else:
				window[position] = coverage

		if not args.deletion1:	# Ensure that the final event is reported
			print_deletion(deletion, del_size)

	if args.deletionx:
		exon_mutations(mutation_list)

	if args.homohetero:
		HomoDel_or_Hetero(mutation_list)

#######################################################################
# EXON MUTATIONS
#	Find deletions occurring within exons
#######################################################################

def exon_mutations(mutation_list):

	frameshifts = 0
	exon_list = []

	def list_append(argument, list):
		with open(argument, 'r') as input:
			for line in input:
				list.append(line)

	list_append(args.exons, exon_list)

	if not args.verbose:
		print("Contig","Start","bp","Exon",sep="\t")

	for mut in mutation_list:
		for ex in exon_list:
			ex = ex.split("\t")
			if (mut[0] == ex[0]) and (int(ex[1]) <= int(mut[1]) <= int(ex[2])):
				if int(mut[2]) % 3 != 0:
					frameshifts += 1
				if args.verbose:
					print(mut[2],"bp mutation at",mut[0],mut[1],"hits exon",ex[3])
				else:
					print(mut[0],mut[1],mut[2],ex[3],sep="\t")
				break

	print("Total number of frameshifts in exons:",frameshifts)

#######################################################################
# HOMOZYGOUS DELETION OR HETEROZYGOTE
#	Calculate percentage coverage difference between first base in a mutation
#	and the base before it, to determine homo/heterozygosity
#######################################################################

## DevNote - Any way to remove the need for a temporary file?

def HomoDel_or_Hetero(mutation_list):

	# Read in each line of args.mutations, then print out an altered version to a temporary .bed file

	temp_bed = "temporary_bed_file"

	if os.path.isfile(temp_bed) == True:
		print("Temporary file can't be written; please ensure",temp_bed,"is not a file.")
		sys.exit()

	for item in mutation_list:
		with open(temp_bed, "a") as output_file:
			output_file.write(item[0] + "\t" + str(int(item[1]) - 2) + "\t" + str(item[1]) + "\n")
	output_file.close()

	# Compare the coverage 

	if args.contig:
		cmd = ["samtools depth -aa -b %s -r %s %s" % (temp_bed, args.contig, args.sortbam)]
	else:
		cmd = ["samtools depth -aa -b %s %s" % (temp_bed, args.sortbam)]

	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

	coverage = []

	with process.stdout as result:
		rows = (line.decode().split('\t') for line in result)
		for row in rows:
			coverage.append(row)

	print("Contig","Position","% cov. difference", sep="\t")

	while len(coverage) >= 2:
		before = coverage[0]
		after = coverage[1]
		if len(coverage) >= 3:
			next = coverage[2]
		if int(after[1]) - int(before[1]) == 1:
			if (int(before[2]) == 0 and int(after[2]) == 0) or (int(before[2]) == 0):
				print(after[0],after[1],"N/A",sep="\t")
			else:
				print(after[0],after[1],round(((100.0/int(before[2]))*int(after[2])), 2),sep="\t")

			del coverage[0]
			if int(next[1]) - int(after[1]) != 1:
				del coverage[0]
		else:
			print("Whoops, something went wrong!")
			sys.exit()

	os.remove(temp_bed)	# Delete the temporary file

#######################################################################
# MEDIAN DEVIATION
#	Complex and contig flags: Calculate median coverage of contig,
#	and identify regions deviating by +/- 50%; output in bed format
#
#	Complex flag only: Calculate median coverage of each contig in assembly,
#	and identify regions deviating by +/- 50%; output in bed format
#
#	Simple and contig flags: Obtain a per-contig median average coverage for the specified contig
#
#	Simple flag only: Obtain a per-contig median average coverage
#######################################################################

## DevNote - Shows funky behaviour at stretches hovering around the threshold...

def median_deviation(args):

	def make_bed(contig_lib,this_contig):	# Generate a bed file of results
		median_cov = median(contig_lib.values())
		lower = median_cov * 0.5
		upper = median_cov * 1.5

		FirstHigh = 0
		LastHigh = 0
		FirstLow = 0
		LastLow = 0

		for key in contig_lib:
			if contig_lib[key] > upper and FirstHigh == 0:
				FirstHigh = key
			if contig_lib[key] > upper and FirstHigh != 0:
				LastHigh = key
			if contig_lib[key] < lower and FirstLow == 0:
				FirstLow = key
			if contig_lib[key] < lower and FirstLow != 0:
				LastLow = key
			if contig_lib[key] < upper and LastHigh != 0:
				print(this_contig,FirstHigh - 1,LastHigh,"HighCoverage",sep="\t")
				FirstHigh = 0
				LastHigh = 0
			if contig_lib[key] > lower and LastLow != 0:
				print(this_contig,FirstLow - 1,LastLow,"LowCoverage",sep="\t")
				FirstLow=0
				LastLow=0
		# Print last high event, if contig ends on a high
		if LastHigh != 0:
			print(this_contig,FirstHigh - 1,LastHigh,"HighCoverage",sep="\t")
		# Print last low event, if contig ends on a low
		if LastLow != 0:
			print(this_contig,FirstLow - 1,LastLow,"LowCoverage",sep="\t")

	if args.contig:

		cmd = ["samtools depth -aa %s -r %s" % (args.sortbam, args.contig)]
		process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

		cov_stats = {}
		current_contig = args.contig
		with process.stdout as result:
			rows = (line.decode().split('\t') for line in result)
			for row in rows:
				position = int(row[1])
				coverage = int(row[2])
				cov_stats[position] = coverage
		if args.complex:
			print("track name=WeirdCoverage","description='Areas +/- 50% of median coverage'",sep="\t")
			make_bed(cov_stats,current_contig)
		elif args.simple:
			print(current_contig,median(cov_stats.values()),sep="\t")

	else:

		cmd = ["samtools depth -aa %s" % args.sortbam]
		process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

		cov_stats = {}
		current_contig = "None"
		if args.complex:
			print("track name=WeirdCoverage","description='Areas +/- 50% of median coverage'",sep="\t")

		with process.stdout as result:
			rows = (line.decode().split('\t') for line in result)
			for row in rows:
				ctg = str(row[0])
				position = int(row[1])
				coverage = int(row[2])
				if current_contig == "None":
					current_contig = ctg
				if ctg == current_contig:
					cov_stats[position] = coverage
				elif ctg != current_contig:
					if args.complex:
						make_bed(cov_stats,current_contig)
					elif args.simple:
						print(current_contig,median(cov_stats.values()),sep="\t")
					cov_stats = {}
					current_contig = ctg
					cov_stats[position] = coverage

			if args.complex:
				make_bed(cov_stats,current_contig)	# Print stats for the final contig
			elif args.simple:
				print(current_contig,median(cov_stats.values()),sep="\t")

#######################################################################
# COVERAGE LIMITS
#	Identify the longest continuous region of a contig where
#	all positions fall between defined coverage limits
#######################################################################

def coverage_limits(args):

	command = ["samtools depth -aa %s -r %s" % (args.sortbam, args.contig)]
	process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

	if args.limits[0] < args.limits[1]:
		lower = args.limits[0]
		upper = args.limits[1]
	else:
		lower = args.limits[1]
		upper = args.limits[0]

	longest = 0
	start = "N/A"
	stop = "N/A"

	current = 0
	current_start = "N/A"
	current_stop = "N/A"

	recording = False

	with process.stdout as result:
		rows = (line.decode().split('\t') for line in result)
		for row in rows:
			position = int(row[1])
			coverage = int(row[2])
			if lower <= coverage <= upper:
				if not recording:
					current_start = position
					recording = True
				current += 1
			else:
				current_stop = position - 1
				if current > longest:
					longest = current
					start = current_start
					stop = current_stop
				current = 0
				recording = False

	print("Longest stretch between " + str(lower) + "x and " + str(upper) + \
	"x coverage on " + args.contig + "\t" + str(longest) + "\t" + str(start) + "-" + str(stop))

#######################################################################
# EXTRACT SEQUENCE
#	This function extracts the sequence of the mapped reads from
#	a part of the reference sequence specified by args.range
#######################################################################

## DevNote - This function needs fixing

def extract_sequence(args):
	command = ("samtools mpileup -uf %s %s -r %s:%s | bcftools view -cg -") \
	% (args.ref, args.sortbam, args.contig, args.range)
	bam = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
	header = ">" + args.contig + ":" + args.range
	seq = ""

	for line in bam.stdout:
		row = line.decode(encoding="utf-8", errors="strict")
		if row[0] == "#":
			pass
		else:
			if row.split("\t")[4] == ".":
				seq += row.split("\t")[3]
			else:
				seq += row.split("\t")[4]
		print(header)
		print(seq)
