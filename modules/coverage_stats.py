#!/usr/bin/env python3


#	Calculate percentage of positions in assembly/contig with
#	read coverage >= a given threshold (default: 20x).
#
#	If output file (and optional GFF) specified, outputs a BED file
#	defining regions with coverage about the threshold (and any
#	gene models overlapping these regions)
#
#	Part of bamparser.py in earlier versions of Bamboozle
#
#	Copyright (C) 2018 Matthew Pinder. matt_pinder13@hotmail.com
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
import os.path

#######################################################################
# PRINT WINDOW TO BED FILE
#######################################################################

#def print_to_bed(start_coord, stop_coord, contig, gff_file, output_file):
def print_to_bed(start_coord, stop_coord, contig, gff_loci, output_file):
	print_genes = []

	# Determine whether the window overlaps with a feature in the gff file
	# DevNote - Any way to make this less wordy?
#	if gff_file:
	if gff_loci:
		if contig in gff_loci.keys():
			for gene in gff_loci[contig]:
				if start_coord <= int(gff_loci[contig][gene][0]) <= stop_coord \
				or start_coord <= int(gff_loci[contig][gene][1]) <= stop_coord \
				or int(gff_loci[contig][gene][0]) <= start_coord <= int(gff_loci[contig][gene][1]) \
				or int(gff_loci[contig][gene][0]) <= stop_coord <= int(gff_loci[contig][gene][1]):
					print_genes.append(gene)

	# Define BED file entry depending on whether any genes were found
	if print_genes:
		output_line = contig + "\t" + str(start_coord - 1) + "\t" + str(stop_coord) + "\t" + ','.join(print_genes) + "\n"
	else:
		output_line = contig + "\t" + str(start_coord - 1) + "\t" + str(stop_coord) + "\n"

	# Print to output file
	with open(output_file, "a") as output_bed:
		output_bed.write(output_line)
	output_bed.close()

#######################################################################
# MAIN
#######################################################################

def main(args):

	# If an output file of the desired name already exists, print warning and quit
	if args.outprefix:
		outfile = args.outprefix + ".bed"
		if os.path.isfile(outfile) == True:
			print("Warning: Output file",outfile,"already exists. Please choose another output prefix.")
			sys.exit(1)

	# Generate a dictionary of dictionaries of lists; contig -> gene name -> coordinates
	in_gff_loci = {}
	if args.gff:
		if args.verbose:
			print("Parsing GFF file...")

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
		with open(outfile, "a") as output_bed:
			output_bed.write("track name=Coverage description=\"Coverage above " + str(args.threshold) + "x\"\n")
		output_bed.close()

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
						print_to_bed(start_coord, (position - 1), current_contig, in_gff_loci, outfile)
					recording = False

				# For cases where recorded window continues to end of contig
				if args.outprefix and (position == contig_lengths[current_contig]):
					if recording:
						print_to_bed(start_coord, position, current_contig, in_gff_loci, outfile)
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

if __name__ == "__main__":
	main()
