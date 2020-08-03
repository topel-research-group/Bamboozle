#!/usr/bin/env python3


#	Scan for potential heterozygous/deletion sites -
#	per base, discrete events, or frameshifts
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

import sys
import subprocess
import os.path

#######################################################################
# PRINT DELETION
#######################################################################

def print_deletion(m, n, version, mutation_list, writeout):
	m.append(n)
	if version == "deletion2" or (version == "deletion3" and n % 3 != 0):
		del_event = str(m[0]) + "\t" + str(m[1]) + "\t" + str(m[2]) + "\n"
		writeout.write(del_event)
	elif version in {"deletionx", "homohetero"}:
		mutation_list.append(m)

#######################################################################
# FIND DELETIONS OCCURRING WITHIN EXONS
#######################################################################

def exon_mutations(exons, mutation_list, verbosity, writeout):

	frameshifts = 0
	exon_list = []

	with open(exons, 'r') as input:
		for line in input:
			if not line.startswith("track name"):
				exon_list.append(line)

	verbose_header = "Contig" + "\t" + "Start" + "\t" + "bp" + "\t" + "Exon" + "\n"
	writeout.write(verbose_header)

	for mut in mutation_list:
		for ex in exon_list:
			ex = ex.split("\t")
			if (mut[0] == ex[0]) and (int(ex[1]) <= int(mut[1]) <= int(ex[2])):
				if int(mut[2]) % 3 != 0:
					frameshifts += 1
				record_out = str(mut[0]) + "\t" + str(mut[1]) + "\t" + str(mut[2]) + "\t" + ex[3] + "\n"
				writeout.write(record_out)
				break

	if verbosity:
		print("Total number of frameshifts in exons:",frameshifts)

#######################################################################
# CALCULATE PERCENTAGE COVERAGE DIFFERENCE BETWEEN FIRST BASE IN A
# MUTATION AND THE BASE BEFORE IT, TO DETERMINE HOMO/HETEROZYGOSITY
#######################################################################

## DevNote - Any way to remove the need for a temporary file?

def HomoDel_or_Hetero(infile, mutation_list, contig, writeout):

	# Print mutation_list to a temporary file

	temp_bed = "temporary_bed_file"

	if os.path.isfile(temp_bed) == True:
		sys.exit("[Error] Temporary file temporary_bed_file already exists.")

	for item in mutation_list:
		with open(temp_bed, "a") as temp_output_file:
			temp_output_file.write(item[0] + "\t" + str(int(item[1]) - 2) + "\t" + str(item[1]) + "\n")
	temp_output_file.close()

	# Compare the coverage

	if contig != "assembly":
		cmd = ["samtools", "depth", "-aa", "-b", temp_bed, "-r", contig, infile]
	else:
		cmd = ["samtools", "depth", "-aa", "-b", temp_bed, infile]

	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=False)

	coverage = []

	with process.stdout as result:
		rows = (line.decode().split('\t') for line in result)
		for row in rows:
			coverage.append(row)

	out_header = "Contig" + "\t" + "Position" + "\t" + "% cov. difference" + "\n"
	writeout.write(out_header)

	while len(coverage) >= 2:
		before = coverage[0]
		after = coverage[1]
		if len(coverage) >= 3:
			next = coverage[2]
		if int(after[1]) - int(before[1]) == 1:
			if (int(before[2]) == 0 and int(after[2]) == 0) or (int(before[2]) == 0):
				out_record = after[0] + "\t" + after[1] + "\t" + "N/A" + "\n"
				writeout.write(out_record)
			else:
				out_record = after[0] + "\t" + after[1] + "\t" + str(round(((100.0/int(before[2]))*int(after[2])), 2)) + "\n"
				writeout.write(out_record)

			del coverage[0]
			if int(next[1]) - int(after[1]) != 1:
				del coverage[0]
		else:
			sys.exit("[Error] Whoops, something went wrong!")

	# Delete the temporary file
	os.remove(temp_bed)

#######################################################################
# MAIN
#######################################################################

## DevNote - currently skips first position

def main(args):

	if not args.outprefix:
		args.outprefix = os.path.basename(args.sortbam[:-4])

	# Determine which variant of the deletions script is being run
	if args.command == "deletion1":
		version = "deletion1"
	elif args.command == "deletion2":
		version = "deletion2"
	elif args.command == "deletion3":
		version = "deletion3"
	elif args.command == "deletionx":
		version = "deletionx"
	elif args.command == "homohetero":
		version = "homohetero"

	output_file = args.outprefix + "." + version + ".txt"

	mutation_list = []

	if args.verbose == True:
		if not args.contig:
			args.contig = "assembly"
		if args.command == "deletion1":
			print("Finding individual deletions in",args.contig)
		elif args.command == "deletion2":
			print("Finding deletion events in",args.contig)
		elif args.command == "deletion3":
			print("Finding frameshift deletion events in",args.contig)
		elif args.command == "deletionx":
			print("Finding deletion events within exons in",args.contig)
		elif args.command == "homohetero":
			print("Determining potentially homozygous/heterozygous deletions in",args.contig)

	if args.contig and args.contig != "assembly":
		cmd = ["samtools", "depth", "-aa", args.sortbam, "-r", args.contig]
	else:
		cmd = ["samtools", "depth", "-aa", args.sortbam]

	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=False)
	
	old_position = 0
	deletion = []
	del_size = 1

	with open(output_file, 'a') as outfile:
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
								if args.command == "deletion1":
									del_pos = ctg + "\t" + str(x) + "\n"
									outfile.write(del_pos)
								else:
									if (int(x) - int(old_position)) != 1:
										if len(deletion) != 0:
											print_deletion(deletion, del_size, version, mutation_list, outfile)
											deletion = []
										deletion.extend([ctg,x])
										del_size = 1
									else:
										del_size +=1
								old_position = x
				else:
					window[position] = coverage

			# Ensure that the final event is reported
			if args.command != "deletion1":
				print_deletion(deletion, del_size, version, mutation_list, outfile)

		if args.command == "deletionx":
			exon_mutations(args.exons, mutation_list, args.verbose, outfile)

		if args.command == "homohetero":
			HomoDel_or_Hetero(args.sortbam, mutation_list, args.contig, outfile)

	outfile.close()

#######################################################################

if __name__ == "__main__":
	main()
