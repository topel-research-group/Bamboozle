#!/usr/bin/env python3


#	Complex and contig flags: Calculate median coverage of contig,
#	and identify regions deviating by +/- 50%; output in bed format
#
#	Complex flag only: Calculate median coverage of each contig
#	in assembly, and identify regions deviating by +/- 50%; output
#	in bed format
#
#	Simple and contig flags: Obtain a per-contig median average
#	coverage for the specified contig
#
#	Simple flag only: Obtain a per-contig median average coverage
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

import os
import subprocess
from statistics import median

#######################################################################
# GENERATE A BED FILE OF RESULTS
#######################################################################

def make_bed(contig_lib,this_contig,openfile):
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
			result1 = this_contig + "\t" + str(FirstHigh - 1) + "\t" + str(LastHigh) + "\t" + "HighCoverage" + "\n"
			openfile.write(result1)
			FirstHigh = 0
			LastHigh = 0
		if contig_lib[key] > lower and LastLow != 0:
			result2 = this_contig + "\t" + str(FirstLow - 1) + "\t" + str(LastLow) + "\t" + "LowCoverage" + "\n"
			openfile.write(result2)
			FirstLow=0
			LastLow=0
	# Print last high event, if contig ends on a high
	if LastHigh != 0:
		final_high = this_contig + "\t" + str(FirstHigh - 1) + "\t" + str(LastHigh) + "\t" + "HighCoverage" + "\n"
		openfile.write(final_high)
	# Print last low event, if contig ends on a low
	if LastLow != 0:
		final_low = this_contig + "\t" + str(FirstLow - 1) + "\t" + str(LastLow) + "\t" + "LowCoverage" + "\n"
		openfile.write(final_low)

#######################################################################
# MAIN
#######################################################################

## DevNote - Shows funky behaviour at stretches hovering around the threshold...

def main(args):

	if not args.outprefix:
		args.outprefix = os.path.basename(args.sortbam[:-4])

	if args.simple:
		outfile = args.outprefix + ".txt"
	elif args.complex:
		outfile = args.outprefix + ".bed"

	with open(outfile, "a") as output_file:

		if args.contig:
			cmd = ["samtools", "depth", "-aa", args.sortbam, "-r", args.contig]
			process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=False)

			cov_stats = {}
			current_contig = args.contig
			with process.stdout as result:
				rows = (line.decode().split('\t') for line in result)
				for row in rows:
					position = int(row[1])
					coverage = int(row[2])
					cov_stats[position] = coverage
			if args.complex:
				bed_header = "track name=Weird Coverage" + "\t" + "description='Areas +/- 50% of median coverage'" + "\n"
				output_file.write(bed_header)
				make_bed(cov_stats,current_contig,output_file)
			elif args.simple:
				txt_record = current_contig + "\t" + str(median(cov_stats.values())) + "\n"
				output_file.write(txt_record)

		else:
			cmd = ["samtools", "depth", "-aa", args.sortbam]
			process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=False)

			cov_stats = {}
			current_contig = "None"
			if args.complex:
				bed_header = "track name=WeirdCoverage" + "\t" + "description='Areas +/- 50% of median coverage'" + "\n"
				output_file.write(bed_header)

			# DevNote = using samtools faidx there may be a more elegant way to report the last contig

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
							make_bed(cov_stats,current_contig,output_file)
						elif args.simple:
							txt_record = current_contig + "\t" + str(median(cov_stats.values())) + "\n"
							output_file.write(txt_record)
						cov_stats = {}
						current_contig = ctg
						cov_stats[position] = coverage
				if args.complex:
					# Print stats for the final contig
					make_bed(cov_stats,current_contig,output_file)
				elif args.simple:
					last_record = current_contig + "\t" + str(median(cov_stats.values())) + "\n"
					output_file.write(last_record)
	output_file.close()

#######################################################################

if __name__ == "__main__":
	main()
