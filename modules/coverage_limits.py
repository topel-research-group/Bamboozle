#!/usr/bin/env python3


#	Identify the longest continuous region of a contig where
#	all positions fall between defined coverage limits
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

import subprocess

#######################################################################
# MAIN
#######################################################################

def main(args):

	if not args.outprefix:
		args.outprefix = os.path.basename(args.sortbam[:-4])
	output_file = args.outprefix + ".cov_limits.txt"

	command = ["samtools", "depth", "-aa", args.sortbam, "-r", args.contig]
	process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=False)

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

	results = "Longest stretch between " + str(lower) + "x and " + str(upper) + "x coverage on " + args.contig + ": " + str(longest) + "bp" + "\n" + \
			args.contig + ":" + str(start) + "-" + str(stop) + "\n"

	if args.verbose:
		print(results)

	with open(output_file, 'a') as result_file:
		result_file.write(results)
	result_file.close()

#######################################################################

if __name__ == "__main__":
	main()
