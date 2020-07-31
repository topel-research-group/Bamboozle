#!/usr/bin/env python3


#	Identify regions of 0x coverage in specified contig,
#	then print reference sequence and GC content at these coordinates
#	Fasta parsing from Biopython, fp.py and
#	http://stackoverflow.com/questions/7654971/parsing-a-fasta-file-using-a-generator-python
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
import os
import subprocess

#######################################################################
# READ THE FASTA FILE
#######################################################################

def read_fasta(fasta):
	name, seq = None, []
	for line in fasta:
		line = line.rstrip()
		if line.startswith(">"):
			if name: yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))

#######################################################################
# GET THE GC CONTENT
#######################################################################

def get_gc(input):
	count = 0
	gc_list = ['G', 'C', 'g', 'c']
	for base in input:
		if base in gc_list:
			count += 1
	gc_content = 100.0 / len(input) * count
	return gc_content

#######################################################################
# PRINT RESULTS
#######################################################################

def zero_print(reference, contig, zero_dict, outfile):
	with open(reference) as fasta, open(outfile, 'a') as output:
		for name, seq in read_fasta(fasta):
			if name[1:] == contig:
				print("GC% for contig:",round(get_gc(seq), 3))
				output.write("Contig\tPositions\tGC%\tSequence\n")
				for key in zero_dict:
					if key + 1 == (zero_dict[key]):
						write_me = contig + "\t" + str(zero_dict[key]) + "\t-\t" + seq[key:zero_dict[key]] + "\n"
					else:
						zero_range = str(key + 1) + "-" + str(zero_dict[key])
						write_me = contig + "\t" + zero_range + "\t" + str(round(get_gc(seq[key:zero_dict[key]]), 3)) + "\t" + seq[key:zero_dict[key]] + "\n"
					output.write(write_me)
	output.close()
	sys.exit()

#######################################################################
# MAIN
#######################################################################

def main(args):
	if args.verbose == True:
		print("Finding zero coverage areas in contig",args.contig)

	if not args.outprefix:
		args.outprefix = os.path.basename(args.sortbam[:-4])
	output_file = args.outprefix + ".zero_regions.txt"

	cmd = ["bedtools", "genomecov", "-bga", "-ibam", args.sortbam]
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=False)
	zeroes = {}
	correct_contig = 0

	with process.stdout as result:
		rows = (line.decode().split('\t') for line in result)
		for row in rows:
			if str(row[0]) == args.contig:
				correct_contig = 1
				coverage = int(row[3])
				if coverage == 0:
					zeroes[int(row[1])] = int(row[2])
			elif str(row[0]) != args.contig and correct_contig == 1:
				zero_print(args.ref, args.contig, zeroes, output_file)
	zero_print(args.ref, args.contig, zeroes, output_file)

#######################################################################

if __name__ == "__main__":
	main()
