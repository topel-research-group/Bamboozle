#!/usr/bin/env python3


#	Identify regions of 0x coverage in specified contig,
#	then print reference sequence and GC content at these coordinates
#	Fasta parsing from Biopython, fp.py and
#	http://stackoverflow.com/questions/7654971/parsing-a-fasta-file-using-a-generator-python
#
#	Part of bamparser.py in earlier versions of Bamboozle
#
#       Copyright (C) 2018 Matthew Pinder. matt_pinder13@hotmail.com
#
#       This program is free software: you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation, either version 3 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program.  If not, see <https://www.gnu.org/licenses/>.

#######################################################################
# IMPORTS
#######################################################################

import sys
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

def zero_print(reference, contig, zero_dict):
	with open(reference) as fasta:
		for name, seq in read_fasta(fasta):
			if name[1:] == contig:
				print("GC% for contig:",round(get_gc(seq), 3))
				print("Contig\tPositions\tGC%\tSequence")
				for key in zero_dict:
					if key + 1 == (zero_dict[key]):
						print(contig,zero_dict[key],"-",seq[key:zero_dict[key]],sep="\t")
					else:
						zero_range = str(key + 1) + "-" + str(zero_dict[key])
						print(contig,zero_range,round(get_gc(seq[key:zero_dict[key]]), 3),seq[key:zero_dict[key]],sep="\t")
	sys.exit()

#######################################################################
# MAIN
#######################################################################

## DevNote - Add try/except statement for bedtools
## DevNote - How to ensure that bam is sorted?

def main(args):
	if not args.contig or not args.ref or not args.sortbam:
		print("Please ensure all required flags are specified; see readme file")
		sys.exit()

	if args.verbose == True:
		print("Finding zero coverage areas in contig",args.contig)

	cmd = ["bedtools genomecov -bga -ibam %s" % args.sortbam]
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
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
				zero_print(args.ref, args.contig, zeroes)
	zero_print(args.ref, args.contig, zeroes)

#######################################################################

if __name__ == "__main__":
	main()
