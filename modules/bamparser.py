#!/usr/bin/env python3


#	Pipeline for retrieving coverage-related statistics from BAM files.
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


import sys
import subprocess
import argparse
import os.path
from statistics import median

#######################################################################

parser = argparse.ArgumentParser(prog="BamParser")
parser.add_argument("-f", "--ref", \
                        help="Reference")
parser.add_argument("-F", "--forward", \
                        nargs='*', \
                        help="Forward reads")
parser.add_argument("-R", "--reverse", \
                        nargs='*', \
                        help="Reverse reads")
parser.add_argument("-b", "--bamfile", \
                        help="BAM infile")
parser.add_argument("--sortbam", \
                        help="Sorted BAM infile")
parser.add_argument("--gff", \
                        help="gff infile")
parser.add_argument("--contigsizes", \
                        help="Contig sizes for gff parser")
parser.add_argument("--feature", \
                        help="Feature for gff parser")
parser.add_argument("-t", "--threads", \
                        default=1, \
                        help="Threads")
parser.add_argument("-e", "--snpeff", \
                        nargs='*', \
                        help="Input options for snpeff, without the '-' before")
parser.add_argument("-s", "--snpsift", \
                        action="store_true", \
                        help="Run snpSift")
parser.add_argument("-r", "--clean", \
                        action="store_true", \
                        help="Removes the SAM and BAM files")
parser.add_argument("-p", "--done", \
                        action="store_true", \
                        help="Add an empty file to mark the directory as done")

#######################################################################
parser.add_argument('--bamparse', action="store_true", help ="Run bamparser")
parser.add_argument('--coverage', action="store_true", help='Print a statistic for what percentage of bases in an assembly have >=Nx coverage')
parser.add_argument('--consensus', action="store_true", help='Extract the consensus sequence of aligned reads from a specific region of the reference sequence (WIP)')
parser.add_argument('--zero', action="store_true", help='Find areas of zero coverage and print the reference sequence, along with a GC percentage')
parser.add_argument('--deletion1', action="store_true", help='Find deletions')
parser.add_argument('--deletion2', action="store_true", help='Find deletion events')
parser.add_argument('--deletion3', action="store_true", help='Find frameshift deletion events')
parser.add_argument('--deletionx', action="store_true", help='Find deletions occurring within exons')
parser.add_argument('--homohetero', action="store_true", help='Attempt to determine whether a deletion is homozygous or heterozygous')
parser.add_argument('--median', action="store_true", help='Find regions differing from contig median by +/- 50%%, or just contig medians')
parser.add_argument('--long_coverage', action="store_true", help='Find the longest region between given coverage limits for a given contig')

parser.add_argument('--complex', action="store_true", help='Print full bed output for median')
parser.add_argument('--simple', action="store_true", help='Print median coverage only for median')
parser.add_argument('-c', '--contig', help='Gives per-contig coverage stats')
parser.add_argument('-d', '--threshold', type=int, nargs='?', const='1', default='20', help='Threshold for calculating coverage percentage; default 20')
parser.add_argument("-a", "--range", help="somethingsomsing")
parser.add_argument("-m", "--mutations", help="List of mutation events; currently requires output from bamboozle deletion function")
parser.add_argument("-x", "--exons", help="Bed file containing exon coordinates (0-based); -m also required")
parser.add_argument("-l", "--limits", type=int, nargs=2, help="Specify lower and upper limits for long_coverage function; two arguments required")
parser.add_argument("-v", "--verbose", action="store_true", help="Be more verbose")
parser.add_argument('--dev', help=argparse.SUPPRESS, action="store_true")

args = parser.parse_args()

if args.feature and args.gff is None:
        parser.error("--feature requires --gff")
elif args.gff and args.feature is None:
        parser.error("--feature requires --gff")

#######################################################################

# Ensure correct version of samtools
def check_samtools():
	try:
		subprocess.check_output('samtools depth 2>&1 | grep -- "-aa"', stderr=subprocess.PIPE, shell=True)
	except subprocess.CalledProcessError:
		print("This version of samtools does not support the `depth -aa` option; please update samtools.")
		exit()

#######################################################################

sortbam = args.sortbam

#######################################################################

# Calculate percentage of positions in assembly/contig with read coverage >= a given threshold (default: 20x)
def coverage_stats(sortbam):

	check_samtools()

	if args.contig:
		cmd = ["samtools depth -aa %s -r %s" % (sortbam, args.contig)]
		sequence = "contig"
	else:
		cmd = ["samtools depth -aa %s" % sortbam]
		sequence = "assembly"

	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

	if args.verbose == True:
		if args.contig:
			print("Obtaining stats for ",args.contig," in ",os.path.basename(sortbam),"; coverage >+",args.threshold,"x.",sep="")
		else:
			print("Obtaining whole-genome stats for ",os.path.basename(sortbam),"; coverage >+",args.threshold,"x.",sep="")
	cov_stats = {}
	num_lines = 0
	with process.stdout as result:
		rows = (line.decode().split('\t') for line in result)
		for row in rows:
			coverage = int(row[2])
			num_lines += 1
			if coverage >= args.threshold:
				if coverage in cov_stats:
					cov_stats[coverage] += 1
				else:
					cov_stats[coverage] = 1
		print("Length of " + sequence + ":",num_lines)
		value = 100.0 / num_lines * sum(cov_stats.values())
		print(round(value, 3),"% of the " + sequence + " with >=",args.threshold,"x coverage.",sep="")


# Identify regions of 0x coverage in specified contig, then print reference sequence and GC content at these coordinates
# Fasta parsing from Biopython, fp.py and http://stackoverflow.com/questions/7654971/parsing-a-fasta-file-using-a-generator-python
## DevNote - Add try/except statement for bedtools
## DevNote - How to ensure that bam is sorted?
def zero_regions(sortbam):

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

	def get_gc(input):
		count = 0
		gc_list = ['G', 'C', 'g', 'c']
		for base in input:
			if base in gc_list:
				count += 1
		gc_content = 100.0 / len(input) * count
		return gc_content

	def zero_print():
		with open(args.ref) as fasta:
			for name, seq in read_fasta(fasta):
				if name[1:] == args.contig:
					print("GC% for contig:",round(get_gc(seq), 3))
					print("Contig\tPositions\tGC%\tSequence")
					for key in zeroes:
						if key + 1 == (zeroes[key]):
							print(args.contig,zeroes[key],"-",seq[key:zeroes[key]],sep="\t")
						else:
							zero_range = str(key + 1) + "-" + str(zeroes[key])
							print(args.contig,zero_range,round(get_gc(seq[key:zeroes[key]]), 3),seq[key:zeroes[key]],sep="\t")
		exit()

	if not args.contig or not args.ref or not sortbam:
		print("Please ensure all required flags are specified; see readme file")
		exit()

	if args.verbose == True:
		print("Finding zero coverage areas in contig",args.contig)

	cmd = ["bedtools genomecov -bga -ibam %s" % sortbam]
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
				zero_print()
	zero_print()

# Scan for potential heterozygous/deletion sites - per base, discrete events, or frameshifts
## DevNote - currently skips first position
def deletion(sortbam):

	check_samtools()

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
		cmd = ["samtools depth -aa %s -r %s" % (sortbam, args.contig)]
	else:
		cmd = ["samtools depth -aa %s" % sortbam]

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

# Find deletions occurring within exons
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


# Calculate percentage coverage difference between first base in a mutation and the base before it, to determine homo/heterozygosity
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
		cmd = ["samtools depth -aa -b %s -r %s %s" % (temp_bed, args.contig, sortbam)]
	else:
		cmd = ["samtools depth -aa -b %s %s" % (temp_bed, sortbam)]

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

# Complex and contig flags: Calculate median coverage of contig, and identify regions deviating by +/- 50%; output in bed format
# Complex flag only: Calculate median coverage of each contig in assembly, and identify regions deviating by +/- 50%; output in bed format
# Simple and contig flags: Obtain a per-contig median average coverage for the specified contig
# Simple flag only: Obtain a per-contig median average coverage
## DevNote - Shows funky behaviour at stretches hovering around the threshold...
def median_deviation(sortbam):

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

	check_samtools()

	if args.contig:

		cmd = ["samtools depth -aa %s -r %s" % (sortbam, args.contig)]
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

		cmd = ["samtools depth -aa %s" % sortbam]
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

# Identify the longest continuous region of a contig where all positions fall between defined coverage limits
def coverage_limits(sortbam):

	check_samtools()

	command = ["samtools depth -aa %s -r %s" % (sortbam, args.contig)]
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

# This function extracts the sequence of the mapped reads from a part of the reference sequence specified by args.range
## DevNote - This function needs fixing
def extract_sequence(args):
	command = ("samtools mpileup -uf %s %s -r %s:%s | bcftools view -cg -") \
	% (args.ref, sortbam, args.contig, args.range)
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
