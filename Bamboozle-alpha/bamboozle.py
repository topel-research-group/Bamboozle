#!/usr/bin/env python3

import sys
import subprocess
import argparse
import time
import os.path
from statistics import median

parser = argparse.ArgumentParser(description='Obtain statistics from BAM files')

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
parser.add_argument('-t', '--threads', type=int, nargs='?', default='1', help='Number of threads; default 1')
parser.add_argument('-d', '--threshold', type=int, nargs='?', const='1', default='20', help='Threshold for calculating coverage percentage; default 20')
parser.add_argument("-f", "--ref", help="Reference sequence file")
parser.add_argument("--sortbam", help="Sorted bam file")
parser.add_argument("-a", "--range", help="somethingsomsing")
parser.add_argument("-m", "--mutations", help="List of mutation events; currently requires output from bamboozle deletion function")
parser.add_argument("-x", "--exons", help="Bed file containing exon coordinates (0-based); -m also required")
parser.add_argument("-l", "--limits", type=int, nargs=2, help="Specify lower and upper limits for long_coverage function; two arguments required")
parser.add_argument("-v", "--verbose", action="store_true", help="Be more verbose")
parser.add_argument('--dev', help=argparse.SUPPRESS, action="store_true")
args = parser.parse_args()


if args.dev == True:
	start_time = time.time()


# Calculate percentage of positions in assembly/contig with read coverage >= a given threshold (default: 20x)
def coverage_stats():

	check_samtools()

	if args.contig:
		cmd = ["samtools depth -aa %s -r %s" % (args.sortbam, args.contig)]
		sequence = "contig"
	else:
		cmd = ["samtools depth -aa %s" % args.sortbam]
		sequence = "assembly"

	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

	if args.verbose == True:
		if args.contig:
			print("Obtaining stats for ",args.contig," in ",os.path.basename(args.sortbam),"; coverage >+",args.threshold,"x.",sep="")
		else:
			print("Obtaining whole-genome stats for ",os.path.basename(args.sortbam),"; coverage >+",args.threshold,"x.",sep="")
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
def zero_regions():

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

	if not args.contig or not args.ref or not args.sortbam:
		print("Please ensure all required flags are specified; see readme file")
		exit()

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
				zero_print()
	zero_print()

# Scan for potential heterozygous/deletion sites - per base, discrete events, or frameshifts
## DevNote - currently skips first position
def deletion():

	check_samtools()

	if args.verbose == True:
		if not args.contig:
			args.contig = "assembly"
		if args.deletion1:
			print("Finding individual deletions in",args.contig)
		elif args.deletion2:
			print("Finding deletion events in",args.contig)
		elif args.deletion3:
			print("Finding frameshift deletion events in",args.contig)

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
							if args.deletion2 or args.deletion3:
								if (int(x) - int(old_position)) != 1:
									if len(deletion) != 0:
										print_deletion(deletion, del_size)
										deletion = []
									deletion.extend([ctg,x])
									del_size = 1
								else:
									del_size +=1
							else:
								print(ctg,x,sep="\t")
							old_position = x
			else:
				window[position] = coverage

		if args.deletion2 or args.deletion3:	# Ensure that the final event is reported
			print_deletion(deletion, del_size)


def print_deletion(m, n):
	m.append(n)
	if not args.deletion3 or (args.deletion3 and n % 3 != 0):
		print(m[0],m[1],m[2],sep="\t")

# Find deletions occurring within exons
## DevNote - Need to find a way to pass results of deletion function directly into this function
def exon_mutations():

	frameshifts = 0
	exon_list = []
	mutation_list = []

	list_append(args.exons, exon_list)
	list_append(args.mutations, mutation_list)

	for mut in mutation_list:
		mut = mut.split("\t")
		for ex in exon_list:
			ex = ex.split("\t")
			if (mut[0] == ex[0]) and (int(ex[1]) <= int(mut[1]) <= int(ex[2])):
				if int(mut[2]) % 3 != 0:
					frameshifts += 1
				if args.verbose:
					print(mut[2].strip("\n"),"bp mutation at",mut[0],mut[1],"hits exon",ex[3])
				else:
					print(mut[0],mut[1],mut[2].strip("\n"),ex[3],sep="\t")
				break

	print("Total number of frameshifts in exons:",frameshifts)


def list_append(argument, list):
	with open(argument, 'r') as input:
		for line in input:
			list.append(line)


# Calculate percentage coverage difference between first base in a mutation and the base before it, to determine homo/heterozygosity
def HomoDel_or_Hetero():

	check_samtools()

	# Read in each line of args.mutations, then print out an altered version to a temporary .bed file

	temp_bed = "TEMP_" + os.path.basename(args.mutations) 

	if os.path.isfile(temp_bed) == True:
		print("Temporary file can't be written; please ensure",temp_bed,"is not a file.")
		sys.exit()

	rows = (line.split('\t') for line in open(args.mutations))
	for row in rows:
		with open(temp_bed, "a") as output_file:
			output_file.write(row[0] + "\t" + str(int(row[1]) - 2) + "\t" + row[1] + "\n")
	output_file.close()

	# Compare the coverage 

	cmd = ["samtools depth -aa -b %s %s" % (temp_bed, args.sortbam)]
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

	coverage = []

	with process.stdout as result:
		rows = (line.decode().split('\t') for line in result)
		for row in rows:
			coverage.append(row)

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
def median_deviation():

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

# Identify the longest continuous region of a contig where all positions fall between defined coverage limits
def coverage_limits():

	check_samtools()

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

# This function extracts the sequence of the mapped reads from a part of the reference sequence specified by args.range
## DevNote - This function needs fixing
def extract_sequence():
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


# Ensure correct version of samtools
def check_samtools():
	try:
		subprocess.check_output('samtools depth 2>&1 | grep -- "-aa"', stderr=subprocess.PIPE, shell=True)
	except subprocess.CalledProcessError:
		print("This version of samtools does not support the `depth -aa` option; please update samtools.")
		exit()


def main():
	if args.deletion1 or args.deletion2 or args.deletion3:
		deletion()
	elif args.deletionx:
		if args.mutations:
			exon_mutations()
		else:
			print("Both -x and -m must be specified to find mutations in exons")
			exit()
	elif args.homohetero:
		HomoDel_or_Hetero()
	elif args.zero:
		zero_regions()
	elif args.consensus:
		extract_sequence()
	elif args.median:
		if args.simple or args.complex:
			median_deviation()
		else:
			print("Please specify --simple for medians only or --complex for full output")
			exit()
	elif args.coverage:
		coverage_stats()
	elif args.long_coverage:
		coverage_limits()
	else:
		parser.print_help(sys.stderr)
		exit()

if __name__ == "__main__":
	main()

if args.dev == True:
	print("Time taken =",(time.time() - start_time),"seconds.")
