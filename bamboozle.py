#!/usr/bin/env python3

import sys
import subprocess
import argparse
import time
import os.path
from statistics import mode
from collections import Counter

parser = argparse.ArgumentParser(description='Obtain statistics regarding percentage coverage from bam files. \
                                              The script gives percentage of positions in an assembly/contig \
                                              with coverage greater than or equal to a given threshold')
parser.add_argument('-c', '--contig', help='Gives per-contig coverage stats')
#parser.add_argument('-c', '--contig', type=str, nargs='+', help='Gives cov. stats for the specified contigs')
parser.add_argument('-t' ,'--threshold', type=int, nargs='?', const=1, default='20', help='Threshold for calculating coverage percentage; default 20')
parser.add_argument("-r", "--refference", help="Reference sequence file")
parser.add_argument("-b", "--bam", help="Bam file")
parser.add_argument("--range", help="somethingsomsing")
parser.add_argument("-z", "--zero", action="store_true", help="Find regions of 0x coverage")
parser.add_argument("-d", "--deletion", action="store_true", help="Scan for potential deletions")
parser.add_argument("-e", "--events", action="store_true", help="Report deletion events, rather than individual positions")
parser.add_argument("-f", "--frameshift", action="store_true", help="Report frameshift deletions, rather than individual positions")
parser.add_argument("-m", "--mutations", help="List of mutation events; output of bamboozle.py -d -e/-f")
parser.add_argument("-x", "--exons", help="Bed file containing exon coordinates (0-based). -m also required.")
parser.add_argument("-o", "--homohetero", action="store_true", help="Determine whether a given deletion is homo- or heterozygous; WIP")
parser.add_argument("--mode", action="store_true", help="Report regions whose coverage differs by +/- >50% of the contig mode; single contig.")
parser.add_argument("--modeall", action="store_true", help="Report regions whose coverage differs by +/- >50% of the contig mode; whole assembly.")
parser.add_argument("-v", "--verbose", action="store_true", help="Be more verbose")
parser.add_argument('--dev', help=argparse.SUPPRESS, action="store_true")
args = parser.parse_args()


if args.dev == True:
	start_time = time.time()

if args.frameshift == True:
	args.events = True

if args.events == True:
	args.deletion = True


def coverage_stats():
	# This function calculates the percentage of positions in an assembly/contig
	# with read coverage >= a given threshold (default: 20x)

	try:
		subprocess.check_output('samtools depth 2>&1 | grep -- "-aa"', stderr=subprocess.PIPE, shell=True)
	except subprocess.CalledProcessError:
		print("This version of samtools does not support the `depth -aa` option; please update samtools.")
		exit()

	cmd = ["samtools depth -aa %s" % args.bam]
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

	if args.verbose == True:
		if args.contig:
			print("Obtaining stats for ",args.contig," in ",os.path.basename(args.bam),"; coverage >+",args.threshold,"%.",sep="")
		else:
			print("Obtaining whole-genome stats for ",os.path.basename(args.bam),"; coverage >+",args.threshold,"%.",sep="")
	cov_stats = {}
	num_lines = 0
	with process.stdout as result:
		rows = (line.decode().split('\t') for line in result)
		for row in rows:
			coverage = int(row[2])
			if args.contig:
				if str(row[0]) == args.contig:
					num_lines += 1
					if coverage >= args.threshold:
						if coverage in cov_stats:
							cov_stats[coverage] += 1
						else:
							cov_stats[coverage] = 1
			else:
				num_lines += 1
				if coverage >= args.threshold:
					if coverage in cov_stats:
						cov_stats[coverage] += 1
					else:
						cov_stats[coverage] = 1

	print("Length of assembly/contig:",num_lines)
	value = 100.0 / num_lines * sum(cov_stats.values())
	print(round(value, 3),"% of the assembly/contig has >=",args.threshold,"x coverage.",sep="")



def zero_regions():
	# This function identifies regions of 0x coverage in a given contig, then prints
	# the reference sequence and GC content at these coordinates

	# Devel. Add try except statement for bedtools here
        # Also ensure that bam is sorted?

	# Fasta parsing from Biopython, fp.py and
	# http://stackoverflow.com/questions/7654971/parsing-a-fasta-file-using-a-generator-python

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

	if not args.contig:
		print("Please specify contig with the -c flag")
		exit()

	if not args.refference:
		print("Please specify reference with the -r flag")
		exit()

	if args.verbose == True:
		print("Finding zero coverage areas in contig",args.contig)

	cmd = ["bedtools genomecov -bga -ibam %s" % args.bam]
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
	zeroes = {}

	with process.stdout as result:
		rows = (line.decode().split('\t') for line in result)
		for row in rows:
			if str(row[0]) == args.contig:
				coverage = int(row[3])
				if coverage == 0:
					zeroes[int(row[1])] = int(row[2])

	with open(args.refference) as fasta:
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



def deletion():
	# This function scans for potential heterozygous/deletion sites, either per-base or as discrete events.
        # Non-frameshift deletions can also be filtered out

	try:
		subprocess.check_output('samtools depth 2>&1 | grep -- "-aa"', stderr=subprocess.PIPE, shell=True)
	except subprocess.CalledProcessError:
		print("This version of samtools does not support the `depth -aa` option; please update samtools.")
		exit()

	cmd = ["samtools depth -aa %s" % args.bam]
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

	# Currently skips the first position...

			if (args.contig == None) or (args.contig and args.contig == str(row[0])):
				if len(window) == 12:
					del window[position - 12]
					window[position] = coverage
					base1 = window[position - 11]
					if ((base1*0.8) <= window[position] <= (base1*1.25)) and (base1 > 0) and (window[position] >= args.threshold):
						for x, y in window.items():
							if y < (base1*0.6) and x not in reported:
								reported.append(x)
								if args.events:
									if new_mutation(x, old_position):
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

			if args.contig and args.contig == previous_ctg:
				break

		if args.events:					# Ensure that the final event is reported
			print_deletion(deletion, del_size)


def new_mutation(new_position, old_position):
	if (int(new_position) - int(old_position)) != 1:
		return True

def print_deletion(m, n):
	m.append(n)
	if (args.frameshift == False) or (args.frameshift and n % 3 != 0):
		print(m[0],m[1],m[2],sep="\t")



def exon_mutations():
	# Need to find a way to pass results of deletion function directly into this function

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



def HomoDel_or_Hetero():
	# This function calculates the percentage coverage difference between the first base in a mutation and the
	# base before it; this allows the user to determine whether a given deletion is homozygous or heterozygous

	try:
		subprocess.check_output('samtools depth 2>&1 | grep -- "-aa"', stderr=subprocess.PIPE, shell=True)
	except subprocess.CalledProcessError:
		print("This version of samtools does not support the `depth -aa` option; please update samtools.")
		exit()

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

	cmd = ["samtools depth -aa -b %s %s" % (temp_bed, args.bam)]
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


def mode_deviation():
	try:
		subprocess.check_output('samtools depth 2>&1 | grep -- "-aa"', stderr=subprocess.PIPE, shell=True)
	except subprocess.CalledProcessError:
		print("This version of samtools does not support the `depth -aa` option; please update samtools.")
		exit()

	# Calculate mode coverage of contig, then identify regions which deviate from this by +/- 50%
	# Output in bed format

	cmd = ["samtools depth -aa %s" % args.bam]
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

	cov_stats = {}
	check_me = 0
	current_contig = args.contig
	with process.stdout as result:
		rows = (line.decode().split('\t') for line in result)
		for row in rows:
			ctg = str(row[0])
			position = int(row[1])
			coverage = int(row[2])
			if ctg == args.contig:
				check_me = 1
				cov_stats[position] = coverage
			elif check_me == 1:
				break
	try:
		mode_cov = mode(cov_stats.values())
	except:
		print("# Note: There is not a single mode in this contig; please check results")
		value_counts = Counter(cov_stats.values())
		for item, frequency in value_counts.most_common(1):
			mode_cov = item

	low_threshold = mode_cov * 0.5	# This is open to change as needed
	high_threshold = mode_cov * 1.5
#	print("Mode of coverage =",mode_cov)
#	print("Low threshold =",low_threshold)
#	print("High threshold =",high_threshold)


	# Shows funky behaviour at stretches hovering around the threshold...

	print("track name=WeirdCoverage","description='Areas +/- 50% of the mode coverage'",sep="\t")

	make_bed(cov_stats,current_contig,low_threshold,high_threshold)


def mode_deviation_all():
	# Shows funky behaviour at stretches hovering around the threshold...

	try:
		subprocess.check_output('samtools depth 2>&1 | grep -- "-aa"', stderr=subprocess.PIPE, shell=True)
	except subprocess.CalledProcessError:
		print("This version of samtools does not support the `depth -aa` option; please update samtools.")
		exit()

	cmd = ["samtools depth -aa %s" % args.bam]
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

	cov_stats = {}
	problem_contigs = []
	current_contig = "None"
	print("track name=WeirdCoverage","description='Areas +/- 50% of the mode coverage'",sep="\t")

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
				try:
					mode_cov = mode(cov_stats.values())
				except:
					problem_contigs.append(current_contig)
					value_counts = Counter(cov_stats.values())
					for item, frequency in value_counts.most_common(1):
						mode_cov = item

				low_threshold = mode_cov * 0.5
				high_threshold = mode_cov * 1.5

				make_bed(cov_stats,current_contig,low_threshold,high_threshold)

				cov_stats = {}
				current_contig = ctg
				cov_stats[position] = coverage
	print("# Note: the following contigs do not have a single mode; please review them")
	for item in problem_contigs:
		print("#",item)

def make_bed(contig_lib,this_contig,lower,upper):
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


def extract_sequence():
	# This function extracts the sequence of the mapped reads 
	# from a part of the reference sequence specified by args.range
	command = ("samtools mpileup -uf %s %s -r %s:%s | bcftools view -cg -") \
	% (args.refference, args.bam, args.contig, args.range)
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



def main():
	if args.deletion:
		deletion()
	elif args.exons:
		if args.mutations:
			exon_mutations()
		else:
			print("Both -x and -m must be specified to find mutations in exons")
			exit()
	elif args.homohetero:
		HomoDel_or_Hetero()
	elif args.zero:
		zero_regions()
	elif args.range:
		extract_sequence()
	elif args.mode:
		mode_deviation()
	elif args.modeall:
		mode_deviation_all()
	else:
		coverage_stats()


if __name__ == "__main__":
	main()

if args.dev == True:
	print("Time taken =",(time.time() - start_time),"seconds.")
