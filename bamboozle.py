#!/usr/bin/env python3

import sys
import subprocess
import argparse
import time

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
parser.add_argument("-e", "--events", action="store_true", help="Report number of deletion events, rather than individual positions")
parser.add_argument("-f", "--frameshift", action="store_true", help="Used with -d (equivalent to -e -f); only report mutations causing a frameshift")
parser.add_argument("-x", "--exons", help="Bed file containing exon coordinates (0-based). -m also required.")
parser.add_argument("-m", "--mutations", help="List of mutation events; output of bamboozle.py -e or -f. -x also required.")
parser.add_argument("-v", "--verbose", action="store_true", help="Be more verbose")
parser.add_argument('--dev', help=argparse.SUPPRESS, action="store_true")
args = parser.parse_args()


if args.dev == True:
	start_time = time.time()

if args.frameshift == True and args.events == False:
	args.events = True

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
			print("Obtaining stats for " + args.contig + " in " + args.bam + "; coverage >+" + str(args.threshold) + "%.")
		else:
			print("Obtaining whole-genome stats for " + args.bam + "; coverage >+" + str(args.threshold) + "%.")
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

	print("Length of assembly/contig: " + str(num_lines))
	value = 100.0 / num_lines * sum(cov_stats.values())
	print(str(round(value, 3)) + "% of the assembly/contig has >=" + str(args.threshold) + "x coverage.")



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
		print("Finding zero coverage areas in contig " + args.contig)

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
				print("GC% for contig: " + str(round(get_gc(seq), 3)))
				print("Contig\tPositions\tGC%\tSequence")
				for key in zeroes:
					if key + 1 == (zeroes[key]):
						print(args.contig + "\t" + str(zeroes[key]) + "\t-\t" + \
						str(seq[key:zeroes[key]]))
					else:
						print(args.contig + "\t" + str(key + 1) + "-" + \
						str(zeroes[key]) + "\t" + \
						str(round(get_gc(seq[key:zeroes[key]]), 3)) + "\t" + \
						str(seq[key:zeroes[key]]))



def deletion():
	# This function scans for potential heterozygous/deletion sites; it is still in development so the
	# current output is still messy

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
					if ((window[position - 11]*0.8) <= window[position] <= (window[position - 11]*1.25)) and (window[position - 11] > 0) and (window[position] >= args.threshold):
						for x, y in window.items():
							if y < (window[position - 11]*0.6) and x not in reported:
								reported.append(x)
								if args.events:
									if new_mutation(x, old_position):
										if len(deletion) != 0:
											deletion.append(del_size)
											if (args.frameshift == False) or (args.frameshift and del_size % 3 != 0):
												print(deletion[0] + "\t" + str(deletion[1]) + "\t" + str(deletion[2]))
											deletion = []
# Merge the two lines below using deletion.extend
										deletion.append(ctg)
										deletion.append(x)
										del_size = 1
									else:
										del_size +=1
								else:
									print(ctg + "\t" + str(x))
								old_position = x
				else:
					window[position] = coverage

			if args.contig and args.contig == previous_ctg:
				break

# Ensure the final event is printed

		if args.events:
			deletion.append(del_size)
			if (args.frameshift == False) or (args.frameshift and del_size % 3 != 0):
				print(deletion[0] + "\t" + str(deletion[1]) + "\t" + str(deletion[2]))

def new_mutation(new_position, old_position):
	if (int(new_position) - int(old_position)) != 1:
		return True

def exon_mutations():

	# Need to find a way to pass results of deletion function directly into this function

	frameshifts = 0
	exon_list = []
	mutation_list = []

	with open(args.exons,'r') as bed_file:
		for line in bed_file:
			exon_list.append(line)

	with open(args.mutations,'r') as in_file:
		for line in in_file:
			mutation_list.append(line)

	for mut in mutation_list:
		mut = mut.split("\t")
		for ex in exon_list:
			ex = ex.split("\t")
			if (mut[1] == ex[0]) and (int(ex[1]) <= int(mut[2]) <= int(ex[2])):
				frameshifts += 1
				print(mut[3].strip("\n") + "bp mutation at " + mut[1] + " " + mut[2] + \
				      " hits exon " + ex[3] + ".")
				break

	print("Total number of frameshifts: " + str(frameshifts))


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
	elif args.exons or args.mutations:
		if args.exons and args.mutations:
			exon_mutations()
		else:
			print("Both -x and -m must be specified to find mutations in exons")
			exit()
	elif args.zero:
		zero_regions()
	elif args.range:
		extract_sequence()
	else:
		coverage_stats()


if __name__ == "__main__":
	main()

if args.dev == True:
	print("Time taken = " + str(time.time() - start_time) + " seconds.")
