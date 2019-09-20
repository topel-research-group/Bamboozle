#!/usr/bin/env python3

# With a gff and the output of coverage_stats, identify which genes get mapped to
# Intended for use with the S. subsalsum data
#
# Will require altering coverage_stats

#######################################################################
# IMPORTS
#######################################################################

import argparse
import subprocess
import re

#######################################################################
# ARGUMENTS
#######################################################################

parser = argparse.ArgumentParser(prog="Bamboozle")

parser.add_argument("--sortbam", \
			nargs='*', \
			help="Sorted BAM infile (N.B. only the BarcodeSearch function accepts multiple inputs)")
parser.add_argument('-c', '--contig', \
			help='Gives per-contig coverage stats')
parser.add_argument("--gff", \
			help="gff infile")
barcode.add_argument("-o", "--outfile", \
			help="Output filename")

parser.add_argument('--dev', \
			help=argparse.SUPPRESS, action="store_true")

args = parser.parse_args()

#######################################################################
# TESTING RUNTIME (--dev)
#######################################################################

if args.dev == True:
	start_time = time()

#######################################################################

def coverage_stats_2():

# Generate a dictionary of dictionaries of lists; contig -> gene name -> coordinates

	if args.gff:
		in_gff_loci = {}

		with open(args.gff, 'r') as input:
			lines = (entry.split("\t") for entry in input)
			for line in lines:

				# Is this standard in gff files?
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

#######################################################################

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

	with open(args.outfile, "a") as output_file:
		output_file.write("track name=Coverage description=\"Coverage above " + args.threshold + "x\"\n")



	cov_stats = {}
	num_lines = 0

# USE CODE FROM BAMPARSER'S COVERAGE_LIMITS AS GUIDELINE FOR SAVING LOCI

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

		output_file.close()



help="Output filename")

#######################################################################

if args.dev == True:
	print("Time taken =",(time() - start_time),"seconds.")
