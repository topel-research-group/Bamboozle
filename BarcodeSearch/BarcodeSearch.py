#!/usr/bin/env python3

# Look for ~5,000 bp regions with conserved ends and a variable centre versus the reference
# Then check whether this region is also variable in subsequent samples
# This will allow identification of barcodes flanked by conserved primer sites

## Imports

import argparse
import subprocess

#######################################################################

## Arguments - define reference, input bams, window size, step size, and primer site size

parser = argparse.ArgumentParser(prog="BarcodeSearch")

parser.add_argument("-f", "--ref", \
			help="Reference")
parser.add_argument("-B", "--BAMs", \
			nargs="+", \
			help="BAM files of samples")
parser.add_argument("--window", \
			type=int, \
			default="5000", \
			help="Window size for barcode search")
parser.add_argument("--step", \
			type=int, \
			default="1", \
			help="Step size for barcode search")
parser.add_argument("--primer_size", \
			type=int, \
			default="21", \
			help="Desired size of conserved regions at beginning and end of barcode")
parser.add_argument("--coverage", \
			type=int, \
			default="100", \
			help="Minimum coverage required in the potential barcode region")
parser.add_argument("--minvar", \
			type=int, \
			default="1", \
			help="Minimum number of desired variant bases per barcode")

parser.add_argument("-t", "--threads", \
			default=1, \
			help="Threads")

parser.add_argument("-c", "--contig", \
			help="Specify contig to investigate")

args = parser.parse_args()

#######################################################################

print("Searching for potential barcodes in",len(args.BAMs),"file(s).\n")

## ENSURE CORRECT MODULES ARE LOADED

# Record length of contig of interest
cmd2 = ["samtools idxstats %s" % (args.BAMs[0])]
process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell=True)

with process2.stdout as result2:
	rows2 = (line2.decode() for line2 in result2)
	for row2 in rows2:
		if row2.split("\t")[0] == args.contig:
			contig_length = int(row2.split("\t")[1])

# Parse BAM files
for n in range(len(args.BAMs)):
	SNP_loci = []
	cmd = ["bcftools mpileup --threads %s --fasta-ref %s -r %s %s | bcftools call --threads %s -mv" % \
		(args.threads, args.ref, args.contig, args.BAMs[n], args.threads)]

	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

	# Generate a list of loci where variants occur
	with process.stdout as result:
		rows = (line.decode() for line in result)
		for row in rows:
			if not row.startswith("#"):
				SNP_loci.append(row.split("\t")[1])
	print(SNP_loci)

	# Start stepping through the file

	suitable_window = {}
	window_coords = []
	saved_window = (0,0,0,0)

	for window_start in range(1,(contig_length - args.window),args.step):
		window_stop = int(window_start + args.window)
		primer1_stop = int(window_start + args.primer_size)
		primer2_start = int(window_stop - args.primer_size)
		window_coords = (window_start,primer1_stop,primer2_start,window_stop)
		unsuitable = 0
		window_SNPS = 0

		# If any variants fall within primer sites, skip the window
		for SNP in SNP_loci:
			if int(SNP) in range(window_start,(primer1_stop + 1)):
				unsuitable += 1
				break
			elif int(SNP) in range(primer2_start,(window_stop + 1)):
				unsuitable += 1
				break

	# If the forward primer site overlaps with the previous entry, extend previous entry
		if unsuitable == 0:
			if window_coords[0] in range(saved_window[0],saved_window[1]):
				saved_window = (saved_window[0],window_coords[1],saved_window[2],window_coords[3])

	# Otherwise start new window
			else:
				if saved_window != (0,0,0,0):
					suitable_window[saved_window[0]] = saved_window
				saved_window = window_coords
	print(suitable_window)


## INDELS COULD PRESENT A PROBLEM



## Does the region have adequate coverage *and* adequate variation vs. reference?
## Coverage should be dealt with in BCFtools step

## Move on to next sample and repeat process, checking only coordinates saved in the first round
## Save coordinates to new list, and use this to check sample 3

## Repeat the above until all samples have been checked

## Start with a single contig, but extend to full assembly in future
