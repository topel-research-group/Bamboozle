#!/usr/bin/env python

import sys
import subprocess
import argparse

# This version of the script was designed to look at per-contig per-base coverage; this is being altered to
# per-base coverage of the entire genome in v3, with per-contig per-base coverage being added as an extension
# later in development

###################################################################################

parser = argparse.ArgumentParser(description='This script is a work in progress.')
#parser.add_argument('InputBam', metavar='input', type=str, help='Bam file to take stats from.')
#parser.add_argument('Threshold', metavar='threshold', type=int, help='Threshold for calculating coverage percentage.')

args = parser.parse_args()

###################################################################################

samtools = str('/usr/local/packages/samtools-1.3.1/samtools')
#input_bam = sys.argv[1]
input_bam = str('/proj/data17/Skeletonema_marinoi_adaptation_to_warming_project/01_mapping/P8352_101/P8352_101_sorted.bam')
#threshold = sys.argv[2]
threshold = 20

# samtools idxstats gives a tab-delimited output of [contig, seq length, mapped reads, unmapped reads] to STDOUT
# First two columns are saved to dictionary contig_stats

contig_stats = {}

cmd1 = [samtools, 'idxstats', input_bam]
process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
with process1.stdout as result1:
	rows1 = (line1.split('\t') for line1 in result1)
	contig_stats = {row1[0]:row1[1] for row1 in rows1}
del contig_stats['*']

#print contig_stats
#for key in contig_stats:
#	print key


# Generate a per-base coverage for each contig (samtools can apparently achieve this, so use this rather than bedtools [reduces no. of dependencies by one])
# Output = Contig-Position-Coverage, tab-separated
# Save output to new dictionary coverage_stats

coverage_stats = {}

cmd2 = [samtools, 'depth', '-aa', input_bam]
process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE)

# For key in contig_stats, how many instances of that EXACT key have a value of 'threshold' or above in column 3 ([2])?

with process2.stdout as result2:
	rows2 = (line2.split('\t') for line2 in result2)
	for key in contig_stats:
		coverage_stats[key] = 0
	for row2 in rows2:
		if (int(row2[2]) >= threshold):
			print row2
			coverage_stats[str(row2[0])] += 1
	for key in contig_stats:
		print key + ' : ' + str(coverage_stats[key])

#		print key
#		for row2 in rows2:
#			if (str(row2[0]) == key) and (row2[2] >= threshold):
#			if (str(row2[0]) == key):
#				print key + " has been found!"
#				coverage_stats[key] += 1
#		print key + ' : ' + str(coverage_stats[key])


# Grep each contig name, cut to extract the coverage line, then
# sum the number of instances of coverage below the threshold

#IFS=$'\n'; for z in $(cat ${INPUT1} | cut -d' ' -f2);
#do X=$(grep -P "^${z}\t" ${INPUT2} | cut -f3 | awk -v a="$THRESHOLD" '$1>=a{c++} END{print c+0}');
#N=$(grep ${z} ${INPUT1} | cut -d' ' -f1);
#Y=$(echo - | awk -v X="$X" -v N="$N" '{print 100 / N * X}');
#echo "Percentage of ${z} with >=${THRESHOLD}x coverage = ${Y}" >> output.txt;
#done
