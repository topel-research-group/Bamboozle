#!/usr/bin/env python3

exon = open('time_exon_sort_flt0.2_noheader.table', 'r').read().splitlines()
gff = open('../../../../../../CDS.gff', 'r').read().splitlines() 
loop = 0
num_lines = sum(1 for line in exon)
while loop < num_lines:
	for row in exon:
		r = row.split(None)
		r = [string.strip() for string in r]
		for ro in gff:
			o = ro.split(None)
			o = [string.strip() for string in o]
			if r[0] == o[0]:
				if int(r[1]) >= int(o[3]) and int(r[1]) <= int(o[4]):
					print(o[8])
				loop += 1




