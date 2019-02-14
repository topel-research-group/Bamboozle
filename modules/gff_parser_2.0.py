#!/usr/bin/env python3

import gffutils
import os
import fnmatch
import subprocess

# Extracts feature from given gff file, output is a gff file 
# which will work as input in the annoation function 
# in the pipeline if args.gff and args.feature are used.
#######################################################
gff = '../../../Skeletonema_marinoi_Ref_v1.1_Primary.all.gff'
feature = 'CDS'
contigsizes = 'chrom.genome'
#sorted_gff =
intergenic_sorted_bed = 'intergenic_regions_sorted.bed'
exon_sorted_bed = 'exon_sorted.bed'
exon_intergenic_sorted_bed = 'tmp.exon_intergenic_sorted.txt'
intron_sorted_bed = 'introns_sorted.bed'
out = 'out.gff'

#######################################################
def main(gff, feature, contigsizes):
	if feature == 'exon':
		# Sort gff file.
#		cmd = ('''cat %s \
#			| awk '$1 ~ /^#/ {print $0;next} {print $0 \
#			| "sort -k1,1 -k4,4n -k5,5n"}' \
#			> %s''') \
#		% (gff, sorted_gff)
#		process = subprocess.Popen(cmd, \
#			stdout=subprocess.PIPE, \
#			shell=True)
#		while process.wait() is None:
#			pass
#		process.stdout.close()

		cmd2 = ('''awk '{if ($3 == "exon") print $1, $4, $5}' %s \
			| tr " " "\\t" > %s''') \
		% (gff, \
		out)
		process2 = subprocess.Popen(cmd2, \
			stdout=subprocess.PIPE, \
			shell=True)
		while process2.wait() is None:
			pass
		process2.stdout.close()

	elif feature == 'intergenic':
		# Sort gff file.
#		cmd = ('''cat %s \
#			| awk '$1 ~ /^#/ {print $0;next} {print $0 \
#			| "sort -k1,1 -k4,4n -k5,5n"}' \
#			> %s''') \
#		% (gff, sorted_gff)
#		process = subprocess.Popen(cmd, \
#			stdout=subprocess.PIPE, \
#			shell=True)
#		while process.wait() is None:
#			pass
#		process.stdout.close()
		
		# Intergenic regions
		cmd1 = ('bedtools complement -i %s -g %s > %s') \
		% (gff, \
		contigsizes, \
		out)
		process1 = subprocess.Popen(cmd1, \
			stdout=subprocess.PIPE, \
			shell=True)
		while process1.wait() is None:
			pass
		process1.stdout.close()

	elif feature == 'intron':
		# Sort gff file.
#		cmd = ('''cat %s \
#			| awk '$1 ~ /^#/ {print $0;next} {print $0 \
#			| "sort -k1,1 -k4,4n -k5,5n"}' \
#			> %s''') \
#		% (gff, sorted_gff)
#		process = subprocess.Popen(cmd, \
#			stdout=subprocess.PIPE, \
#			shell=True)
#		while process.wait() is None:
#			pass
#		process.stdout.close()

		# Intergenic regions.
		cmd1 = ('bedtools complement -i %s -g %s > %s') \
		% (gff, \
		contigsizes, \
		intergenic_sorted_bed)
		process1 = subprocess.Popen(cmd1, \
			stdout=subprocess.PIPE, \
			shell=True)
		while process1.wait() is None:
			pass
		process1.stdout.close()

		# Exons.
		cmd2 = ('''awk '{if ($3 == "exon") print $1, $4, $5}' %s \
			| tr " " "\\t" > %s''') \
		% (gff, \
		exon_sorted_bed)
		process2 = subprocess.Popen(cmd2, \
			stdout=subprocess.PIPE, \
			shell=True)
		while process2.wait() is None:
			pass
		process2.stdout.close()

		# Cat intergenic and exons.
		cmd3 = ('cat %s %s | sort -k1,1 -k2,2n > %s') \
		% (exon_sorted_bed, \
		intergenic_sorted_bed, \
		exon_intergenic_sorted_bed)
		process3 = subprocess.Popen(cmd3, \
			stdout=subprocess.PIPE, \
			shell=True)
		while process3.wait() is None:
			pass
		process3.stdout.close()
		
		# Introns.
		cmd4 = ('bedtools complement -i %s -g %s > %s') \
		% (exon_intergenic_sorted_bed, \
		contigsizes, \
		out)
		process4 = subprocess.Popen(cmd4, \
			stdout=subprocess.PIPE, \
			shell=True)
		while process4.wait() is None:
			pass
		process4.stdout.close()
		
		# Remove exon and intergenic regions files when done.
		os.remove('intergenic_regions_sorted.bed')
		os.remove('exon_sorted.bed')

	else:
		gff_path = os.path.dirname(gff)  
		file_list = os.listdir(gff_path)
		for f in file_list:
			# Open db if it exists, looks for it in the same dir as the gff file.
			if fnmatch.fnmatch(f, 'gff.db'):
				db = gffutils.FeatureDB(gff_path+'/gff.db')
				with open('out.gff', 'w') as fout:
					for i in db.features_of_type(feature, order_by='ID'):
						fout.write(str(i) + '\n')
				fout.close()
				break
		else:
			# If db does not exists a new is created.
			db = gffutils.create_db(gff, \
						dbfn=gff_path+'/gff.db', \
						force=False, \
						keep_order=True, \
						merge_strategy='merge', \
						sort_attribute_values=True)

			db = gffutils.FeatureDB(gff_path+'/gff.db')
			with open('out.gff', 'w') as fout:
				for i in db.features_of_type(feature, order_by='ID'):
					fout.write(str(i) + '\n')
			fout.close()


if __name__ == "__main__":
	main(gff, feature, contigsizes)

