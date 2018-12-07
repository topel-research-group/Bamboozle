#!/usr/bin/env python3

import gffutils
import os
import fnmatch

# Extracts feature from given gff file, output is a gff file that is input in annoation function 
# (snpEff -interval out.gff ...) in the pipeline if args.gff and args.feature
def main(gff, feature=''):
	gff_path = os.path.dirname(gff)  
	file_list = os.listdir(gff_path)
	for f in file_list:
		if fnmatch.fnmatch(f, 'gff.db'):
			db = gffutils.FeatureDB(gff_path+'/gff.db')
			with open('out.gff', 'w') as fout:
				for i in db.features_of_type(feature, order_by='ID'):
					fout.write(str(i) + '\n')
			fout.close()
			break
	else:
		db = gffutils.create_db(gff, dbfn=gff_path+'/gff.db', force=False, keep_order=True,\
					merge_strategy='merge', sort_attribute_values=True)

		db = gffutils.FeatureDB(gff_path+'/gff.db')
		with open('out.gff', 'w') as fout:
			for i in db.features_of_type(feature, order_by='ID'):
				fout.write(str(i) + '\n')
		fout.close()


if __name__ == "__main__":
	main(gff, feature='')
