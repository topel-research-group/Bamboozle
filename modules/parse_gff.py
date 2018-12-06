#!/usr/bin/env python3

import gffutils
import os

def main(gff, feature=''):
	for file in os.listdir('.'):
		if not os.path.exists("gff.db"):
			fn = gffutils.example_filename(gff)
			db = gffutils.create_db(fn, dbfn='gff.db', force=False, keep_order=True,\
						merge_strategy='merge', sort_attribute_values=True)
			for i in db.features_of_type(feature, order_by='start'):
				return i

	db = gffutils.FeatureDB('gff.db')
	for i in db.features_of_type(feature, order_by='start'):
		return i



if __name__ == "__main__":
	main(gff, feature='')
