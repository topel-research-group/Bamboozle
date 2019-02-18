#!/usr/bin/env python3


import unittest
from os import sys, path
sys.path.append('../scripts')
import eq_system 
import subprocess
import os
import fnmatch 

# Tests the equation system made to use the fst values 
# from angsd_fst.py to calculate degree of environmental 
# impact due to warming.

class TestProcess(unittest.TestCase):

	def setUp(self):
		# Run script.
		cmd = ('../scripts/eq_system.py', \
			'-1', '../example_data/ex_1.csv', \
			'-2', '../example_data/ex_2.csv', \
			'-3', '../example_data/ex_3.csv')
		process = subprocess.Popen(cmd, \
			stdout=subprocess.PIPE)	
		while process.wait() is None:
			pass
		process.stdout.close()

	# Test output files.
	def test001_output_temp(self):
		myProcess = subprocess.check_call('diff \
			temp.csv \
			../example_data/temp.csv', \
			shell=True)
		self.assertIs(myProcess, 0)

	def test002_output_loc(self):
		myProcess = subprocess.check_call('diff \
			local.csv \
			../example_data/local.csv', \
			shell=True)
		self.assertIs(myProcess, 0)

	def test003_output_time(self):
		myProcess = subprocess.check_call('diff \
			time.csv \
			../example_data/time.csv', \
			shell=True)
		self.assertIs(myProcess, 0)

	def tearDown(self):
		pass

	@classmethod
	# Remove output files.
	def tearDownClass(cls):
		for file in os.listdir():
			if fnmatch.fnmatch(file, '*.csv'): 
				os.remove(file)	
			if fnmatch.fnmatch(file, '*.png'):
				os.remove(file)

if __name__ == '__main__':
	unittest.main()
