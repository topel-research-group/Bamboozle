#!/usr/bin/env python3

import unittest
from os import sys, path 
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import vilmas_pipeline  
import subprocess 
import os

class TestClass:
	# Select specific test-files for testing functions
	def setup_args(self):
#		print(vilmas_pipeline.args)
		args = vilmas_pipeline.parser.parse_args(['-f', \
					'../example_data/reference.txt', \
					'-F', \
					'../data1_R1.fastq.gz', \
					'-R', \
					'../data1_R2.fastq.gz'])
		return args
		
class TestProcess(TestClass, unittest.TestCase):

	def setUp(self):
		pass

	# Test the bowtie2 function in vilmas_pipeline
	def test_bowtie2(self):
		args = self.setup_args()
		self.assertIsNone(vilmas_pipeline.bowtie2(args))
#		self.assertRaises(TypeError)

	def test_reference_argument(self):
		self.args = self.setup_args()
		self.assertEqual(self.args.ref, '../example_data/test_Skeletonema_marinoi_Ref.txt')

	def tearDown(self):
		# Remove an empty directory
		if 'Bowtie2' in os.listdir():
			os.rmdir('Bowtie2')


if __name__ == '__main__':
        unittest.main()

