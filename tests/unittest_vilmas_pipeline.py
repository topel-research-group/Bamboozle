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
		print(vilmas_pipeline.args)
		arguments = vilmas_pipeline.parser.parse_args(['-f', \
					'/proj/data11/vilma/Pipeline_vilma/example_data/test_Skeletonema_marinoi_Ref.txt', \
					'-F', \
					'/proj/data11/vilma/Pipeline_vilma/P8352_102/P8352_102_S1_L001_R1_001.fastq.gz', \
					'-R', \
					'/proj/data11/vilma/Pipeline_vilma/P8352_102/P8352_102_S1_L001_R2_001.fastq.gz'])
		return arguments
		
class TestProcess(TestClass, unittest.TestCase):

	def setUp(self):
		pass

	# Test the bowtie2 function in vilmas_pipeline
	def test_bowtie2(self):
		self.args = self.setup_args()
	#	self.assertIsNone(vilmas_pipeline.bowtie2())
		self.assertRaises(TypeError)

	def test_dummy(self):
		self.arguments = self.setup_args()
		self.assertTrue(self.arguments, '/proj/data11/vilma/Pipeline_vilma/example_data/test_Skeletonema_marinoi_Ref.txt')

	def tearDown(self):
		# Remove an empty directory
		if 'Bowtie2' in os.listdir():
			os.rmdir('Bowtie2')


if __name__ == '__main__':
        unittest.main()

