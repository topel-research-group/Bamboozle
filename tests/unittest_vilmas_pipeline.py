#!/usr/bin/env python3

import unittest
from os import sys, path 
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import vilmas_pipeline  
import subprocess 
import argparse 

class TestClass:
	# Select specific test-files for testing functions
	def setup_args(self):
		print(vilmas_pipeline.args)
		arguments = vilmas_pipeline.parser.parse_args(['-f', '/proj/data11/vilma/Pipeline_vilma/example_data/test_Skeletonema_marinoi_Ref.txt', \
		'-F', '/proj/data11/vilma/Pipeline_vilma/P8352_102/P8352_102_S1_L001_R1_001.fastq.gz', '-R', '/proj/data11/vilma/Pipeline_vilma/P8352_102/P8352_102_S1_L001_R2_001.fastq.gz'])
		return arguments
		
class TestProcess(TestClass, unittest.TestCase):

	def setUp(self):
		pass

	# Should test the bowtie2 function in vilmas_pipeline
	def test_bowtie2(self):
		self.args = self.setup_args()
	#	print (vilmas_pipeline.bowtie2())
	#	self.assertIsNone(vilmas_pipeline.bowtie2())
		self.assertRaises(TypeError)



if __name__ == '__main__':
        unittest.main()

