#!/usr/bin/env python3

import unittest
from os import sys, path 
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import vilmas_pipeline  
import subprocess 
import os
#import shutil


class TestClass:
	# Select specific test-files for testing functions
	def setup_args(self):
#		print(vilmas_pipeline.args)
		args = vilmas_pipeline.parser.parse_args(['-f', \
					'../example_data/reference.txt', \
					'-F', \
					'../example_data/data1_R1_00.fastq.gz', \
					'-R', \
					'../example_data/data1_R2_00.fastq.gz'])
		return args

	def setup_args2(self):
		args = vilmas_pipeline.parser.parse_args(['-f', \
					'../example_data/reference.txt', \
					'-b', \
					'../example_data/data1.bam'])
		return args
		
	def setup_args3(self):
		args = vilmas_pipeline.parser.parse_args(['-f', \
					'../example_data/reference.txt'])
		return args


class TestProcess(TestClass, unittest.TestCase):

	def setUp(self):
		current_directory = os.getcwd()
		bcftools_directory = os.path.join(current_directory, r'Bcftools')
		if not os.path.exists(bcftools_directory):
                	os.makedirs(bcftools_directory)	

	# Test the bowtie2 function in vilmas_pipeline
	def test_bowtie2(self):
		args = self.setup_args()
		self.assertIsNone(vilmas_pipeline.bowtie2(args))
#		self.assertRaises(TypeError)

	def test_reference_argument(self):
		self.args = self.setup_args()
		self.assertEqual(self.args.ref, '../example_data/reference.txt')

	def test_samtools_view(self):
		self.assertIsNone(vilmas_pipeline.samtools_view())

	def test_samtools_sort(self):
		self.assertIsNone(vilmas_pipeline.samtools_sort())

	def test_bam_input(self):
		args = self.setup_args2()
		self.assertIsNone(vilmas_pipeline.bam_input(args))

	def test_samtools_index(self):
		self.assertIsNone(vilmas_pipeline.samtools_index())

	def test_bcftools(self):
		args = self.setup_args3()
		self.assertIsNone(vilmas_pipeline.bcftools(args))

	def test_annotation(self):
		self.assertIsNone(vilmas_pipeline.annotation())

	def test_snpsift(self):
		self.assertIsNone(vilmas_pipeline.snpsift())

#	def tearDown(self):
#		# Remove an directory
#		if 'Bowtie2' in os.listdir():
#			os.rmdir('Bowtie2')
#		if 'Bcftools' in os.listdir():
#			os.rmdir('Bcftools')

if __name__ == '__main__':
        unittest.main()

