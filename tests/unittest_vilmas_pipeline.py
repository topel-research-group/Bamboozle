#!/usr/bin/env python3

import unittest
from os import sys, path 
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import vilmas_pipeline  
import subprocess 
import os
import shutil


class TestClass:
	# Select specific test-files for testing functions
	def setup_args(self):
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
	
	# Make directories
	def setUp(self):
		current_directory = os.getcwd()
		bowtie2_directory = os.path.join(current_directory, r'Bowtie2')
		if not os.path.exists(bowtie2_directory):
			os.makedirs(bowtie2_directory)

		bcftools_directory = os.path.join(current_directory, r'Bcftools')
		if not os.path.exists(bcftools_directory):
                	os.makedirs(bcftools_directory)	

	# Test that functions output is None in pipeline 
	def test010_bowtie2(self):
		args = self.setup_args()
		self.assertIsNone(vilmas_pipeline.bowtie2(args))

	def test020_reference_argument(self):
		self.args = self.setup_args()
		self.assertEqual(self.args.ref, '../example_data/reference.txt')

	def test030_samtools_view(self):
		self.assertIsNone(vilmas_pipeline.samtools_view())

	def test040_samtools_sort(self):
		self.assertIsNone(vilmas_pipeline.samtools_sort())

	def test050_bam_input(self):
		args = self.setup_args2()
		self.assertIsNone(vilmas_pipeline.bam_input(args))

	def test060_samtools_index(self):
		self.assertIsNone(vilmas_pipeline.samtools_index())

	def test070_bcftools(self):
		args = self.setup_args3()
		self.assertIsNone(vilmas_pipeline.bcftools(args))
		
	def test080_annotation(self):
		self.assertIsNone(vilmas_pipeline.annotation())

	def test090_snpsift(self):
		self.assertIsNone(vilmas_pipeline.snpsift())

	# Test output files
	def test100_sam_output(self):
		test_sam = open('Bowtie2/tests.sam')
		expected_sam = open('../example_data/data1.sam')
		self.assertEqual(test_sam.readlines(), expected_sam.readlines())

	def test110_output_bam(self):
		myProcess = subprocess.check_call('diff Bowtie2/tests.bam ../example_data/data1.bam', shell=True)
		self.assertIs(myProcess, 0)
		
	def test120_snpsift_output(self):
		test_snpeff = open('Bcftools/tests.snpsift_table.txt', 'r')
		expected_snpeff = open('../example_data/data1.snpsift_table.txt', 'r')
		self.assertEqual(test_snpeff.readlines(), expected_snpeff.readlines())


	def tearDown(self):
		pass

	@classmethod
	def tearDownClass(cls):
		# Remove the directories after the test is done
		if 'Bowtie2' in os.listdir():
			shutil.rmtree('Bowtie2')
		if 'Bcftools' in os.listdir():
			shutil.rmtree('Bcftools')
		

if __name__ == '__main__':
        unittest.main()

