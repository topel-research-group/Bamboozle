#!/usr/bin/env python3

import unittest
from os import sys, path 
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import vilmas_pipeline  
import subprocess 
import os
import shutil
import gzip 
import pickle

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
		bowtie2_directory = os.path.join(current_directory, r'Bowtie2')
		if not os.path.exists(bowtie2_directory):
			os.makedirs(bowtie2_directory)

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

	# Test output
	def test_sam_output(self):
		test_sam = open('Bowtie2/tests.sam')
		expected_sam = open('../example_data/data1.sam')
		self.assertEqual(test_sam.readlines(), expected_sam.readlines())

	def setup_bam_output(self):
		test_bam = ('Bowtie2/tests_bam_md5sum.txt')
		expected_bam = ('../example_data/data1_bam_md5sum.txt')
		myProcess = subprocess.Popen('md5sum -b Bowtie2/tests.bam > Bowtie2/tests_bam_md5sum.txt', stdout=subprocess.PIPE, shell=True)

	def test_bam_output(self):
#		myProcess2 = subprocess.Popen('md5sum -c Bowtie2/tests_bam_md5sum.txt ../example_data/data1_bam_md5sum.txt', stdout=subprocess.PIPE, shell=True)
		myProcess2 = subprocess.check_call('diff Bowtie2/tests.bam ../example_data/data1.bam', shell=True)
		self.assertIs(myProcess2, 0)
		

#	def test_bcftools_output(self):
#		test_bcftools = ('Bcftools/tests.bcftools_filtered.vcf.gz')
#		expected_bcftools = ('../example_data/data1.bcftools.vcf.gz')
#		myProcess = subprocess.run('zcat ../example_data/data1.bcftools.vcf.gz', shell=True, stdout=subprocess.PIPE)
#		myProcess2 = subprocess.run('zcat Bcftools/tests.bcftools_filtered.vcf.gz', shell=True, stdout=subprocess.PIPE)
#		self.assertEqual(myProcess.stdout.readlines(), myProcess2.stdout.readlines())


#	def tearDown(self):
#		# Remove an directory
#		if 'Bowtie2' in os.listdir():
#			shutil.rmtree('Bowtie2')
#		if 'Bcftools' in os.listdir():
#			shutil.rmtree('Bcftools')

if __name__ == '__main__':
        unittest.main()

