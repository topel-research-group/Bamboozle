#!/usr/bin/env python3


# Usage: unittest_pipeline.py -f <any key> 


import unittest
from os import sys, path 
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import bamboozle
from modules import pipeline 
import subprocess 
import os
import shutil


class TestClass:
    # Select specific test files for testing functions.
#   def setup_args(self):
#       args = bamboozle.parser.parse_args(['-f', \
#                   '../example_data/reference.txt', \
#                   '-F', \
#                   '../example_data/data1_R1_00.fastq.gz', \
#                   '-R', \
#                   '../example_data/data1_R2_00.fastq.gz'])
#       return args

    def setup_args2(self):
        args = (bamboozle.input_commands.parse_args(['-f', \
                    '../example_data/reference.txt', \
                    '-b', \
                    '../example_data/data1.bam']))
        args = (bamboozle.other_commands.parse_args(['-t', "1"]))
        print(args)
#        return args

    def setup_args3(self):
        args = bamboozle.input_commands.parse_args(['-f', \
                    '../example_data/reference.txt'])
        return args

#   def setup_args4(self):
#       args = bamboozle.input_commands.parse_args(['-f', \
#                   '../example_data/reference.txt', \
#                   '-b', \
#                   '../example_data/data1.bam', \
#                   '--gff', \
#                   '../example_data/data1.gff', \
#                   '--feature', \
#                   'exon'])
#       return args


class TestProcess(TestClass, unittest.TestCase):
    
    def setUp(self):
        pass

    # Test reference input.
    def test070_bcftools(self):
        args = self.setup_args2()
        self.assertIsNone(pipeline.bcftools(args))

    def test080_annotation(self):
        args = self.setup_args2()
        self.assertIsNone(pipeline.annotation(args))

#   def test081_annotation(self):
#       args = self.setup_args4()
#       self.assertIsNone(pipeline.annotation(args))

    def test090_snpsift(self):
        args = self.setup_args2()
        self.assertIsNone(pipeline.snpsift(args))

    # Test output files.
    def test110_snpsift_output(self):
        test_snpsift = open('Bcftools/tests.snpsift_table.txt', 'r')
        expected_snpsift = open('../example_data/data1.snpsift_table.txt', 'r')
        self.assertEqual(test_snpsift.readlines(), expected_snpsift.readlines())
        test_snpsift.close()
        expected_snpsift.close()

#   def test120_gff_parser_output(self):
#       myProcess = subprocess.check_call('diff \
#           out.gff \
#           ../example_data/out.gff', \
#           shell=True)
#       self.assertIs(myProcess, 0)

    def tearDown(self):
        pass

    @classmethod
    def tearDownClass(cls):
    # Remove the directories and files after the test is done.
            pass


if __name__ == '__main__':
    options = [sys.argv[0]] + [a for a in sys.argv if a.startswith("-")]
    unittest.main(argv=options)
