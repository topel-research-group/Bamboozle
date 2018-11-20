
import unittest
from os import sys, path 
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import vilmas_pipeline  
import subprocess #import Popen, PIPE
import argparse 



def parse_args(args):	
	parser = argparse.ArgumentParser(prog="ADD-SCRIPT-NAME-HERE")
	parser.add_argument("-f", "--ref")#, required=True, help="Reference")
	parser.add_argument("-F", "--forward", nargs='*', help="Forward reads")
	parser.add_argument("-R", "--reverse", nargs='*', help="Reverse reads")
	parser.add_argument("-b", "--bamfile", help="BAM infile")
	parser.add_argument("-t", "--threads", default=1, help="Threads")
	parser.add_argument("-s", "--snpsift", action="store_true", help="Run snpSift")
	parser.add_argument("-r", "--clean", action="store_true", help="Removes the SAM and BAM files")
	parser.add_argument("-p", "--done", action="store_true", help="Add an empty file to mark the directory as done")
#	args = parser.parse_args()
	return parser.parse_args(args)

class TestClass:
	def setup_args(self):
		return parse_args(['-f', '/proj/data11/vilma/Pipeline_vilma/example_data/test_Skeletonema_marinoi_Ref.txt', '-F', '/proj/data11/vilma/Pipeline_vilma/P8352_102/P8352_102_S1_L001_R1_001.fastq.gz', '-R', '/proj/data11/vilma/Pipeline_vilma/P8352_102/P8352_102_S1_L001_R2_001.fastq.gz'])
#		print vilmas_pipeline.args
#		pass
		
class TestProcess(TestClass, unittest.TestCase):

	def setUp(self):
		pass

	def test_bowtie2(self):
		self.args = self.setup_args()
#		self.bowtie2 = vilmas_pipeline.bowtie2
#		self.assertIsNone(self.bowtie2)
		print self.args 

#	def test_samtools_view(self):
#		cmd3 = ('samtools view -Sb %s > %s') % (sam, bam)
#		process3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE, shell=True)
#		stdout = process3.communicate()
#                self.assertIsNone(stdout)
#
#	def test_samtools_sort(self):
#		cmd4 = ['samtools', 'sort', bam, '-o', sorted_bam_out]
#		process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE)
#		stdout = process4.communicate()
#                self.assertIsNone(stdout)


#	def test_bcftools(self):
#		self.assertIsNotNone(vilmas_pipeline.bcftools)	

#	def test_annotation(self):
#		self.assertIsNotNone(vilmas_pipeline.annotation)	

#	def test_filtering(self):
#		self.assertIsNotNone(vilmas_pipeline.filtering)	

if __name__ == '__main__':
        unittest.main()

