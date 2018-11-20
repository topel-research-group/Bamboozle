
import unittest
import vilmas_pipeline  
import subprocess #import Popen, PIPE
import os

class TestClass:
	
	def setup_args(self):
		#return vilmas_pipeline.args(['-f', '/proj/data11/vilma/Pipeline_vilma/example_data/test_Skeletonema_marinoi_Ref.txt', '-F', 'P8352_102_S1_L001_R1_001.fastq.gz', '-R', 'P8352_102_S1_L001_R2_001.fastq.gz'])
		print vilmas_pipeline.args
		
class TestProcess(TestClass, unittest.TestCase):

	def setUp(self):
		pass

	def test_bowtie2(self):
		self.args = self.setup_args()
		self.assertIsNone(vilmas_pipeline.bowtie2(self.args))	

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

