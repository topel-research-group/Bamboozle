#André Soares - 15/04/2020
#Pseudocode for a script calling SVs from BAM/SAM files
#very much based on Matt's(?) pipeline.py in Bamboozle/modules

#Input: BAM/SAM alignments of diatom (S. marinoi and others in the future?) genomes against a reference
#Output: 
#	1) VCF files of alterations in gene models, with LOF prediction.
#	2) summary tables of alterations that should include type (deletion, insertion, inversion, duplication, other genomic rearrangmements),
#size in bp, quality of call, LOF potential, gene product, genomic location, zygosity, ?...

#1. Take in mapped genomes, verify they exist, are sorted BAM files
#1a. Take in fasta, GFF, Pilon-corrected assembly references as arguments
#2. Run GRIDSS, make parameters available or define ideal/minimal parameters (do they even exist?)
#  (test run GRIDSS first, make sure outputs are immediately compatible with next step)
#3. Mask VCF files of SV calls with bedtools using Pilon assembly
#4. Run snpEff on masked SV calls to obtain predicted alterations in gene models
#5. Run snpEff snpSift to obtain summary tables

