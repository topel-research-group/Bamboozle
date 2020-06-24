#André Soares - June 2020

#Pseudocode to take in multiple .vcf inputs and metadata associated to them and plot:
#	0. check metadata table for requirements
#		a. check that vcf file sample name matches metadata table column
#		b. check if other columns exist
#		c. check for integrity of other columns in metadata table
#	1. associate vcfs to metadata (format? tsv, etc?)
#		a. take sample names from vcf filename
#		b. take sample name from metadata table column
#	2. produce table(s) of SVs, metadata per sample (used to plot 4,5)
#		a. number of SV types, avg length, per sample, per chromosome
#	4. produce plots of SV types per sample (w metadata?)
#	5. produce circos plot against ref genome?
