#André Soares - June 2020

#Pseudocode/notes to write code to take in multiple .vcf inputs and metadata associated to them and plot useful stuff:
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

#	How about an interactive shiny app? or html being made available with all this in it?
#		Sample-specific analysis via circos, other plots of interesting stuff
#		General stats
#	Could be a server-based solution, if circos plots or other results hard to generate
#	
#	OUTPUTS:
#		1. TABLE
#			COLLATE? VCF and add sample-specific metadata?
#			a.
#		2. INTERACTIVE HTML OUTPUT
#
#			
#			DASH BIO FOR CIRCOS
#
#		PLOT 1 - SV type length violin/point plot:
#		 - OPTION to see by sample and/or by group
#		 - OPTION to generate sample- or faceted group-level plot
#			a. x = chr
#			b. y = length
#			c. colour = sv_type
#			d. face
#		PLOT 2 - Circos plot:
#		 - OPTION to see by sample
#		 - OPTION to see x many groups faceted
#		(consider only first ANN in VCF?)
#			L1: chromosome
#			L2: indels (<=50bp), points?
#			L3: SVs (>50bp), bars?
#			L4: LOF Potential (bars, colours)
#			L5: LOF Annotation (bars, colours)
#coverage?
#			Tracks:
