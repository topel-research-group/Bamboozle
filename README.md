# Pipeline_vilma

This is a pipeline that will take raw data, such as fastq files, or a BAM file as input and give a filtered human readable output. The purpose is to simplify and make the bioinformatic analysis easier to reproduce.

The steps include mapping, SNP calling, Fst statistics, filtering and plotting of the results.  

## Features  
* This pipeline will give you a summary of the analysis in table format with effect prediction and annotated variants  
* Fst statistics on the results can be made with the additional program `fst.py`, this will give you a table with Fst statistics and a plot of the Fst results
* If you want to plot a specific contig use the additional program `make_plot.py`

## Usage  
`vilmas_pipeline.py -f <REFERENCE> <INPUT>`     

Input files:  
`-F <FORWARD READS> -R <REVERSE READS>` or `-b <BAMFILE>` 
 
####Example:  
Fastq files as input:
`vilmas_pipeline.py -f Skeletonema_marinoi_Ref_v1.1_Primary.all.fst -F file_R1.fastq.gz -R file_R2.fastq.gz`  

BAM file as input:   
`vilmas_pipeline.py -f Skeletonema_marinoi_Ref_v1.1_Primary.all.fst -b file.bam`

### Dependencies
```
Versions used:  
Bowtie2/v2.3.3.1  
samtools/v1.9  
bcftools/v1.9  
snpEff/v.4.3t

Versions used in fst.py and make_plot.py:
VCFtools/v0.1.13  
Highcharts/v6.2.0
```  

