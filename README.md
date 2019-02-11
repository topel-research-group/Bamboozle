# Pipeline_vilma

This is a pipeline that will take raw data, such as FASTQ files, or a BAM file as input and give a filtered human readable output. The purpose is to simplify and make the bioinformatic analysis easier to reproduce.

The steps include mapping, SNP calling, Fst statistics, filtering and plotting of the results.  

## Features  
* This pipeline will give you a summary of the analysis in table format with effect prediction and annotated variants  
* Fst statistics on the results can be made with the additional program `fst.py`, this will give you a table with Fst statistics and a plot of the Fst results
* If you want to plot a specific contig use the additional program `make_plot.py`

## Usage  
`vilmas_pipeline.py -f <REFERENCE> <INPUT>`     

Input files:  
`-F <FORWARD READS> -R <REVERSE READS>` or `-b <BAMFILE>` 
   
#### Example:  
FASTQ files as input:  
`vilmas_pipeline.py -f Skeletonema_marinoi_Ref_v1.1_Primary.all.fst -F file_R1.fastq.gz -R file_R2.fastq.gz`  

BAM file as input:   
`vilmas_pipeline.py -f Skeletonema_marinoi_Ref_v1.1_Primary.all.fst -b file.bam`
***  
| Utility script | Description |  
|---|---|  
|[`angsd_fst.py`](https://github.com/topel-research-group/Pipeline_vilma/wiki/Manual)| Fst statistics using BAM file from the SAMtools step as input. Uses genotype likelihoods to calculate Fst values, preferred if coverage is low or medium.|  
|[`vcftools_fst.py`](https://github.com/topel-research-group/Pipeline_vilma/wiki/Manual)| Fst statistics using VCFtools --weir-fst-pop takes vcf files from the SnpEff step as input.|  
|[`make_plot.py`](https://github.com/topel-research-group/Pipeline_vilma/wiki/Manual)| Makes an interactive plot of either the ANGSD csv results or the VCFtools csv results by specifying the contig you want to plot. Outputs a html file.|

## Dependencies

**Versions used:**  
Bowtie2/v2.3.3.1  
samtools/v1.9  
bcftools/v1.9  
snpEff/v.4.3t  
gffutils/v0.9  

**Versions used in utility scripts:**  
VCFtools/v0.1.13  
Highcharts/v6.2.0  
angsd/v0.918  

## More information  
* [Wiki](https://github.com/topel-research-group/Pipeline_vilma/wiki)  

