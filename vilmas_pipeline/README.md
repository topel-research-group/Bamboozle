# Pipeline_vilma

**This is a pipeline that will take raw data, such as FASTQ files, or a BAM file as input and give a filtered human readable output. The purpose is to simplify and make the bioinformatic analysis easier to reproduce.**

The steps include mapping, SNP calling, Fst statistics, filtering and plotting of the results.  

<img src="https://user-images.githubusercontent.com/42669709/52559213-fbb93580-2df4-11e9-96f6-a8e6877352b6.png" width="400">

## Features  
* This pipeline will give you a summary of the analysis in table format with effect prediction and annotated variants  
* Fst can be calculated using the utility scripts `vcftools_fst.py` or `angsd_fst.py`, this will give you a table with Fst values and a plot of the results
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

|Program|Versions used|
|---|---| 
|`vilmas_pipeline` | Bowtie2/v2.3.3.1 <br/> samtools/v1.9 <br/> bcftools/v1.9 <br/> snpEff/v.4.3t <br/> gffutils/v0.9 <br/> Bedtools2/v2.27.1|
|`angsd_fst` |angsd/v0.918 |
|`vcftools_fst` | VCFtools/v0.1.13 |  
|`make_plot` | Highcharts/v6.2.0 |  

## More information  
* [Wiki](https://github.com/topel-research-group/Pipeline_vilma/wiki)  

