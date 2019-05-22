# Utility scripts  
These scripts will do Fst statistics and plot the results from the pipeline.  
See additional options and information in the [Manual](https://github.com/topel-research-group/Bamboozle/wiki/Manual/)  in the Wiki.

## Usage    
`vcftools_fst.py --pop1 <INPUT POPULATION1> --pop2 <INPUT POPULATION2>`   
Where INPUT is a directory with subdirectories corresponding to the samples. The script uses [VCFtools](https://vcftools.github.io), vcf.gz files as input.   

`angsd_fst.py --pop1 <INPUT POPULATION1> --pop2 <INPUT POPULATION2> --ref <FASTA>`  
Where INPUT is a directory with subdirectories corresponding to the samples and FASTA a reference genome in FASTA format. The script uses [ANGSD](https://github.com/ANGSD/angsd), bam files as input.  

`make_plot.py --inputfile <INPUTFILE> --grep <"NAME">`   
Where INPUTFILE is the resulting csv file from `vcftools_fst.py` or `angsd_fst.py` and "NAME" is the name of the contig you want to plot.  

`kmer_plot.py --inputfile <INPUTFILE> --grep <"NAME"> --length <int>`  
Where INPUTFILE is a csv file with information about kmer peaks with Contig,Position,Value as header.  
NAME is the name of the contig you want to plot and length is the length of the contig.
  
### Examples:  
`vcftools_fst.py --pop1 Pop_1/ --pop2 Pop_2/`   

Ouput: a table and a csv file with the names, positions and Fst value and a summary of the results in a pdf.  

`angsd_fst.py --pop1 Pop_1/ --pop2 Pop_2/ --ref example.fst`  

Ouput: a table and a csv file with the names, positions and Fst value and a summary of the results in a png image.

`make_plot.py --inputfile pop1_pop2_results.csv --grep "contig_name"`   

Ouput: a html file with an interactive plot of the selected contig.  

`make_plot_wholefile_3col.py --inputfile pop1_pop2_results.csv --yaxis Value`   

Ouput: a html file with an interactive from the csv file provided (If the file is big, use `make_plot.py` or `pandas_fst_plot.py` instead).  

`pandas_fst_plot.py --inputfile pop1_pop2_results.csv --yaxis Value`

Output: a png with the results from the csv file provided.

`kmer_plot.py --inputfile KmerPeakPositions.csv --grep "contig_name" --length 1592421`  

Ouput: a html file with an interactive plot of the selected contig.   


### Dependencies  
VCFtools/v0.1.13   
Highcharts/v6.2.0   
angsd/v0.918  

