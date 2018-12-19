# Utility scripts  
These scripts will do Fst statistics and plot the results from the pipeline.  

## Usage    
`fst.py --pop1 <INPUT POPULATION1> --pop2 <INPUT POPULATION2>`   
Where INPUT is a directory with subdirectories corresponding to the samples.  

`angsd.py --pop1 <INPUT POPULATION1> --pop2 <INPUT POPULATION2> --ref <FASTA>`  
Where INPUT is a directory with subdirectories corresponding to the samples and FASTA a reference genome in FASTA format.

`make_plot.py --inputfile <INPUTFILE> --grep <"NAME">`   
Where INPUTFILE is the resulting csv file from `fst.py` or `angsd.py` and "NAME" is the name of the contig you want to plot.  

`kmer_plot.py --inputfile <INPUTFILE> --grep <"NAME"> --length <int>`  
Where INPUTFILE is a csv file with information about kmer peaks with Contig,Position,Value as header.  
NAME is the name of the contig you want to plot and length is the length of the contig.
  
#### Example:  
`fst.py --pop1 Pop_1/ --pop2 Pop_2/`   

Ouput: a table and a csv file with the names, positions and Fst value and a summary of the results in a pdf.  

`angsd.py --pop1 Pop_1/ --pop2 Pop_2/ --ref example.fst`  

Ouput: a table and a csv file with the names, positions and Fst value (pdf summary planned).

`make_plot.py --inputfile pop1_pop2_results.csv --grep "contig_name"`   

Ouput: a html file with an interactive plot of the selected contig.  

`kmer_plot.py --inputfile KmerPeakPositions.csv --grep "contig_name" --length 1592421`  

Ouput: a html file with an interactive plot of the selected contig.   


### Dependencies  
VCFtools/v0.1.13   
Highcharts/v6.2.0   
angsd/v0.918  

