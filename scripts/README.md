# Utility scripts  
These scripts will do Fst statistics and plot the results from the pipeline.  

## Usage    
`fst.py --pop1 <INPUT POPULATION1> --pop2 <INPUT POPULATION2>`   
Where INPUT is a directory with subdirectories corresponding to the samples.  

`make_plot.py --inputfile <INPUTFILE> --grep <"NAME">`   
Where INPUTFILE is the resulting csv file from the `fst.py` script and "NAME" is the name of the contig you want to plot.  

  
#### Example:  
`fst.py --pop1 Pop_1/ --pop2 Pop_2/`   

Ouput: a table and a csv file with the names, positions and Fst value and a summary of the results in a pdf.  


`make_plot.py --inputfile pop1_pop2_results.csv --grep "contig_name"`   

Ouput: a html file with an interactive plot of the selected contig.   


### Dependencies  
VCFtools/v0.1.13   
Highcharts/v6.2.0   
angsd/v0.918  

