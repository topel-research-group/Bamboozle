# Quantification of Bamboozle-derived amplicon-sequenced barcodes

1.) Open `Quantification.sh` in your text editor and enter the required information at the start of the script:
  * Absolute path to the directory containing your input data
  * Location of your reference FASTA file (formatted for dada2's assignSpecies function](http://benjjneb.github.io/dada2/training.html))
  * Adapter/primer sequences (and other parameters) for data trimming
  * Read merging program

If running this on the test dataset, the following variables should be set:
* `INDATA="${PWD}/sample_data/RawData/"`
* `REFFILE="sample_data/reference.fasta"`
* `ADAPTER="AGATCGGAAGAG"`
* `FWDPRIMER="AGGYTTCGCCTCCTCAAAC"`
* `REVPRIMER="GGCACGATGCACACGCAAAG"`
* `QUALITY=28`
* `MINLENGTH=180`
* `MERGE="bbmerge"`

2.) Ensure that this directory contains the required additional files
  * `index_file.lst` - A file where the first column contains the names of all samples to be analysed
  * `Allele_indexing.txt` - A tab-separated file containing information required for separating shared alleles
  * `Indexing.txt` - A tab-separated file containing information on the experimental setup

3.) Run `Quantification.sh`

The `sample_data` directory contains the following files, to demonstrate the functionality of this quantification pipeline:
* `RawData/` - raw reads from three *S. marinoi* mixed samples
* `reference.fasta` - reference FASTA for the *S. marinoi* experiment (see below)
* `Allele_indexing.txt` (see above)
* `index_file.lst` (see above)
* `Indexing.txt` (see above)

########

2023-02-16
Author: Björn Andersson

# Analysis of intra-specific metabarcoding using amplicon sequencing of C12W1 locus in *Skeletonema marinoi*

## Script and experimental summary

The `Barcode_analysis.R` script is a custom-made analysis pipeline of amplicon sequences of a hypervariable locus in *Skeletonema marinoi* (SM). 
The locus was bioinformaticaly identified based on analysis of whole genome sequences of 55 strains of SM from two Baltic Sea locations (Gåsfjärden: VG, and Gropviken: GP). 
It was predicted to have at least one unique allele for every strain used in the experiment.  
In summary the script requires sample data file(s) that have already been mapped against a reference database, i.e. numerical counts of observations of specific sequences. 
In our case this corresponds to counts of alleles of 59 different strains of SM (e.g. `C12W1_abundances/Input/P21502_101.tsv`). The `Indexing.txt` file contains metadata on each `.tsv` 
datafiles sample which comes from both control samples of sequences from individual genotypes, as well as mixed samples during a selection experiment. This experiment is described in 
detail in (Andersson et al. 2023, in prep for Molecular Ecology) but briefly, strains from the two populations were mixed at equal cell densities, grown semi-continuously in the 
exponental growth phase for 42 days (50-100 generations) with toxic copper stress (Cu) or without (C), and DNA samples where collected on day 0 (MasterMix, or MM), 9 (t9), or 42 (t42). 
The data files are then processed according to allelic information in `Allele_indexing.txt` file. This file is complex and not all information is essential to the current script 
function (See section *Input files and Metadata* below for more information). 

In short, the script reads the input file (`strain_table.tsv` in the dada2 output directory), extracts the sample information and matches it with the correct metadata. It then parses out 
non-unique allele observations using differential equations and the abundance of the secondary, unique, allele of strains with shared alleles. It then normalizes the allele observations 
to relative abundance using the sum of all amplicons matching known alleles (consequently, excluding all non-perfect matches). While doing this, it generates QC figures for each sample 
into the `02_plots/Plots` subdirectory, and outputs the normalized data into individual files in the `02_plots/Indexed` directory. It also outputs the total amplicon count into 
`02_plots/Read_counts`. The rest of the script explores the data by generating a number of graphs that are exported into the `02_plots/Plots` subdirectory.  
Note: As a warning, the script seemingly randomly crashes on my machine in some of of these loops, which is why Quartz is sometime called on and generating pop-up windows. 
A workaround is to disable the visualization of plots.

Note that all ASV sequences that do not by name match any of the list of alleles - `Indexing.txt` - will be ignored in the analysis. 
This can be modified by adding them to the indexing files.

![Experimetal Design1](https://github.com/Bearstar85/Cu_evolution/blob/master/ExperimentalDesign1.jpg)

## Input files and metadata

### **Indexing.txt**
This file contains metadata on the experimental samples, and contains the following columns:
* "ID" corresponds to the sample name and is used to join with "Barcode" in `strain_table.tsv`. 
* "Sample" has unique identifier information about the sample, e.g. if it's from a single genotype (e.g. *GP2-4_26*), or from a mixed selection sample (e.g. *PIPT_GP_Cu5_t9_C12W1*). 
* "Experiment" has category information on if it's a Strain/Genotype. 
* "Population" has category information on what population the sample is from (GP or VG). 
* "Treatment" has category information on treatment (Cu or C, NA for strains). 
* "Timepoint" has numerical information on time point in days (NA for strains). 
* "Bottel" has category information on unique identifier of each experimental bottle (e.g. VG_C1, NA for strains). 
* "Replicate" has integer information on replicate number per treatment-population combination (1-5, NA for strains).

### **Allele_indexing.txt**
Note that currently the **Barcode_analysis.R** script uses only information in columns: "Barcode", "Diffrential_ID", "Strain" and "Differential_equation". 
* "Barcode" contains unique identifier of a strain's alleles, with asterisks signifying when an allele is shared between 2 (one asterisk), or 3 strains (two asterisks), in the database. 
Used to join with `strain_table.tsv` file so this column must be identical to that for proper matching. 
* "Diffrential_ID" contains a short unique identifier for each strain's alleles (eg. ac, or aa1 and aa2 for a shared allele). 
* "Strain" contains information about the strains that contain the allele (e.g. strain VG1-2_94, or GP2-4_45and46 which are clones, two strains VG1-2_99or65 lack genotype sample but have 
been putatively linked to one out of two individuals). 
* "Differential_equation" contains the arithmetic formula that informs the script how observations of non-unique alleles are parsed based on abundance of other alleles in the sample. 
(e.g. ac for the unique allele ac, or aa1 * (an / (an + bt)) for the shared allele aa1, where an is the second allele of GP2-4_26 and bt the second one for GP2-4_71).

Other columns are non-essential but:
* "Uniqe_all" - Binary indexing if allele is unique amongst all alleles in database
* "Uniqe_pop" - Binary indexing if allele is unique amongst all alleles in individual population (GP or VG)
* "Homologs" - Alphabetical indexing of which alleles have homologs or are unique.
* "Index" - Conserved row number for sorting.

## **sample_data/reference.fasta**
Contains the database of known allel sequences


## Structure of scripts section by section

### Barcode_analysis.R, the main script

#### Houskeeping the loading packages
Lines 1-26 does housekeeping, loading packages, defining the input file - `strain_table.tsv` - and creating the relevant subdirectories.

#### Data normalization
Lines 29-37 read in the data.
* Note: certain Indexing names in Allels (column "Diffrential_ID") will not comply with differential equation function (notably x and df cannot be used).  
Lines 40-58 sets up the differential equation command.  
Lines 63-154 loop through each sample in the input file to append indexing, solve the differential equation, normalize the absolute amplicon observations to the relative abundance, 
output the indexed data to the `Index` subdirectory (note that CSV format is used with `.txt` file extension), write out a header, append the data to the dataframe `All_reads`, makes
QC plots for each file, and save them as PDFs in the `Plots` subdirectory.

#### Graphical interpretation of data
From here on the data only gets subselected, restructured, and plotted.  
Lines 202-300 focus on analysis of the time zero mastermix samples, and any false positive observations in them. Most relevant results are outputted as figures to the `Plots` subdirectory.  
Lines 337-357 plot how the individual strains' alleles are behaving during the selection experiment, and in the individual genotype samples. This loop is complex and prone to crashing, 
so there are multiple places where it can be cut short to avoid this by commenting or uncommenting a loop-ending `}`.  
However, the first lines of code in the loop **must** be executed for downstream parts of the script to function (specifically, the dataframe AlleleRatios needs to be filled on lines 313-352).

#### Allele ratios and correlations
This part of the script explores how the two alleles within a strain covary. This script works best if all strains have 2 alleles, but it will try and plot single alleles and should 
work on haploid or homologous strain data.  
Lines 433-575 contain lots of exploration into how the alleles covary, and what false positive signals we see in the data. Graphs are outputted into the `Plots/Strains` subdirectory. 
My personal favorite graph is `FigReg.pdf`.

#### Clustering of samples
This part of the script was used to track down deviating alleles, and cluster occurence patterns of unknown ASVs during trobleshooting. It is somewhat unfinished in part since it is  
largely redundant now.  
What it is doing is restructuring tha data to a Matrix (barcode x sample) that can be used for heatmap clustering, PCAs and correlation analysis. Given how clearly the strain-specific 
alleles covary now, and how agreeable replication samples are, this type of analysis did not yield much further insight for us, but it may be useful for troubleshooting new datasets in 
an iterative fashion in the known alleles ASVs, and also see how unknown ASVs where linked with certain strains (chimeras or common artifactual ASVs that slipped through the pipeline).

#### StrainCounts
This is a modified way of summing up the strain counts by their alleles. It can silence alleles prone to artefacts, and compensate for this through observations of the second allele. 
The input and output is smilar to the Data normalization step above, but uses other parts of the indexing file, e.g. the equations in the "Differential_equation_strains" column.
