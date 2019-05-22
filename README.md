# Bamboozle

**This is a pipeline that will take raw data, such as FASTQ files, or a BAM file as input and give a filtered human readable output. The purpose is to simplify and make the bioinformatic analysis easier to reproduce.**

The steps include mapping, SNP calling, Fst statistics, filtering and plotting of the results.  

<img src="https://user-images.githubusercontent.com/42669709/52559213-fbb93580-2df4-11e9-96f6-a8e6877352b6.png" width="400">

**The pipeline is also able to obtain various coverage-related statistics from read mapping results**

## Features  
* The pipeline will give you a summary of the analysis in table format with effect prediction and annotated variants  
* Fst can be calculated using the utility scripts `vcftools_fst.py` or `angsd_fst.py`, this will give you a table with Fst values and a plot of the results
* If you want to plot a specific contig use the additional program `make_plot.py`

## Usage examples
|                   Function                   |                                        Example syntax                                         |                                              Notes                                              |
|----------------------------------------------|-----------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
| FASTQ input                                  | `bamboozle.py -f <REFERENCE> -F <FORWARD READS> -R <REVERSE READS>`                           | Run pipeline from beginning and align reads using Bowtie2
| (Unsorted) BAM input                         | `bamboozle.py -f <REFERENCE> -b <BAMFILE>`                                                    | Skips aligning step and starts from the SAMtools step
| (Sorted) BAM input                           | `bamboozle.py -f <REFERENCE> --sortedbam <SORTED BAMFILE>`                                    | Skips Bowtie2 and SAMtools, starts at BCFtools for SNP calling
| Coverage statistics                          | `bamboozle.py --coverage --sortbam <BAMFILE> [-c <CONTIG>] [-d <THRESHOLD>]`                  | Retrieving a statistic for what percentage of bases in an assembly have >= Nx coverage          |
| Consensus sequence                           | `bamboozle.py --consensus -f <REFERENCE> --sortbam <BAMFILE> -c <CONTIG> -a <RANGE>`          | Extracting consensus sequence of aligned reads from a specific region of the reference sequence |
| Zero coverage                                | `bamboozle.py --zero -f <REFERENCE> --sortbam <BAMFILE> -c <CONTIG>`                          | Finding areas of zero coverage and printing the reference sequence, along with a GC percentage  |
| Identify deletions 1                         | `bamboozle.py --deletion1 --sortbam <BAMFILE> [-c <CONTIG>]`                                  | Identify deletions and print every deletion position                                            |
| Identify deletions 2                         | `bamboozle.py --deletion2 --sortbam <BAMFILE> [-c <CONTIG>]`                                  | Identify deletions and combine adjacent positions into discrete events                          |
| Identify deletions 3                         | `bamboozle.py --deletion3 --sortbam <BAMFILE> [-c <CONTIG>]`                                  | Identify deletions and only report events resulting in a frameshift (i.e. not divisible by 3)   |
| Idenitfy deletions within exons              | `bamboozle.py --deletionx --sortbam <BAMFILE> -x <EXONS> [-c <CONTIG>]`                       | Identify deletions events occuring within exons                                                 |
| Identify homozygous/heterozygous deletions   | `bamboozle.py --homohetero --sortbam <BAMFILE> [-c <CONTIG>]`                                 | Compares coverage between a deletion and the previous position to try and determine whether a deletion is homozygous or heterozygous |
| Identify deviation from median (complex)     | `bamboozle.py --median --complex [-c <CONTIG>] --sortbam <BAMFILE> > <output.bed>`            | Find regions which differ from the contig median by +/- 50%, and output them in .bed format     |
| Identify deviation from median (simple)      | `bamboozle.py --median --simple [-c <CONTIG>] --sortbam <BAMFILE> > <output.txt>`             | Output the median for one or all contigs in .txt format                                         |
| Identify longest stretch of a given coverage | `bamboozle.py --long_coverage --sortbam <BAMFILE> -c <CONTIG> -l <LOWER LIMIT> <UPPER LIMIT>` | Identify the longest stretch in a given contig which has coverage between two specified limits |

***  
| Utility script | Description |  
|---|---|  
|[`angsd_fst.py`](https://github.com/topel-research-group/Bamboozle/wiki/Manual)| Fst statistics using BAM file from the SAMtools step as input. Uses genotype likelihoods to calculate Fst values, preferred if coverage is low or medium.|  
|[`vcftools_fst.py`](https://github.com/topel-research-group/Bamboozle/wiki/Manual)| Fst statistics using VCFtools --weir-fst-pop takes vcf files from the SnpEff step as input.|  
|[`make_plot.py`](https://github.com/topel-research-group/Bamboozle/wiki/Manual)| Makes an interactive plot of either the ANGSD csv results or the VCFtools csv results by specifying the contig you want to plot. Outputs a html file.|

## Dependencies

|     Program      |                                                                  Versions used                                                                   |
|------------------|--------------------------------------------------------------------------------------------------------------------------------------------------| 
|`Bamboozle` | Bowtie2/v2.3.3.1 <br/> samtools/v1.9 <br/> bcftools/v1.9 <br/> snpEff/v.4.3t <br/> gffutils/v0.9 <br/> Bedtools2/v2.27.1 <br/> Java (for snpEff) |
|`angsd_fst`       | angsd/v0.918                                                                                                                                     |
|`vcftools_fst`    | VCFtools/v0.1.13                                                                                                                                 |  
|`make_plot`       | Highcharts/v6.2.0                                                                                                                                |  

## More information  
* [Wiki](https://github.com/topel-research-group/Bamboozle/wiki)  


## Additional notes
* `--coverage`
  * `-c` flag is optional; if not specified, analysis will be run on the whole assembly
  * `-d` flag is optional; if not specified, coverage of >= 20x will be reported
  * Requires a version of samtools which allows the `depth -aa` option
  * Percentage values are rounded to 3 decimal places

* `--consensus`
  * This function needs to be fixed...

* `--zero`
  * Percentage values are rounded to 3 decimal places

* `--deletion1`, `--deletion2`, `--deletion3`, `--deletionx` and `--homohetero`
  * `-c` flag is optional; if not specified, analysis will be run on the whole assembly
  * These functions use the default threshold of 20 for assessing the surrounding coverage;
    `-t` can therefore be used to try and find deletions in lower-coverage areas,
    however this hasn't been extensively tested

* `--deletionx`
  * `-x` flag - a .bed file of exons - is required

* `--homohetero`
  * Currently generates an (automatically-deleted) intermediate file


## Features to add/things to address

* Allow -c flag to accept multiple contigs (one dictionary per contig?)
  * Would therefore need to add additional arguments into both coverage-related functions...

* Move print statement out of coverage_stats function

* Refine deletion function to allow it to output differences between datasets, e.g. where deletions
  occur in one subset of the data but not another (e.g between two different environmental conditions)

* Refine handling of borderline deletion cases
  * e.g. a three bp heterozygous event where the first base is just outside the boundary of being reported

* In HomoDel_or_Hetero function, find a way to eliminate the need for an intermediate file

* Any way to speed things up when searching a contig later in the .bam file?
  * Is speed also an issue when reading through the output of bedtools?
    * Try `mmap`? (https://www.quora.com/What-are-the-fastest-ways-to-read-text-lines-in-large-files-by-Python)

* SPEED-UP CAN BE ACHIEVED USING `samtools depth -aa -r`
  * Can't speed up functions which use `bedtools genomecov`, e.g. `--zero`
    * Speed is definitely an issue here...
* Code can be simplified thanks to this flag...

* Add function for giving MEAN average, instead of/as well as median

* Should the try/except statements for Java and snpEff in `pipeline.py` be moved into `annotation()`?

* Does the `samtools_view()` function get `args.threads` passed to it?

* `--homohetero` function bugs out if there is only one event in a given contig (`next` never assigned)

* Play with groupings in help file

* Add a flag to include the `snpEff` functions in the `bamparse` parts of the pipeline?

## Known bugs
* In the consensus sequence function, it is possible for two different comma-separated alternatives to be printed
  together in the output

* The homohetero function errors out if only one event occurs in a contig (as the `next` variable is never assigned)
