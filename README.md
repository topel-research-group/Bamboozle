# Bamboozle

**This is a pipeline that will take raw data, such as FASTQ files, or a BAM file as input and give a filtered human readable output. The purpose is to simplify and make the bioinformatic analysis easier to reproduce.**

The steps include mapping, SNP calling, Fst statistics, filtering and plotting of the results.  

<img src="https://user-images.githubusercontent.com/42669709/52559213-fbb93580-2df4-11e9-96f6-a8e6877352b6.png" width="400">

**The pipeline is also able to obtain various coverage-related statistics from read mapping results**

## Features  
* The pipeline will give you a summary of the analysis in table format with effect prediction and annotated variants  
* Fst can be calculated using the utility scripts `vcftools_fst.py` or `angsd_fst.py`, this will give you a table with Fst values and a plot of the results
* If you want to plot a specific contig use the additional program `make_plot.py`

## Basic usage examples

For all optional arguments, see the [Bamboozle wiki](https://github.com/topel-research-group/Bamboozle/wiki) or use `bamboozle.py <command> -h`

|                   Function                   |                                          Example syntax                                           |                                              Notes                                              |
|----------------------------------------------|---------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
| FASTQ input                                  | `bamboozle.py -f <REFERENCE> -F <FORWARD READS> -R <REVERSE READS>`                               | Run pipeline from beginning and align reads using Bowtie2                                                                               |
| (Sorted) BAM input                           | `bamboozle.py -f <REFERENCE> -b <BAMFILE>`                                                        | Skips Bowtie2 and SAMtools, starts at BCFtools for SNP calling                                                                          |
| Coverage statistics                          | `bamboozle.py coverage {-f <ref> -F <fwd> -R <rev> | -b <bam>} [options]`                         | Find the percentage of bases in an assembly with >= Nx coverage; can also output this information (and overlapping genes) in a BED file |
| Consensus sequence                           | `bamboozle.py consensus -f <ref> {-F <fwd> -R <rev> | -b <bam>} -c <contig> -a <range> [options]` | Extracting consensus sequence of aligned reads from a specific region of the reference sequence                                         |
| Zero coverage                                | `bamboozle.py zero -f <ref> {-F <fwd> -R <rev> | -b <bam>} -c <contig> [options]`                 | Finding areas of zero coverage and printing the reference sequence, along with a GC percentage                                          |
| Identify deletions 1                         | `bamboozle.py deletion1 {-f <ref> -F <fwd> -R <rev> | -b <bam>} [options]`                        | Identify deletions and print every deletion position                                                                                    |
| Identify deletions 2                         | `bamboozle.py deletion2 {-f <ref> -F <fwd> -R <rev> | -b <bam>} [options]`                        | Identify deletions and combine adjacent positions into discrete events                                                                  |
| Identify deletions 3                         | `bamboozle.py deletion3 {-f <ref> -F <fwd> -R <rev> | -b <bam>} [options]`                        | Identify deletions and only report events resulting in a frameshift (i.e. not divisible by 3)                                           |
| Idenitfy deletions within exons              | `bamboozle.py deletionx {-f <ref> -F <fwd> -R <rev> | -b <bam>} -x <bed> [options]`               | Identify deletions events occuring within exons                                                                                         |
| Identify homozygous/heterozygous deletions   | `bamboozle.py homohetero {-f <ref> -F <fwd> -R <rev> | -b <bam>} [options]`                       | Compares coverage between a deletion and the previous position to try and determine whether a deletion is homozygous or heterozygous    |
| Identify deviation from median (complex)     | `bamboozle.py median {-f <ref> -F <fwd> -R <rev> | -b <bam>} --complex [options]`                 | Find regions which differ from the contig median by +/- 50%, and output them in .bed format                                             |
| Identify deviation from median (simple)      | `bamboozle.py median {-f <ref> -F <fwd> -R <rev> | -b <bam>} --simple [options]`                  | Output the median for one or all contigs in .txt format                                                                                 |
| Identify longest stretch of a given coverage | `bamboozle.py long_coverage {-f <ref> -F <fwd> -R <rev> | -b <bam>} -c <contig> -l <upper limit> <lower limit> [options]` | Identify the longest stretch in a given contig which has coverage between two specified limits                                          |
| Identify potential barcode regions           | `bamboozle.py barcode -f <ref> -b <bam1> <bam2> ... <bamN> -o <prefix> [options]`                 | Output BED and TXT files of coordinates for potential barcoding regions in the specified BAM files                                      |
| Run loss-of-function pipeline                | `bamboozle.py lof [ADD OPTIONS HERE]`                                                             | [ADD DESCRIPTION HERE]                                                                                                                  |

***  
| Utility script | Description |  
|---|---|  
|[`angsd_fst.py`](https://github.com/topel-research-group/Bamboozle/wiki/Manual)| Fst statistics using BAM file from the SAMtools step as input. Uses genotype likelihoods to calculate Fst values, preferred if coverage is low or medium.|  
|[`vcftools_fst.py`](https://github.com/topel-research-group/Bamboozle/wiki/Manual)| Fst statistics using VCFtools --weir-fst-pop takes vcf files from the SnpEff step as input.|  
|[`make_plot.py`](https://github.com/topel-research-group/Bamboozle/wiki/Manual)| Makes an interactive plot of either the ANGSD csv results or the VCFtools csv results by specifying the contig you want to plot. Outputs a html file.|

## Dependencies

|      Program       |                                                                  Versions used                                                                   |
|--------------------|--------------------------------------------------------------------------------------------------------------------------------------------------| 
|`Bamboozle`         | Bowtie2/v2.3.3.1 <br/> samtools/v1.9 <br/> bcftools/v1.9 <br/> snpEff/v.4.3t <br/> gffutils/v0.9 <br/> Bedtools2/v2.27.1 <br/> Java (for snpEff) |
|(`barcodesearch.py`)| NumPy and Levenshtein Python libraries                                                                                                           |
|`angsd_fst`         | angsd/v0.918                                                                                                                                     |
|`vcftools_fst`      | VCFtools/v0.1.13                                                                                                                                 |  
|`make_plot`         | Highcharts/v6.2.0                                                                                                                                |  

## More information  
* [Wiki](https://github.com/topel-research-group/Bamboozle/wiki)  


## Additional notes
* `coverage`
  * `-c` flag is optional; if not specified, analysis will be run on the whole assembly
  * `-d` flag is optional; if not specified, coverage of >= 20x will be reported
  * If the `-o` flag is specified, a BED file of regions with coverage above the threshold will be generated
    * The `--gff` flag can also be specified, which labels the above BED file with any gene models overlapping the high-coverage windows
  * Requires a version of samtools which allows the `depth -aa` option
  * Percentage values are rounded to 3 decimal places

* `consensus`
  * This function needs to be fixed...

* `zero`
  * Percentage values are rounded to 3 decimal places

* `deletion1`, `deletion2`, `deletion3`, `deletionx` and `homohetero`
  * `-c` flag is optional; if not specified, analysis will be run on the whole assembly
  * These functions use the default threshold of 20 for assessing the surrounding coverage;
    `-t` can therefore be used to try and find deletions in lower-coverage areas,
    however this hasn't been extensively tested

* `deletionx`
  * `-x` flag - a .bed file of exons - is required

* `homohetero`
  * Currently generates an (automatically-deleted) intermediate file

* `barcode`
  * Note that the bed file output uses 0-based coordinates, so the start positions are one lower
    than the actual start positions; the coordinates in the txt file output are correct
  * In addition to the bed and txt file output, this function also generates bgzipped vcf files
    and their respective index files

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
  * Can't speed up functions which use `bedtools genomecov`, e.g. `zero`
    * Speed is definitely an issue here...
* Code can be simplified thanks to this flag...

* Add function for giving MEAN average, instead of/as well as median

* `homohetero` function bugs out if there is only one event in a given contig (`next` never assigned)

* Add mutual exclusivity to FASTQ/BAM input

## Known bugs
* In the consensus sequence function, it is possible for two different comma-separated alternatives to be printed
  together in the output

* The homohetero function errors out if only one event occurs in a contig (as the `next` variable is never assigned)
