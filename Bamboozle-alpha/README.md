# Bamboozle

Script for retrieving statistics or other types of sequence information from .bam files. 

## Coverage statistics

Retrieving a statistic for what percentage of bases in an assembly have >= Nx coverage

```bash
./bamboozle.py --mode coverage -b <BAMFILE> -c <CHROMOSOME/CONTIG> -t <THRESHOLD>
```
```bash
./bamboozle.py --mode coverage -b P8352_101_sorted.bam -c 000343F -t 25
```
* `-c` flag is optional; if not specified, analysis will be run on the whole assembly
* `-t` flag is optional; if not specified, coverage of >= 20x will be reported
* Note: Requires a version of samtools which allows the `depth -aa` option
* Note: percentage values are rounded to 3 decimal places


## Consensus sequence

Extracting consensus sequence of aligned reads from a specific region of the reference sequence.

```bash
./bamboozle.py --mode consensus -r <REFERENCE> -b <BAMFILE> -c <CHROMOSOME/CONTIG> -a <RANGE>
```
```bash
./bamboozle.py --mode consensus -r Skeletonema_marinoi_Ref_v1.1.1.fst -b P8352_150_sorted.bam -c 000028F -a 686188-691148
```


## Zero coverage

Finding areas of zero coverage and printing the reference sequence, along with a GC percentage

```bash
./bamboozle.py --mode zero -r <REFERENCE> -b <BAMFILE> -c <CHROMOSOME/CONTIG>
```
```bash
./bamboozle.py --mode zero -r Skeletonema_marinoi_Ref_v1.1.1.fst -b P8352_101_sorted.bam -c 000343F
```
* Note: percentage values are rounded to 3 decimal places


## Identify deletions

Finding deletions; can print each deletion position, groups of adjacent deletion positions into events, or only
those events resulting in a frameshift
* Note: this can now also be done for a specific contig
* Note: this uses the default threshold of 20 for assessing the surrounding coverage; `-t` can therefore be used
  to try and find deletions in lower-coverage areas, however this hasn't been extensively tested

* Mode 1 - prints every deletion position
```bash
./bamboozle.py --mode deletion-1 -b <BAMFILE> [-c <CHROMOSOME/CONTIG>]
```

* Mode 2 - combine adjacent deletion positions into discrete events
```bash
./bamboozle.py --mode deletion-2 -b <BAMFILE> [-c <CHROMOSOME/CONTIG>]
```

* Mode 3 - only report events which represent a frameshift (i.e. not divisible by 3)
```bash
./bamboozle.py --mode deletion-3 -b <BAMFILE> [-c <CHROMOSOME/CONTIG>]
```

## Identify deletions within exons

Finding deletions occurring within exons; requires a bed file of exons and a text file of deletions (the output
of the above `deletion-2` and `deletion-3` modes)
* Note: Both `-x` and `-m` are required
* Note: This function will eventually be merged with the above `deletions` function

```bash
./bamboozle.py --mode deletion-x -x <EXONS> -m <MUTATIONS>
```

## Identify homozygous or heterozygous deletions

Helping to identify whether a deletion is homozygous or heterozygous, by comparing coverage between the deletion
and the previous base
* Note: This currently generates an (automatically-deleted) intermediate file
* Note: This function will eventually be merged with the above `deletions` function

```bash
./bamboozle.py --mode homohetero -b <BAMFILE> -m <MUTATIONS>
```

## Find regions which differ from the contig median by +/- 50%, and output them in .bed format

```bash
./bamboozle.py --mode median-one -c <CONTIG> -b <BAMFILE> > <output.bed>
```

## As above, but for each contig in the assembly

```bash
./bamboozle.py --mode median-all -b <BAMFILE> > <output.bed>
```

## Obtain a median per-contig coverage for the assembly

This is a stripped-down version of the `median-all` function, providing a simplified list of only the contig medians

```bash
./bamboozle.py --mode median_pc_coverage -b <BAMFILE> > <output.txt>
```

## Identify the longest stretch in a given contig which has coverage between two specified limits

```bash
./bamboozle.py --mode long-coverage -b <BAMFILE> -c <CONTIG> -l <LOWER LIMIT> <UPPER LIMIT>
```


## Features to add

* Add capacity to first run Bowtie2 mapping?
* Allow -c flag to accept multiple contigs (one dictionary per contig?)
  * Would therefore need to add additional arguments into both coverage-related functions...
* Move print statement out of coverage_stats function
* Refine deletion function to allow it to output differences between datasets, e.g. where
  deletions occur in one subset of the data but not another (e.g Warm vs. Cold results)
* Refine handling of borderline deletion cases, i.e. a three bp heterozygous event where the first
  base is just outside the boundary of being reported (see P8352_101 - 000202F:6,205-6,207)
* In exon_mutations function, find a way to pass the results of deletion function directly into
  exon_mutations, rather than requiring an intermediate file
* In HomoDel_or_Hetero function, find a way to pass the results of deletion function directly into
  HomoDel_or_Hetero, rather than requiring two intermediate files

* Any way to speed things up when searching a contig later in the .bam file?
  * Is speed also an issue when reading through the output of bedtools?
    * Try `mmap`? (https://www.quora.com/What-are-the-fastest-ways-to-read-text-lines-in-large-files-by-Python)

* SPEED-UP CAN BE ACHIEVED USING `samtools depth -aa -r`
  * Can't speed up functions which use `bedtoools genomecov`, e.g. `--mode zero`
    * Speed is definitely an issue here...
* Code can be simplified thanks to this flag...

* Add function for giving MEAN average, instead of/as well as median

## Known bugs
In the consensus sequence function, it is possible for two different comma-separated alternatives to be printed
together in the output
