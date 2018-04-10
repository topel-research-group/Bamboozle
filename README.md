# Bamboozle

Script for retrieving statistics or other types of sequence information from .bam files. 


Retrieving a statistic for what percentage of bases in an assembly have >= Nx coverage

```bash
./bamboozle.py -b <BAMFILE> -c <CHROMOSOME/CONTIG> -t <THRESHOLD>
```
```bash
./bamboozle.py -b P8352_101_sorted.bam -c 000343F -t 25
```
* `-c` flag is optional; if not specified, analysis will be run on the whole assembly
* `-t` flag is optional; if not specified, coverage of >= 20x will be reported
* Note: Requires a version of samtools which allows the `depth -aa` option
* Note: percentage values are rounded to 3 decimal places

Extracting consensus sequence of aligned reads from a specific region of the referense sequence.

```bash
./bamboozle.py -r <REFERENCE> -b <BAMFILE> -c <CHROMOSOME/CONTIG> -a <RANGE>
```
```bash
./bamboozle.py -r Skeletonema_marinoi_Ref_v1.1.1.fst -b P8352_150_sorted.bam -c 000028F -a 686188-691148
```

Finding areas of zero coverage and printing the reference sequence, along with a GC percentage

```bash
./bamboozle.py -z -r <REFERENCE> -b <BAMFILE> -c <CHROMOSOME/CONTIG>
```
```bash
./bamboozle.py -z -r Skeletonema_marinoi_Ref_v1.1.1.fst -b P8352_101_sorted.bam -c 000343F
```
* Note: percentage values are rounded to 3 decimal places

Finding deletions; can print each deletion position, groups of adjacent deletion positions into events, or only
those events resulting in a frameshift
* Note: this can now also be done for a specific contig
* Note: this uses the default threshold of 20 for assessing the surrounding coverage; `-t` can therefore be used
  to try and find deletions in lower-coverage areas, however this hasn't been extensively tested

* Default usage - prints every deletion position
```bash
./bamboozle.py -d -b <BAMFILE> [-c <CHROMOSOME/CONTIG>]
```

* With `-e` flag - combine adjacent deletion positions into discrete events
```bash
./bamboozle.py -d -e -b <BAMFILE> [-c <CHROMOSOME/CONTIG>]
```

* With `-f` flag (or `-e -f`) - only report events which represent a frameshift (i.e. not divisible by 3)
```bash
./bamboozle.py -d [-e] -f -b <BAMFILE> [-c <CHROMOSOME/CONTIG>]
```

Finding deletions occuring within exons; requires a bed file of exons and a text file of deletions (the output
of the above `-d -e/-f` function)
* Note: Both `-x` and `-m` are required
* Note: This function will eventually be merged with the above `deletions` function

```bash
./bamboozle.py -x <EXONS> -m <MUTATIONS>
```



## File to point tests at

example/P8352_101_sorted.bam

Expected results (20x) - 97393434, 60.78250203191316%

Contigs to try: 000343F, 000111F-001-01
(Stats for 000343F (20x): 8932, 87.3152709359606%)

## Features to add

* Allow -c flag to accept multiple contigs (one dictionary per contig?)
  * Would therefore need to add additional arguments into both coverage-related functions...
* Move print statement out of coverage_stats function
* Refine deletion function to allow it to output differences between datasets, e.g. where
  deletions occur in one subset of the data but not another (e.g Warm vs. Cold results)
* Refine handling of borderline deletion cases, i.e. a three bp heterozygous event where the first
  base is just outside the boundary of being reported (see P8352_101 - 000202F:6,205-6,207)
* In exon_mutations function, find a way to pass the results of deletion function directly into
  exon_mutations, rather than requiring an intermediate file

## Existing codes to incorporate

* HomoDel_or_Hetero.py
  * Compares coverage at the first base of the deletion with the previous base; from the percentage difference
    one can determine whether this is a homozygous or heterozygous deletion
  * Input is the result of `samtools depth -aa -b [frameshift locations].bed [sample].bam > [CODE INPUT]`
    * `samtools depth` accepts as input a list of bams instead, try this for working with multiple samples,
       e.g. all Warm samples
    * This file generates a tab-separated file in the format - contig, position, coverage depth
      * Note: If two side-by-side mutations are noted in the initial bed file, entries in the new file won't
        be duplicated, i.e. if there is a mutation at pos. 2 in some samples and at pos. 3 in others, instead
        of showing coverage 1, 2, 2, 3, it will just show 1, 2, 3; this has been addressed in the code, however
  * Output is a list of deletion positions and their percentage vs. the previous position
  * `./HomoDel_or_Hetero.py -i [input].txt > [output]`
