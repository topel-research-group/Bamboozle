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



## File to point tests at

example/P8352_101_sorted.bam

Expected results (20x) - 97393434, 60.78250203191316%

Contigs to try: 000343F, 000111F-001-01
(Stats for 000343F (20x): 8932, 87.3152709359606%)

## Features to add

* Allow -c flag to accept multiple contigs (one dictionary per contig?)
  * Would therefore need to add additional arguments into both coverage-related functions...
* Move print statement out of coverage_stats function
* Limit number of decimal places in GC% output to make the results tidier