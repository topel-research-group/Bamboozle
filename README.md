# Bamboozle

Script for retrieving statistics or other types of sequence information from .bam files. 



Retrieving a statistic for what percentage of bases in an assembly have >= Nx coverage

```bash
bamboozle.py -r <REFERENCE> -b <BAMFILE>
```

Extracting consensus sequence of aligned reads from a specific region of the referense sequence.

```bash
./bamboozled.py -r <REFERENCE> -b <BAMFILE> -c <CHROMOSOME/CONTIG> -a <RANGE>
```
```bash
./bamboozled.py -r Skeletonema_marinoi_Ref_v1.1.1.fst -b P8352_150_sorted.bam -c 000028F -a 686188-691148
```

## Features to add

* A flag for ignoring specified contigs (for example, ignoring organellar contigs so that it doesn't
  distort results with super-high coverage)
* Add try/except statements
* Alter samtools variable so that it points to the samtools installation regardless of where it is



## File to point tests at

/proj/data17/Skeletonema_marinoi_adaptation_to_warming_project/01_mapping/P8352_101/P8352_101_sorted.bam
