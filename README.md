# Bamboozle

Script for retrieving statistics from .bam files; in the first instance, for retrieving a statistic for what
percentage of bases in an assembly have >= Nx coverage

## v1

Unfinished implementation in Bash, before conversion to Python 2.x

## v2

Designed to look at each contig individually

## v3

Designed to look at the assembly as a whole; per contig functionality can be added later

## v4

Added per-contig functionality, with a flag to specify whether whole-genome or per-contig stats are required

## v5

Improved help screen; started adding helpful(?) error messages
* Attempting to implement ignored-contigs option, but experiencing some difficulties

## v5 Basic

Altering the script to, by default, assume a whole-genome coverage analysis with a threshold of 20 (unless
otherwise specified), to facilitate inclusion of a 'contigs to be ignored' flag
* Having problems implementing the 'default' value for threshold...


## Features to add

* A flag for ignoring specified contigs (for example, ignoring organellar contigs so that it doesn't
  distort results with super-high coverage)
* Add try/except statements
* Alter samtools variable so that it points to the samtools installation regardless of where it is
