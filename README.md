# Bamboozle

Script for retrieving statistics from .bam files; in the first instance, for retrieving a statistic for what
percentage of bases in an assembly have >= Nx coverage

# v2

Designed to look at each contig individually

# v3

Designed to look at the assembly as a whole; per contig functionality can be added later

# v4

Added per-contig functionality, with a flag to specify whether whole-genome or per-contig stats are required

# Features to add

* A flag for ignoring specified contigs (for example, ignoring organellar contigs so that it doesn't
  distort results with super-high coverage)
* Add try/except statements
* Alter samtools variable so that it points to the samtools installation regardless of where it is
