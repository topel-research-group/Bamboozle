# Bamboozle

Script for retrieving statistics from .bam files; in the first instance, for retrieving a statistic for what
percentage of bases in an assembly have >= Nx coverage

## Features to add

* Allow -c flag to accept multiple contigs
* Add try/except statement for samtools



## Version history:
* v1 - Unfinished implementation in Bash, before conversion to Python 2.x
* v2 - Designed to look at each contig individually
* v3 - Designed to look at the assembly as a whole; per contig functionality can be added later
* v4 - Added per-contig functionality, with a flag to specify whether assembly or per-contig stats are required
* v5 - Improved help screen; started adding helpful(?) error messages
  * Attempting to implement ignored-contigs option, but experiencing some difficulties
* v5 Basic - Altering the script to, by default, assume a whole-genome coverage analysis with a threshold of 20
 (unless otherwise specified), to facilitate inclusion of a 'contigs to be ignored' flag
  * Having problems implementing the 'default' value for threshold...



## File to point tests at

/proj/data17/Skeletonema_marinoi_adaptation_to_warming_project/01_mapping/P8352_101/P8352_101_sorted.bam

Expected results (20x) - 97393434, 60.78250203191316%

Contigs to try: 000343F, 000111F-001-01
(Stats for 000343F (20x): 8932, 87.3152709359606%)
