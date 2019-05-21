# Merging pipelines

## Alterations

* Adjusted the time-keeping for the `--dev` flag by using `time()` instead of `time.time()`
  * Also moved the function from `bamparser.py` to `main.py`

* Added a check for `bamparser.py` flags to eliminate the need for a `--bamparse` flag
  * There's probably a more elegant method, but it's functional
  * This also removes the ambiguity of the `-b` flag for (unsorted) bam input

* Add try/except statements to `pipeline.py` to ensure that, when snpEff is called,
  both Java and snpEff itself are loaded

* A few minor alterations to some of the `bamparser.py` functions, for clarity

## To do

* Implement additional checks for presence of snpEff and Java
  * See similar implementation for samtools

* Is it necessary to run snpEff for every function of `bamparser.py` that uses the `reference` flag?

* `Bowtie2` and `Bcftools` directories don't need to be created if the script doesn't explicitly require them

* `--homohetero` function bugs out if there is only one event in a given contig (`next` never assigned)

* Play with groupings in help file

* Enable `bamparser.py` functions to accept fastq or unsorted bam input

* May need to add `shell=True` to get some of the mapping functions to work
