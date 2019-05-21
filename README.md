# Merging pipelines

## Alterations

* Adjusted the time-keeping for the `--dev` flag by using `time()` instead of `time.time()`
  * Also moved the function from `bamparser.py` to `main.py`


## To do

* Implement additional checks for presence of snpEff and Java
  * See similar implementation for samtools

* Is it necessary to run snpEff for every function of `bamparser.py` that uses the `reference` flag?

* `Bowtie2` and `Bcftools` directories don't need to be created if the script doesn't explicitly require them

* `--homohetero` function bugs out if there is only one event in a given contig (`next` never assigned)

* Way to eliminate the need for the `--bamparse` flag?
  * Specify list of flags, combined with an `if ... in ...` argument?

* Play with groupings in help file
