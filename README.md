# Merging pipelines

## Alterations

* Adjusted the time-keeping for the `--dev` flag by using `time()` instead of `time.time()`
  * Also moved the function from `bamparser.py` to `main.py`

* Added a check for `bamparser.py` flags to eliminate the need for a `--bamparse` flag
  * There's probably a more elegant method, but it's functional
  * This also removes the ambiguity of the `-b` flag for (unsorted) bam input

* Add try/except statements to `pipeline.py` to ensure that, when snpEff is called,
  both Java and snpEff itself are loaded
  * This may need to be moved into the `annotation()` function?

* A few minor alterations to some of the `bamparser.py` functions, for clarity

* `bamparser` can now accept either fastq or unsorted bam files as input

* `args.ref` requirement removed from the middle clause of `main.py`, so that snpEff doesn't *always* run
  * Return it if required for the non-`bamparse` functions

* Added a `glob` function to determine newly-sorted bam names
  * Improves results display of `bamparse` functions

* Added a test to ensure no bam files are in Bowtie2 before starting the process,
  as this could cause file name confusion

## To do

* Does the `samtools_view()` function get `args.threads` passed to it?

* `--homohetero` function bugs out if there is only one event in a given contig (`next` never assigned)

* Is `sortbam = args.sortbam` redundant? Can `args.sortbam` be passed directly to functions?

* Play with groupings in help file
