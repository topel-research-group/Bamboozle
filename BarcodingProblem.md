# Problem
Current implementation of Bamboozle barcode function doesn't work

* Only takes homozygotes into account?
1. Phasing of Illumina reads?
  * Leverage paired-end data somehow...
  * `samtools phase`?
    * Test output - will need to rerun BCFtools to generate consensuses...
      * View parent and two child BAMs in IGV
      * Try a few values of the `-k` flag?
    * Phased alleles in dictionary, and identify unique alleles across all strains
    * Compare hashes of the variants? E.g. hash of '35T 467G'

`samtools phase` appears to be giving good results based on the first test
* Trying to generate VCF files and test out a couple more samples

* Completely disregard sites which are heterozygous in *any* sample?
* Ignore phasing and work off of base frequency per site?
2. Kmer table per strain per allele, based on reads mapping to the site (`barcodesearch_v2.py`)
  * Identify a unique Kmer within the barcoding locus for each strain/allele
  * Kmers extracted from reads, so *would* be phased

* Coverage - normalise by the count of a common between-strains Kmer?



Will PCR introduce bias when counting?
* Replicates to ensure we get the same results each time?


Need to generate test data to verify that the pipeline works







~

# Problem

1. Looking at the Live2Tell data, the barcoding function has failed to return the desired results

2. When viewed in IGV, many of the mapped reads are highlighted in red, meaning that the insert size is larger than expected
  * Is this a quirk of the view?
  * If the insert size is restricted, will this result in certain variants not being found?

3. When assembling regions identified as unique between strains in Bamboozle, the regions have turned out to be identical in some strains
  * 112 and 113 in particular seem to be very similar

  * Use bowtie2's --maxins/--minins flags?
    * Will certain alleles be unreported in this instance?
  * Use bcftools consensus's -I flag to report variants using ambiguity codes?
    * Levenshtein distance couldn't be used in hom vs. het situations
      * Max/min number of differences could be done later if needed

  * How does a reference-based assembly look in this case? Would a longer repeat region be assembled as intended?


# Solutions?

1. Enforce a fragment size limit on the mapping (e.g. +/- 10% of expected?)
  * Add a --fragment flag, which influences the Bowtie2 -I and -X flags

2. Apply a depth filter to ensure that DP in the VCF file is sufficiently high
  * This in conjunction with identifying average coverage post-mapping?

3. Does the quality threshold need adjusting?


# Proposed new workflow

1. Map the reads to the reference, but limit the fragment size to avoid repeat region issues
* If mapping has already been done, check the BAM header and flash a warning if the details don't add up

Map only concordant pairs, i.e. no singletons or discordant pairs?

2. Get an average (median?) coverage for the genome, for use with the DP field later

3. Run variant calling, but include a depth filter based on the genome coverage calculated above,  
   in addition to the existing quality filter (does this need adjusting?)

4. Assemble both alleles when calculating the consensus distances
