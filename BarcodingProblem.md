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
