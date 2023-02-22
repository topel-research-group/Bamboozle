#!/bin/bash

# Pipeline for quantifying amplicon-sequenced Bamboozle-derived barcodes

################################################################################
# 0) Let's set some variables!
################################################################################

# Absolute path to root directory of your input data
INDATA=""

# Location of your reference fasta
REFFILE=""


# Trimming settings

## 3' adapter for Cutadapt to trim off (Illumina universal adapter: "AGATCGGAAGAG")
ADAPTER=""

## Forward primer for Cutadapt to trim off
FWDPRIMER=""

## Reverse primer for Cutadapt to trim off
REVPRIMER=""

## Read quality to trim from 5' end
QUALITY=

## Minimum trimmed read length
MINLENGTH=


# Which program will you use for read merging?
# Either "bbmerge" or "dada2"
MERGE=""

################################################################################
# 1) Set up a data directory to perform quality control (and merging)
################################################################################

# Make a symbolic link to the reference FASTA for later
ln -s $REFFILE reference.fasta

# Make a directory for your data, and move into it
mkdir 00_data
cd 00_data

# Make a subdirectory for each sample, and create links to the raw data in the relevant subdirectory
# This saves space versus copying the data
while read i; do
        SAMPLE=$(echo $i | cut -f1)
        mkdir ${SAMPLE}

        FWDREADS=$(find "${INDATA}/${SAMPLE}" -type f -name "*R1*fastq.gz")
        REVREADS=$(find "${INDATA}/${SAMPLE}" -type f -name "*R2*fastq.gz")

        ln -s $FWDREADS ${SAMPLE}/
        ln -s $REVREADS ${SAMPLE}/
done < ../index_file.lst

################################################################################
# 2) Trim your input data
################################################################################

# Trim each pair of read files, and produce a FastQC summary file
# Note: the -a/-A flags refer to 3' adapters, the -g/-G flags refer to 5' primers,
# and -n 2 tells Cutadapt to remove two sequences (i.e. adapter AND primer)
# Parameters for quality (-q) and minimum length (-m) can be adjusted according to your requirements)
while read i; do
        FWDREADS=$(ls ${i}/${i}*R1*fastq.gz)
        REVREADS=$(ls ${i}/${i}*R2*fastq.gz)
        cutadapt -a $ADAPTER -A $ADAPTER -g $FWDPRIMER -G $REVPRIMER -n 2 -q $QUALITY -m $MINLENGTH \
                -o ${i}/${i}.R1.trimmed.fastq.gz -p ${i}/${i}.R2.trimmed.fastq.gz $FWDREADS $REVREADS
        fastqc -o ${i} ${i}/${i}.R1.trimmed.fastq.gz ${i}/${i}.R2.trimmed.fastq.gz
done < ../index_file.lst

multiqc -ip -o MultiQC_Report .

################################################################################
# 3) Optional pre-merging of reads
################################################################################

if [ $MERGE == "bbmerge" ]
then
	while read i; do
	        bbmerge.sh in1=${i}/${i}.R1.trimmed.fastq.gz in2=${i}/${i}.R2.trimmed.fastq.gz \
		        out=${i}/${i}.merged.fastq.gz \
			outu1=${i}/${i}.unmerged.R1.fastq.gz outu2=${i}/${i}.unmerged.R2.fastq.gz
	done < ../index_file.lst
else
	echo "Merging will be done in dada2; skipping BBmerge."
fi

################################################################################
# 4) Denoising and quantifying reads with dada2
################################################################################

cd ..

if [ $MERGE == "dada2" ]
then
	Rscript rundada2.PE.R
elif [ $MERGE == "bbmerge" ]
then
	Rscript rundada2.BBmerged.R
fi

################################################################################
# 6) Visualisation
################################################################################

mkdir 02_plots
cd 02_plots

# Plug in the R script here
Rscript ../Barcode_analysis.R
