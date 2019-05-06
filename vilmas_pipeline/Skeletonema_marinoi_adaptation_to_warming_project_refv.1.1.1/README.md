# *Skeletonema marinoi* adaptation to warming analysis  
> Using Bamboozle   
> For more information about the options and manuals: https://github.com/topel-research-group/Bamboozle/wiki    

**This is a log of the analysis, with comments and the commands used to get the results.**  

---

### Files used in analysis: 

|FILE|PATH|
|---|---|
|**FASTA reference**|`/proj/data11/Skeletonema_marinoi_Ref_v1.1.1.fst` |
|**FASTQ files**|`/proj/data21/Skeletonema_marinoi/Genome/Skeletonema_marinoi_adaptation_to_warming_project/00_data/A.Godhe_17_01-P8352/`| 
|**Pipeline output files**|`/proj/data11/vilma/Pipeline_vilma/dev/Skeletonema_marinoi_adaptation_to_warming_project_marinoi_refv.1.1.1/`|  
|**GFF file for GFF parsing and phenotype genotype correlation test**| **Original:** <br/> `/proj/data26/Skeletonema_marinoi_genome_project/03_Annotation/Skeletonema_marinoi_Ref_v1.1_Primary/Unique_models_per_locus_ManualCuration/Sm_OnemRNAPerGene.FinalWithSequence.ManualCuration.gff` <br/> **Filtered:** <br/> `/proj/data11/vilma/Sm_OnemRNAPerGene.FinalWithSequence.ManualCuration_nofasta_sorted_flt.gff`|  
|**GFF file for filtering out the annotations**|`/proj/data11/vilma/Skeletonema_marinoi_Ref_v1.1_Primary.all.gff` <br/> `/proj/data11/vilma/CDS.gff`|  
|**Map file for gene names**|`/proj/data26/Skeletonema_marinoi_genome_project/03_Annotation/Skeletonema_marinoi_Ref_v1.1_Primary/Uniq_models_per_loci/Sm_Uniq_models_per_loci.ID.map`|    
|**FASTA file for protein sequences**|`/proj/data26/Skeletonema_marinoi_genome_project/03_Annotation/Skeletonema_marinoi_Ref_v1.1_Primary/Unique_models_per_locus_ManualCuration/FASTA_files/UniqueModels.pep.fasta`|  

---
### Raw data
|GROUP/POPULATION|SAMPLE|
|---|---|
|**control_2000**| `P8352_145` <br/> `P8352_146` <br/> `P8352_147` </br> `P8352_150`|
|**cold_1960**|`P8352_101` <br/> `P8352_103` <br/> `P8352_108` <br/> `P8352_109` <br/> `P8352_111` <br/> `P8352_112`|
|**warm_2000**|`P8352_113` <br/> `P8352_114` <br/> `P8352_115` <br/> `P8352_116` <br/> `P8352_117` <br/> `P8352_118` <br/> `P8352_119` <br/> `P8352_120` <br/> `P8352_121` <br/> `P8352_122` <br/> `P8352_123` <br/> `P8352_124` <br/>|
---

The sge script that was used for the pipeline is marinoi.sge. With the options:  
`vilmas_pipeline.py --threads $NSLOTS --ref Skeletonema_marinoi_Ref_v1.1.1.fst --forward *R1.Pair.fastq.gz --reverse *R2.Pair.fastq.gz --clean --done`  

The analysis was first made in `/proj/data11/vilma/Pipeline_vilma/dev/Skeletonema_marinoi_adaptation_to_warming_project_marinoi_refv.1.1.1/` before the merge with Bamboozle.
So I rsync the results from running the pipeline to cwd.  

```
rsync -hav /proj/data11/vilma/Pipeline_vilma/dev/Skeletonema_marinoi_adaptation_to_warming_project_marinoi_refv.1.1.1/control_2000/ \
/proj/data11/vilma/Bamboozle/vilmas_pipeline/dev/Skeletonema_marinoi_adaptation_to_warming_project_marinoi_refv.1.1.1/control_2000 

rsync -hav /proj/data11/vilma/Pipeline_vilma/dev/Skeletonema_marinoi_adaptation_to_warming_project_marinoi_refv.1.1.1/cold_1960/ \
/proj/data11/vilma/Bamboozle/vilmas_pipeline/dev/Skeletonema_marinoi_adaptation_to_warming_project_marinoi_refv.1.1.1/cold_1960  

rsync -hav /proj/data11/vilma/Pipeline_vilma/dev/Skeletonema_marinoi_adaptation_to_warming_project_marinoi_refv.1.1.1/warm_2000/ \
/proj/data11/vilma/Bamboozle/vilmas_pipeline/dev/Skeletonema_marinoi_adaptation_to_warming_project_marinoi_refv.1.1.1/warm_2000  
```

I checked that the directories were the same using diff (there's probably a better way, md5sum?)     
`diff -qr DIR1 DIR2`    

### v.1.1.1_fststatistics   

> `angsd_fst.py`    

Then the angsd analysis using `angsd_fst.py`:  
The original GFF file was from:    
`/proj/data26/Skeletonema_marinoi_genome_project/03_Annotation/Skeletonema_marinoi_Ref_v1.1_Primary/Unique_models_per_locus_ManualCuration/Sm_OnemRNAPerGene.FinalWithSequence.ManualCuration.gff`   

```
---Same for all---
#$ -cwd
#$ -q Annotation-3 
#$ -pe mpi 40 
#$ -S /bin/bash

module load angsd/v0.923 

---Difference---
/proj/data11/vilma/Bamboozle/vilmas_pipeline/scripts/angsd_fst.py \
-1 ../../control_2000/ \
-2 ../../cold_1960/ \
--ref /proj/data11/Skeletonema_marinoi_Ref_v1.1.1.fst \
--gff ../../../tmp_pipeline/Sm_OnemRNAPerGene.FinalWithSequence.ManualCuration_nofasta_sorted_flt.gff \
--feature exon \
--contigsizes ../../../tmp_pipeline/chrom.genome \
--clean -w 900 -s 450 

/proj/data11/vilma/Bamboozle/vilmas_pipeline/scripts/angsd_fst.py \
-1 ../../control_2000/ \
-2 ../../warm_2000/ \
--ref /proj/data11/Skeletonema_marinoi_Ref_v1.1.1.fst \
--gff ../../../tmp_pipeline/Sm_OnemRNAPerGene.FinalWithSequence.ManualCuration_nofasta_sorted_flt.gff \
--feature exon \
--contigsizes ../../../tmp_pipeline/chrom.genome \
--clean -w 900 -s 450

/proj/data11/vilma/Bamboozle/vilmas_pipeline/scripts/angsd_fst.py \
-1 ../../warm_2000/ \
-2 ../../cold_1960/ \
--ref /proj/data11/Skeletonema_marinoi_Ref_v1.1.1.fst \
--gff ../../../tmp_pipeline/Sm_OnemRNAPerGene.FinalWithSequence.ManualCuration_nofasta_sorted_flt.gff \
--feature exon \
--contigsizes ../../../tmp_pipeline/chrom.genome \
--clean -w 900 -s 450

```  

### eq_system  

> `eq_system.py`  

Then `eq_system.py` analysis on resulting ANGSD files.  

```
../../../scripts/eq_system.py \
-1 ../v.1.1.1_fststatistics/warm_2000_cold_1960/ANGSD/angsd_fst_headers_exon.csv \
-2 ../v.1.1.1_fststatistics/control_2000_cold_1960/ANGSD/angsd_fst_headers_exon.csv \
-3 ../v.1.1.1_fststatistics/control_2000_warm_2000/ANGSD/angsd_fst_headers_exon.csv \
--mean --angsd
```  

Printing: 
 
```  
Fontconfig warning: ignoring UTF-8: not a valid region tag
--- 12.287084341049194 seconds ---  
```   

Take `temp.csv` and filter out the outliers to a new file to use as search pattern for annotation script.  
First sort on generic numerical value and then filtering.  

`cat temp.csv | sort -t ',' -k3g | tr ',' '\t' | awk '{if ($3 >= 0.2) print}' | tr '\t' ',' > temp_exon_sort_flt0.2.csv`  

The contigs that appeared more often than others in FST values over 0.2 were found using the command:

`cut -d ',' -f1 temp_exon_sort_flt0.2.csv | sort |uniq -c | sort -n`

The contigs that appeared over 100 times were:

```
101 Sm_000016F
102 Sm_000006F
102 Sm_000008F
118 Sm_000010F
139 Sm_000001F
139 Sm_000004F
139 Sm_000014F
140 Sm_000003F
140 Sm_000012F
227 Sm_000000F
```

### eggNOG  

> `annotations_gff.py`   
> eggNOG -> `emapper.py`   

Take filtered `temp_exon_sort_flt0.2.csv` file and use as input in annotation script and the GFF file as second input file.   
First remove header, then sort and convert csv to table. 
  
`echo "$(tail -n +2 ../eq_system/temp_exon_sort_flt0.2.csv)" > temp_exon_sort_flt0.2_noheader.csv`

`cat temp_exon_sort_flt0.2_noheader.csv | tr ',' '\t' | sort -k2n,2 > temp_exon_sort_flt0.2_noheader.table`  

Use `annotations_gff.py` to get the annotations from the GFF file.  
The GFF is a filtered file with only CDS results from 
`/proj/data11/vilma/Skeletonema_marinoi_Ref_v1.1_Primary.all.gff`.

```
Inputs: 
exon = open('temp_exon_sort_flt0.2_noheader.table', 'r').read().splitlines()
gff = open('/proj/data11/vilma/CDS.gff', 'r').read().splitlines()
```  

`./annotations_gff.py > annotations_flt0.2.txt`  

Select only the uniq, grep on "mRNA-1" and filter out only the names from the annotation file.

`grep "mRNA-1" annotations_flt0.2.txt | sort -u > annotations_flt0.2_uniq_mRNA-1.txt`  

Open vim and use:

``` 
1) %s/^ID=[0-9]\{0,3}://gc 
2) %s/^ID=//gc

```

`cut -d ":" -f1 annotations_flt0.2_uniq_mRNA-1.txt > annotations_flt0.2_uniq_mRNA-1_just_names.txt`

Should look something like this:  

```
maker-Sm_000018F-snap-gene-9.271-mRNA-1
maker-Sm_000030F-snap-gene-2.161-mRNA-1
maker-Sm_000030F-snap-gene-4.134-mRNA-1
maker-Sm_000084F-snap-gene-1.174-mRNA-1
```

Then I use grep to filter out the only matching annotations in the file `Sm_Uniq_models_per_loci.ID.map` 
from `/proj/data26/Skeletonema_marinoi_genome_project/03_Annotation/Skeletonema_marinoi_Ref_v1.1_Primary/Uniq_models_per_loci/Sm_Uniq_models_per_loci.ID.map`:  

`grep -f annotations_flt0.2_uniq_mRNA-1_just_names.txt Sm_Uniq_models_per_loci.ID.map | cut -f2 > annotations_flt0.2_uniq_mRNA-1_genenames.txt`

Looks something like this:

```
Sm_t00000006-RA
Sm_t00000007-RA
Sm_t00000009-RA
Sm_t00000022-RA
```
Comparing the `annotations_flt0.2_uniq_mRNA-1_just_names.txt` and `annotations_flt0.2_uniq_mRNA-1_genenames.txt`
shows that the first have 2081 lines and the second have 2068 lines. The different lines were found using:

Where `field1_genenames.txt` is made using the same command as above but cutting field 1:

`grep -f annotations_flt0.2_uniq_mRNA-1_just_names.txt Sm_Uniq_models_per_loci.ID.map | cut -f1 > field1_genenames.txt`

`grep -v -f field1_genenames.txt annotations_flt0.2_uniq_mRNA-1_just_names.txt`

This command print following differences:

```
maker-Sm_000018F-snap-gene-9.271-mRNA-1
maker-Sm_000030F-snap-gene-2.161-mRNA-1
maker-Sm_000030F-snap-gene-4.134-mRNA-1
maker-Sm_000154F-snap-gene-0.112-mRNA-1
maker-Sm_000186F-snap-gene-0.53-mRNA-1
snap_masked-Sm_000021F-processed-gene-1.89-mRNA-1
snap_masked-Sm_000025F-processed-gene-2.81-mRNA-1
snap_masked-Sm_000025F-processed-gene-3.105-mRNA-1
snap_masked-Sm_000069F-processed-gene-0.64-mRNA-1
snap_masked-Sm_000079F-processed-gene-0.98-mRNA-1
snap_masked-Sm_000114F-processed-gene-0.132-mRNA-1
snap_masked-Sm_000125R-processed-gene-0.109-mRNA-1
snap_masked-Sm_000188F-processed-gene-0.49-mRNA-1
```

These names were not found in the mapping file `Sm_Uniq_models_per_loci.ID.map`. 
When the regions were viewed in the Skeletonema marinoi genome browser genes were 
found and when searching the gene name found in the genome browser it was present in the `Sm_Uniq_models_per_loci.ID.map` file. 
The problem seems to be that the annotation names are not linked to the gene names in some cases. 

Then use `fasta2tab` and `tab2fasta` and `grep -f` to get the protein sequences, 
with the FASTA file `UniqueModels.pep.fasta` from 
`/proj/data26/Skeletonema_marinoi_genome_project/03_Annotation/Skeletonema_marinoi_Ref_v1.1_Primary/Unique_models_per_locus_ManualCuration/FASTA_files/UniqueModels.pep.fasta`.  

`fasta2tab UniqueModels.pep.fasta | grep -f annotations_flt0.2_uniq_mRNA-1_genenames.txt | tab2fasta > pep.fasta`  

The `pep.fasta` is the input file to the `eggNOG` script to get the functional annotation of the protein sequences. 

sge script `eggnog.sge`:  

```
#$ -cwd
#$ -q Annotation-3
#$ -pe mpi 10
#$ -S /bin/bash

#module load eggnog-mapper/v1.0.3 

python2 /usr/local/packages/eggnog-mapper-1.0.3/emapper.py -i pep.fasta --output exon_w900_s450_diamond -m diamond --cpu 10

#emapper.py -i pep.fasta --output exon_diamond -m diamond --cpu 10
```  

Filter out the interesting information of the resulting file, contigname, scores, name of gene and annotation of the function.

`cat exon_w900_s450_diamond.emapper.annotations | cut -f1,3,5,13 | awk -F '\t' '$3!=""' > genes.txt`   
`wc -l: 548 genes.txt`  

Looks something like this:   

```
#query_name	seed_ortholog_evalue	predicted_gene_name	eggNOG annot
Sm_t00000006-RA	7.2e-36	YBEQ	Sel1 domain protein repeat-containing protein
Sm_t00000084-RA	5.4e-143	ALLC	Allantoicase
Sm_t00000098-RA	3.2e-20	NAS2	26S Proteasome non-ATPase regulatory subunit
Sm_t00000174-RA	4.9e-108	OSI_18976	CCT motif family protein
```


**NON_SYNONYMOUS MUTATIONS**:  

Merging vcf files from `snpeff_exon.vcf.gz` using `bcftools merge` then filter out only the non synonymous annotated SNPs.  

From sge script:  
`bcftools merge --threads $NSLOTS *merged_exon.vcf.gz -Oz -o cold_control_warm_merged_exon.vcf.gz` 

`grep "NON_SYNONYMOUS" cold_control_warm_merged_exon.vcf.gz > cold_control_warm_merged_exon_non_synonymous.vcf`

(Old version of SnpEff, using the flag `-classic` for EFF annotations. Newer version show `missense_variant` instead of `NON_SYNONYMOUS`)

`wc -l: 436411 cold_control_warm_merged_exon_non_synonymous.vcf`

Use `SnpSift` to filter out the annotation names.

`java -jar /usr/local/packages/snpEff/SnpSift.jar extractFields -e "." -s "," cold_control_warm_merged_exon_non_synonymous.vcf "EFF[*].TRID" > non_syn.snpsift.table`
`wc -l: 436412 non_syn.snpsift.table (with header from snpsift)`

Filter the file so you only have the unique names, same format as previously shown.  

`cat non_syn.snpsift.table | cut -d ':' -f2 | cut -d ',' -f1 | sort -u > just_names_non_syn_uniq.txt` 
`wc -l: 21077 just_names_non_uniq.txt`

Filter out the interesting non-synonymous annotation names by using the file with annotated region names previously obtained from `ANGSD`.

`grep -f annotations_flt0.2_uniq_mRNA-1_just_names.txt just_names_non_syn_uniq.txt > non_synonymous_just_names.txt`  
`wc -l: 1966 non_synonymous_just_names.txt` 

Using the unique non-synonymous annotation names as search pattern against the mapping file.  

`grep -f non_synonymous_just_names.txt Sm_Uniq_models_per_loci.ID.map |cut -f2 > non_synonymous_genenames.txt`  
`wc -l: 1963 non_synonymous_genenames.txt`

When comparing `non_synonymous_just_names.txt` and `non_synonymous_genenames.txt` the first file had 1966 lines and the second has 1963 lines. 
The lines that differed were found using `grep -f`.

Where `field1_genenames_non_syn.txt` was made using the same command as above but cutting field 1 instead:

`grep -f non_synonymous_just_names.txt Sm_Uniq_models_per_loci.ID.map |cut -f1 > field1_genenames_non_syn.txt`  

`grep -v -f field1_genenames_non_syn.txt non_synonymous_just_names.txt`

The command printed following differences:

```
maker-Sm_000018F-snap-gene-9.271-mRNA-1
maker-Sm_000186F-snap-gene-0.53-mRNA-1
snap_masked-Sm_000021F-processed-gene-1.89-mRNA-1
```

These names were not found in the mapping file `Sm_Uniq_models_per_loci.ID.map`. 
When the regions were viewed in the Skeletonema marinoi genome browser the same problem as before occurred. 
There are genes in the regions but there seem to be a problem in the linkage between the annotation name and the gene name for some reason.

Then use `fasta2tab` and `tab2fasta` and `grep -f` to get the protein sequences. 

`fasta2tab UniqueModels.pep.fasta | grep -f non_synonymous_genenames.txt | tab2fasta > pep_non_synonymous.fasta` 

Then using `eggNOG` to get the functional annotation as previous.

Filter out the interesting information of the resulting file, contigname, scores, name of gene and annotation of the function.

`cat non_synonymous_diamond.emapper.annotations |cut -f1,3,5,13 | awk -F '\t' '$3!=""' > genes_non_synonymous.txt`
`wc -l: 511 genes_non_synonymous.txt`

### pheno geno correlation  

> `getgenosfromvcf_clean.py`  

Preparing files for correlation test.  

```
1) Merged all exon parsed files from Bcftools directory using bcftools merge. 
Resulting in three merged files: cold_merged_exon.vcf, control_merged_exon.vcf and warm_merged_exon.vcf.
2) Merge the three files to one large merged file: cold_control_warm_merged_exon.vcf.
3) Used the script getgenosfromvcf_clean.py (Pierre) to get matrix.
```

`python2 ./getgenosfromvcf_clean.py cold_control_warm_merged_exon.vcf genos_exon.vcf rows 1` 

Because my vcf files don't contain the GQ (GATK) I altered the script, 
setting every GQ to 1, in this case it means assuming that every site provided in the vcf files were of acceptable quality, 
just to make it work (need to find a solution for this in the future, 
maybe using the PL score and taking the second smallest value and using that instead).

```
4) Filtering the resulting file replacing 00 -> 0, 01 -> 1, 11 -> 2 and everything else to -1. 
```

``` 
# In vim:
:%s/\t00/\t0/gc
:%s/\t01/\t1/gc
:%s/\t11/\t2/gc
:%s/\.\./-1/gc
:%s/\t[0-9]\{2}\t/\t-1\t/gc
:g/-1/d
```

```
5) Filtering the file so I have one file with only 0,1,2 values. 
This is used as genotype input file to the correlation test. 
genos_exon_flt.txt
```
Unfiltered result file:  
lines = 2715372  

Filtered result file:  
lines = 55839

```
6) The phenotype file (from Kai Lohbeck) was filtered and had the resulting format:
```
```
Strain	treatment	Location	g1_AlignmenttoSm_Refv111	mean_max	mean_mode	mean_shape
P8352_145	temp	1	77,21	0,799147762	13,84614141	-4,07870399
P8352_146	temp	1	81,02	0,804494333	14,63967677	-0,04974438
P8352_147	temp	1	74,78	0,748955673	20,77559596	-0,633426031
P8352_150	temp	1	73,86	0,627595509	6,768969697	-6,355021381
```  

```
7) R scripts were provided by Marina Axelsson-Fisk and made to fit the warming project data instead. 
```  

```
The scripts (total 9):
/proj/data11/vilma/Bamboozle/vilmas_pipeline/dev/Skeletonema_marinoi_adaptation_to_warming_project_refv.1.1.1/pheno_geno_correlation/*_corr_test/s_marinoi_mean_*.R
```  

I used a GFF file with all features to get the gene names as preparation for EggNOG annotation as in previous steps (`Sm_OnemRNAPerGene.FinalWithSequence.ManualCuration_nofasta_sorted_flt.gff`). 

Results:  

**warm_2000**  

```
0.997921002078998	Sm_000005F_1341045
0.997898002101998	Sm_000005F_1341115
0.997977002022998	Sm_000027F_373084
0.997992002007998	Sm_000027F_373089
0.997942002057998	Sm_000027F_373144
0.997987002012998	Sm_000027F_373145
0.997983002016998	Sm_000027F_373262
0.997954002045998	Sm_000027F_373264
0.998004001995998	Sm_000027F_373281
0.998028001971998	Sm_000027F_373284
0.998041001958998	Sm_000000F_2343744
0.997932002067998	Sm_000009F_440051
0.995454004545995	Sm_000041F_65148
0.995451004548995	Sm_000070F_231318
0.995408004591995	Sm_000070F_231353
0.995478004521996	Sm_000070F_231387
0.995435004564995	Sm_000000F_2168161
0.995413004586995	Sm_000000F_2223953
0.995341004658995	Sm_000000F_2223956
0.995498004501996	Sm_000000F_2223976
mean_max adjusted

9.99999000001e-07	Sm_000002F_270669
0.998772001227999	Sm_000006F_859671
0.998694001305999	Sm_000001F_658247
0.997906002093998	Sm_000004F_585489
0.997921002078998	Sm_000004F_585758
0.998034001965998	Sm_000006F_859703
0.997951002048998	Sm_000006F_859977
0.998005001994998	Sm_000000F_1701870
9.99999000001e-07	Sm_000014F_709027
9.99999000001e-07	Sm_000028F_90113
9.99999000001e-07	Sm_000028F_90165
9.99999000001e-07	Sm_000028F_90167
9.99999000001e-07	Sm_000028F_90389
9.99999000001e-07	Sm_000000F_2258098
9.99999000001e-07	Sm_000002F_270684
0.00113499886500114	Sm_000011F_587732
0.995933004066996	Sm_000034F_621441
0.996069003930996	Sm_000000F_222525
0.996067003932996	Sm_000038F_18565
0.00202899797100203	Sm_000004F-009-01_22264
mean_mode adjusted

0.984713015286985	Sm_000000F_1735602
0.985100014899985	Sm_000000F_2265593
0.984653015346985	Sm_000013F_558414
0.984808015191985	Sm_000021F_426236
0.984811015188985	Sm_000021F_426340
0.984705015294985	Sm_000021F_426364
0.985012014987985	Sm_000021F_426430
0.984783015216985	Sm_000021F_426442
0.984862015137985	Sm_000021F_426529
0.984742015257985	Sm_000021F_426542
0.984890015109985	Sm_000021F_426642
0.984765015234985	Sm_000021F_426650
0.984877015122985	Sm_000021F_426700
0.985028014971985	Sm_000021F_426753
0.984669015330985	Sm_000021F_426955
0.984990015009985	Sm_000021F_427003
0.984940015059985	Sm_000021F_427015
0.985084014915985	Sm_000021F_427077
0.985115014884985	Sm_000021F_427181
0.984718015281985	Sm_000021F_427239
mean_shape adjusted
```
**eggNOG warm**

```
# emapper version: emapper-541224e emapper DB: 4.5.1
# command: ./emapper.py  -i warm_2000_corr.fasta --output warm_2000_correlation_test_diamond -m diamond --cpu 10
# time: Fri Apr 26 10:07:37 2019
#query_name     seed_ortholog_evalue    predicted_gene_name     eggNOG annot
Sm_t00000038-RA 2.4e-95         protein kinase kinase kinase
Sm_t00000570-RA 1.1e-37         alpha-galactosidase
Sm_t00000732-RA 1.9e-198        UBA1    Ubiquitin-like modifier activating enzyme
Sm_t00000766-RA 1.2e-117        DNAJC2  Transcription factor
Sm_t00004202-RA 6.2e-120                
Sm_t00002637-RA 1e-62           Ankyrin Repeat
Sm_t00016573-RA 3.5e-81 PSME4   Proteasome (Prosome, macropain) activator subunit 4
Sm_t00006438-RA 1.1e-63 TRMU    tRNA 5-methylaminomethyl-2-thiouridylate methyltransferase
Sm_t00005056-RA 2.8e-09 GNT1    Glycosyltransferase (GlcNAc)
# 9 queries scanned
# Total time (seconds): 58.9582929611
# Rate: 0.15 q/s
```

**cold_1960**  

```
0.833658166341834	Sm_000012F-008-01_155064
0.833545166454834	Sm_000012F-008-01_155424
0.832851167148833	Sm_000012F-008-01_155450
0.833773166226834	Sm_000012F-008-01_155934
0.833219166780833	Sm_000012F-008-01_156000
0.832871167128833	Sm_000012F-008-01_156183
9.99999000001e-07	Sm_000032F-002-01_22369
9.99999000001e-07	Sm_000032F-002-01_22370
9.99999000001e-07	Sm_000032F-002-01_22412
0.832804167195833	Sm_000074F_169107
0.833338166661833	Sm_000074F_169267
0.833608166391834	Sm_000074F_169330
0.833073166926833	Sm_000074F_169392
0.833149166850833	Sm_000074F_169414
0.833843166156834	Sm_000074F_169417
0.832576167423833	Sm_000074F_169420
0.833302166697833	Sm_000074F_169522
0.832753167246833	Sm_000074F_173596
0.832515167484833	Sm_000074F_173649
0.832856167143833	Sm_000010F_478555
mean_max adjusted

0.833617166382834	Sm_000007F-006-01_9130
0.833376166623833	Sm_000025F-004-01_80648
0.833386166613833	Sm_000074F_184302
9.99999000001e-07	Sm_000125R_87901
9.99999000001e-07	Sm_000125R_87903
9.99999000001e-07	Sm_000010F_718918
0.833031166968833	Sm_000010F_817478
9.99999000001e-07	Sm_000010F_817527
9.99999000001e-07	Sm_000010F_817776
9.99999000001e-07	Sm_000010F_824065
9.99999000001e-07	Sm_000010F_824692
9.99999000001e-07	Sm_000010F_824696
9.99999000001e-07	Sm_000010F_824704
9.99999000001e-07	Sm_000010F_824705
9.99999000001e-07	Sm_000010F_824745
9.99999000001e-07	Sm_000010F_826293
9.99999000001e-07	Sm_000010F_826651
9.99999000001e-07	Sm_000010F_826754
9.99999000001e-07	Sm_000010F_827430
9.99999000001e-07	Sm_000010F_827601
mean_mode adjusted

0.833446166553833	Sm_000007F-006-01_9130
0.833730166269834	Sm_000025F-004-01_80648
0.832566167433833	Sm_000074F_184302
9.99999000001e-07	Sm_000125R_87901
9.99999000001e-07	Sm_000125R_87903
9.99999000001e-07	Sm_000010F_718918
0.833157166842833	Sm_000010F_817478
9.99999000001e-07	Sm_000010F_817527
9.99999000001e-07	Sm_000010F_817776
9.99999000001e-07	Sm_000010F_824065
9.99999000001e-07	Sm_000010F_824692
9.99999000001e-07	Sm_000010F_824696
9.99999000001e-07	Sm_000010F_824704
9.99999000001e-07	Sm_000010F_824705
9.99999000001e-07	Sm_000010F_824745
9.99999000001e-07	Sm_000010F_826293
9.99999000001e-07	Sm_000010F_826651
9.99999000001e-07	Sm_000010F_826754
9.99999000001e-07	Sm_000010F_827430
9.99999000001e-07	Sm_000010F_827601
mean_shape adjusted
```   
**eggNOG cold** 

```
# emapper version: emapper-541224e emapper DB: 4.5.1
# command: ./emapper.py  -i cold_1960_corr.fasta --output cold_1960_correlation_test_diamond -m diamond --cpu 10
# time: Fri Apr 26 10:20:15 2019
#query_name     seed_ortholog_evalue    predicted_gene_name     eggNOG annot
Sm_t00018801-RA 8.2e-10 PDE2    Guanylate cyclase
Sm_t00018796-RA 1.3e-83 PIGM    Mannosyltransferase
Sm_t00018798-RA 6.4e-09 OCAR_7315       Sel1 domain protein repeat-containing protein
Sm_t00020094-RA 1.7e-19         
Sm_t00006919-RA 5.4e-112                6-phosphofructokinase
Sm_t00007013-RA 2.7e-102                DUF1295 domain protein
Sm_t00007045-RA 7.3e-14         WW domain
Sm_t00007047-RA 4.3e-16         Tesmin TSO1-like CXC domain containing protein
# 8 queries scanned
# Total time (seconds): 49.2979159355
# Rate: 0.16 q/s
```


**control_2000**  

```
0.933159066840933	Sm_000032F-002-01_22412
0.933264066735933	Sm_000032F-002-01_23777
0.933181066818933	Sm_000032F-002-01_23869
0.933217066782933	Sm_000009F-003-01_42490
0.932974067025933	Sm_000074F_190950
0.933666066333934	Sm_000074F_191578
0.933408066591933	Sm_000001F-008-01_18879
0.933117066882933	Sm_000000F-015-01_698
0.932943067056933	Sm_000010F_869296
0.932473067526932	Sm_000010F_903693
0.933204066795933	Sm_000010F_903879
0.933422066577933	Sm_000010F_907178
0.933484066515934	Sm_000010F_907182
0.933176066823933	Sm_000010F_907183
0.933561066438934	Sm_000010F_907184
0.933632066367934	Sm_000010F_907188
0.933139066860933	Sm_000010F_913633
0.933121066878933	Sm_000010F_914084
0.933157066842933	Sm_000010F_914302
0.933234066765933	Sm_000010F_997492
mean_max adjusted

0.933163066836933	Sm_000074F_187741
0.933376066623933	Sm_000010F_863717
0.933767066232934	Sm_000010F_872308
0.933224066775933	Sm_000010F_912741
0.933137066862933	Sm_000010F_912751
0.933027066972933	Sm_000010F_912913
0.933584066415934	Sm_000010F_912921
0.933514066485933	Sm_000010F_912924
0.933900066099934	Sm_000010F_912930
0.933873066126934	Sm_000010F_912948
0.933663066336934	Sm_000010F_913112
0.933620066379934	Sm_000010F_273785
0.933446066553933	Sm_000054F_146542
0.933258066741933	Sm_000054F_172249
0.933618066381934	Sm_000054F_301651
0.933288066711933	Sm_000111F_9007
0.933370066629933	Sm_000111F_9008
0.932989067010933	Sm_000111F_9017
0.933489066510933	Sm_000111F_9164
0.932987067012933	Sm_000111F_20977
mean_mode adjusted

9.99999000001e-07	Sm_000074F_187741
9.99999000001e-07	Sm_000010F_863717
9.99999000001e-07	Sm_000010F_872308
9.99999000001e-07	Sm_000010F_912741
9.99999000001e-07	Sm_000010F_912751
9.99999000001e-07	Sm_000010F_912913
9.99999000001e-07	Sm_000010F_912921
9.99999000001e-07	Sm_000010F_912924
9.99999000001e-07	Sm_000010F_912930
9.99999000001e-07	Sm_000010F_912948
9.99999000001e-07	Sm_000010F_913112
9.99999000001e-07	Sm_000010F_273785
9.99999000001e-07	Sm_000054F_146542
9.99999000001e-07	Sm_000054F_172249
9.99999000001e-07	Sm_000054F_301651
9.99999000001e-07	Sm_000111F_9007
9.99999000001e-07	Sm_000111F_9008
9.99999000001e-07	Sm_000111F_9017
9.99999000001e-07	Sm_000111F_9164
9.99999000001e-07	Sm_000111F_20977
mean_shape adjusted
```
**eggNOG control**

```
# emapper version: emapper-541224e emapper DB: 4.5.1
# command: ./emapper.py  -i control_2000_corr.fasta --output control_2000_correlation_test_diamond -m diamond --cpu 10
# time: Fri Apr 26 10:15:09 2019
#query_name     seed_ortholog_evalue    predicted_gene_name     eggNOG annot
Sm_t00018398-RA 2.3e-95         lysine (K)-specific demethylase
Sm_t00018802-RA 1.8e-301                FAD binding domain
Sm_t00018803-RA 3.3e-12         
Sm_t00018502-RA 8.5e-20 UBI4    Ubiquitin exists either covalently attached to another protein, or free (unanchored). When covalently bound, it is conjugated to target proteins via an isopeptide bond either as a monomer (monoubiquitin), a polymer linked via different Lys residues of the ubiquitin (polyubiquitin chains) or a linear polymer linked via the initiator Met of the ubiquitin (linear polyubiquitin chains). Polyubiquitin chains, when attached to a target protein, have different functions depending on the Lys residue of the ubiquitin that is linked Lys-6-linked may be involved in DNA repair
Sm_t00007065-RA 2.8e-30 ACTR1A  Actins are highly conserved proteins that are involved in various types of cell motility and are ubiquitously expressed in all eukaryotic cells
Sm_t00007071-RA 1.1e-41 CCNF    cyclin a
Sm_t00007077-RA 2.3e-146        ABCF1   ATP-binding cassette, subfamily F (GCN20), member 1
Sm_t00007083-RA 1.7e-28 BUD32   Component of the EKC KEOPS complex which promotes both telomere uncapping and telomere elongation (By similarity). The complex is required for efficient recruitment of transcriptional coactivators. Important for bud site selection (By similarity)
Sm_t00007082-RA 2.4e-60 JMJD8   jumonji domain containing
Sm_t00007106-RA 1.5e-72 ANAPC7  anaphase promoting complex subunit 7
# 10 queries scanned
# Total time (seconds): 39.3249869347
# Rate: 0.25 q/s
```

**The warm results:**  

|warm_2000|Permutation result|Contig|Position|GFF (feature and gene name)|eggNOG annotation|non_synonymous|   
|---|---|---|---|---|---|---|   
|**mean_max**|0.997921002078998 |`Sm_000005F`| 1341045 | `CDS` **Sm_t00002637-RA** |`Ankyrin Repeat`| yes |
|**mean_max**|0.997898002101998| `Sm_000005F` | 1341115| `CDS` **Sm_t00002637-RA** | `Ankyrin Repeat` | yes| 
|**mean_max**| 0.997932002067998|`Sm_000009F`|440051| `CDS` **Sm_t00006438-RA**|**TRMU** `tRNA 5-methylaminomethyl-2-thiouridylate methyltransferase`| yes | 
|**mean_max** | 0.995435004564995 |`Sm_000000F`|2168161| `CDS` **Sm_t00000732-RA** | **UBA1** `Ubiquitin-like modifier activating enzyme` | yes | 
|**mean_max** |0.995413004586995|`Sm_000000F`| 2223953 | `3' UTR` **Sm_t00000766-RA** | **DNAJC2**  `Transcription factor`| yes | 
|**mean_max** |0.995341004658995 | `Sm_000000F`|  2223956| `3' UTR` **Sm_t00000766-RA**| **DNAJC2**  `Transcription factor` |yes|  
|**mean_mode** |0.996069003930996|`Sm_000000F`|222525|`CDS` **Sm_t00000038-RA**|`protein kinase kinase kinase`|yes|

Comparing the phenotype genotype correlation test results with the `eq_system.py` FST results showed that 10 of 14 gene names were found in the `eq_system.py` results, however the FST values were of varying results:

`grep -f annotations_warm_2000_corr_uniq_genenames.txt ../../eggNOG/Temperature_exon/temp_annotation_genenames.txt`  

```
Sm_t00000038-RA
Sm_t00000552-RA
Sm_t00000570-RA
Sm_t00000732-RA
Sm_t00000808-RA
Sm_t00004202-RA
Sm_t00005056-RA
Sm_t00006438-RA
Sm_t00010656-RA
Sm_t00016573-RA
```
  					
`cat annotations_warm_2000_corr_uniq_genenames.txt`   	
				
```
Sm_t00000038-RA
Sm_t00000552-RA
Sm_t00000570-RA
Sm_t00000732-RA
Sm_t00000766-RA
Sm_t00000808-RA
Sm_t00002637-RA
Sm_t00004202-RA
Sm_t00005056-RA
Sm_t00006007-RA
Sm_t00006438-RA
Sm_t00008013-RA
Sm_t00010656-RA
Sm_t00016573-RA
```  


## Conclusion
The sample sizes were altogether too low for a fair assessment of the environmental impact on the genotypes. Nevertheless the analysis provided an example of how much easier the process become when using a customized pipeline in bioinformatics. The pipeline makes the analysis easy to reproduce and the work that is put in the process was only at one time, hence the human error was lowered as well. The tests were performed as if the sample sizes and data were sufficient and the process, programs and models have shown to be successful. The phenotype genotype correlation test uses an updated GFF file to get the gene names, this has not get been done on the FST results (only an older version) this might give some other results in the annotation. I decided to trust the results from the correlation test more and have therefore put the FST analysis a bit on hold. 


