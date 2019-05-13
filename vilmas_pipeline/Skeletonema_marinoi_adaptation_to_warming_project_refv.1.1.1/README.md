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
|**Genotype and phenotype files**|G:`/proj/data11/vilma/Bamboozle/vilmas_pipeline/dev/Skeletonema_marinoi_adaptation_to_warming_project_marinoi_refv.1.1.1/pheno_geno_correlation/*_corr_test/genotype_vcf012_*_wo_extraNcol.txt` <br/> P:`/proj/data11/vilma/Bamboozle/vilmas_pipeline/dev/Skeletonema_marinoi_adaptation_to_warming_project_marinoi_refv.1.1.1/pheno_geno_correlation/*_corr_test/phenotypes_*.txt`| 

---
### Raw data
|GROUP/POPULATION|SAMPLE|
|---|---|
|**control_1960**| `P8352_136` |
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

The samples were divided in four groups/populations, control_1960, control_2000, cold_1960 and warm_2000 with 1, 4, 6 and 12 individuals. 
Control_1960 was excluded from the analysis due to its extreme low sample size (one individual). 
The control_2000 individuals were from the control, referred to as the cold sediment core in the project description. 
The cold_1960 and the warm_2000 were taken from same core, referred to as the warm sediment core.  

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

`fasta2tab UniqueModels.pep.fasta | grep -f annotations_flt0.2_uniq_mRNA-1_genenames.txt | tab2fasta > temp_pep.fasta`  

The `temp_pep.fasta` is the input file to the `eggNOG` script to get the functional annotation of the protein sequences. 

sge script `eggnog.sge`:  

```
#$ -cwd
#$ -q Annotation-3
#$ -pe mpi 10
#$ -S /bin/bash

#module load eggnog-mapper/v1.0.3 

python2 /usr/local/packages/eggnog-mapper-1.0.3/emapper.py -i temp_pep.fasta --output exon_w900_s450_diamond -m diamond --cpu 10

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


### pheno geno correlation  

> `vcftools --012`  

Preparing files for correlation test.  

```
1) Merged all exon parsed files from Bcftools directory using bcftools merge. 
Resulting in three merged files: cold_merged_exon.vcf, control_merged_exon.vcf and warm_merged_exon.vcf.
2) Merge the three files to one large merged file: cold_control_warm_merged_exon.vcf.
3) Preparation of genotypes file for correlation test
```

```  
Run vcftools --012

$ vcftools --012 —vcf cold_control_warm_merged_exon.vcf --out [OUTFILE]
OUTPUT:
0 May  9 10:50 genotypes_vcftools_012.txt
154M May  9 10:54 genotypes_vcftools_012.txt.012
704 May  9 10:53 genotypes_vcftools_012.txt.012.indv
49M May  9 10:54 genotypes_vcftools_012.txt.012.pos
903 May  9 10:54 genotypes_vcftools_012.txt.log

# Transpose using awk

awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str”\t”a[i,j];
        }
        print str
    }
}' FILE

# Add SNPs as row names

$ paste genotypes_vcftools_012.txt.012.pos transpose_genotypes_vcftools_012.txt.012 > OUTFILE

# Use vim to add column names (Sample names)

```  

```
4) Filtering the file so I have one file with only 0,1,2 values (remove -1 by using vim: ":g/-1/d"). 
This is used as genotype input file to the correlation test. 
With the resulting format:
```

```
#CHROM	POS	P8352_145	P8352_146	P8352_147	P8352_150
Sm_000111F-001-01	4912	2	2	2	2
Sm_000111F-001-01	4926	2	2	2	2
Sm_000111F-001-01	5306	2	2	1	2
Sm_000111F-001-01	27495	2	2	2	2
Sm_000111F-001-01	34864	2	2	2	2
Sm_000111F-001-01	35209	2	1	2	2
Sm_000111F-001-01	41771	2	2	2	2
Sm_000111F-001-01	45270	2	2	2	2
Sm_000111F-001-01	52438	2	2	2	2
```

```
5) The phenotype file (from Kai Lohbeck, University of Gothenburg) was filtered and had the resulting format:
```
```
Strain	treatment	Location	Depth(cm)	mean_max	mean_mode	mean_shape
P8352_145	temp	1	21	0,799147762	13,84614141	-4,07870399
P8352_146	temp	1	21	0,804494333	14,63967677	-0,04974438
P8352_147	temp	1	21	0,748955673	20,77559596	-0,633426031
P8352_150	temp	1	21	0,627595509	6,768969697	-6,355021381
```  

```
6) R scripts were provided by Marina Axelsson-Fisk (Chalmers) and made to fit the warming project data instead. 
```  

```
The scripts (total 9):
/proj/data11/vilma/Bamboozle/vilmas_pipeline/dev/Skeletonema_marinoi_adaptation_to_warming_project_refv.1.1.1/pheno_geno_correlation/*_corr_test/s_marinoi_mean_*.R
```  

```
7) I used the script annotations_gff.py and a GFF file with all features to get the gene names as preparation for EggNOG annotation as in previous steps (`/proj/data11/vilma/Sm_OnemRNAPerGene.FinalWithSequence.ManualCuration_nofasta_sorted_flt.gff`). 
```

Results:  

**warm_2000**  

```
0.997984002015998	Sm_000005F_1341045
0.997953002046998	Sm_000005F_1341115
0.997914002085998	Sm_000027F_373084
0.998040001959998	Sm_000027F_373089
0.998023001976998	Sm_000027F_373144
0.997969002030998	Sm_000027F_373145
0.997972002027998	Sm_000027F_373262
0.997946002053998	Sm_000027F_373264
0.997969002030998	Sm_000027F_373281
0.997924002075998	Sm_000027F_373284
0.998021001978998	Sm_000000F_2343744
0.997906002093998	Sm_000009F_440051
0.995473004526995	Sm_000041F_65148
0.995415004584995	Sm_000070F_231318
0.995430004569995	Sm_000070F_231353
0.995515004484996	Sm_000070F_231387
0.995260004739995	Sm_000000F_2168161
0.995553004446996	Sm_000000F_2223953
0.995380004619995	Sm_000000F_2223956
0.995391004608995	Sm_000000F_2223976
mean_max adjusted	
	
9.99999000001e-07	Sm_000002F_270669
0.998737001262999	Sm_000006F_859671
0.998728001271999	Sm_000001F_658247
0.997962002037998	Sm_000004F_585489
0.997918002081998	Sm_000004F_585758
0.998016001983998	Sm_000006F_859703
0.997981002018998	Sm_000006F_859977
0.997945002054998	Sm_000000F_1701870
9.99999000001e-07	Sm_000014F_709027
9.99999000001e-07	Sm_000028F_90113
9.99999000001e-07	Sm_000028F_90165
9.99999000001e-07	Sm_000028F_90167
9.99999000001e-07	Sm_000028F_90389
9.99999000001e-07	Sm_000000F_2258098
9.99999000001e-07	Sm_000002F_1170684
0.00110599889400111	Sm_000011F_587732
0.995914004085996	Sm_000034F_621441
0.995986004013996	Sm_000000F_222525
0.995925004074996	Sm_000038F_18565
0.002000997999002	Sm_000004F-009-01_22264
mean_mode adjusted	
	
0.985012014987985	Sm_000000F_1735602
0.984846015153985	Sm_000000F_2265593
0.984816015183985	Sm_000013F_558414
0.984919015080985	Sm_000021F_426236
0.984944015055985	Sm_000021F_426340
0.984882015117985	Sm_000021F_426364
0.984826015173985	Sm_000021F_426430
0.984795015204985	Sm_000021F_426442
0.985160014839985	Sm_000021F_426529
0.984900015099985	Sm_000021F_426542
0.985016014983985	Sm_000021F_426642
0.985038014961985	Sm_000021F_426650
0.984838015161985	Sm_000021F_426700
0.984637015362985	Sm_000021F_426753
0.984954015045985	Sm_000021F_426955
0.984988015011985	Sm_000021F_427003
0.984758015241985	Sm_000021F_427015
0.984998015001985	Sm_000021F_427077
0.985038014961985	Sm_000021F_427181
0.984633015366985	Sm_000021F_427239
mean_shape adjusted	
``` 

```
# emapper version: emapper-a3830d3 emapper DB: 4.5.1
# command: ./emapper.py  -i pep_annotations_warm_2000_corr_genenames.fasta --output warm_2000_correlation_test_diamond -m diamond --cpu 10
# time: Mon May 13 10:31:15 2019
#query_name	seed_ortholog_evalue	predicted_gene_name	eggNOG annot
Sm_t00000038-RA	2.4e-95		protein kinase kinase kinase
Sm_t00000570-RA	1.1e-37		alpha-galactosidase
Sm_t00000732-RA	1.9e-198	UBA1	Ubiquitin-like modifier activating enzyme
Sm_t00000766-RA	1.2e-117	DNAJC2	Transcription factor
Sm_t00004202-RA	6.2e-120		
Sm_t00004561-RA	9.2e-51		Syntaxin 6, N-terminal
Sm_t00002637-RA	1e-62		Ankyrin Repeat
Sm_t00016573-RA	3.5e-81	PSME4	Proteasome (Prosome, macropain) activator subunit 4
Sm_t00006438-RA	1.1e-63	TRMU	tRNA 5-methylaminomethyl-2-thiouridylate methyltransferase
Sm_t00005056-RA	2.8e-09	GNT1	Glycosyltransferase (GlcNAc)
# 10 queries scanned
# Total time (seconds): 45.1475229263
# Rate: 0.22 q/s
```  

**cold_1960**  

```
0.833448166551833	Sm_000012F-008-01_155064
0.833004166995833	Sm_000012F-008-01_155424
0.832804167195833	Sm_000012F-008-01_155450
0.832924167075833	Sm_000012F-008-01_155934
0.833077166922833	Sm_000012F-008-01_156000
0.833377166622833	Sm_000012F-008-01_156183
9.99999000001e-07	Sm_000032F-002-01_112369
9.99999000001e-07	Sm_000032F-002-01_112370
9.99999000001e-07	Sm_000032F-002-01_112412
0.833138166861833	Sm_000074F_169267
0.833044166955833	Sm_000074F_169330
0.833072166927833	Sm_000074F_169392
0.833459166540833	Sm_000074F_169414
0.833685166314834	Sm_000074F_169417
0.832971167028833	Sm_000074F_169420
0.833479166520833	Sm_000074F_169522
0.833419166580833	Sm_000074F_173596
0.832953167046833	Sm_000074F_173649
0.833617166382834	Sm_000010F_478555
0.832768167231833	Sm_000010F_480458
mean_max adjusted	
	
0.833769166230834	Sm_000007F-006-01_9130
0.833705166294834	Sm_000025F-004-01_80648
0.833248166751833	Sm_000074F_184302
9.99999000001e-07	Sm_000125R_87901
9.99999000001e-07	Sm_000125R_87903
9.99999000001e-07	Sm_000010F_718918
0.833845166154834	Sm_000010F_817478
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
	
0.833756166243834	Sm_000007F-006-01_9130
0.832779167220833	Sm_000025F-004-01_80648
0.833672166327834	Sm_000074F_184302
9.99999000001e-07	Sm_000125R_87901
9.99999000001e-07	Sm_000125R_87903
9.99999000001e-07	Sm_000010F_718918
0.833227166772833	Sm_000010F_817478
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

```
# emapper version: emapper-a3830d3 emapper DB: 4.5.1
# command: ./emapper.py  -i pep_annotations_cold_1960_corr_genenames.fasta --output cold_1960_correlation_test_diamond -m diamond --cpu 10
# time: Mon May 13 11:25:12 2019
#query_name	seed_ortholog_evalue	predicted_gene_name	eggNOG annot
Sm_t00018801-RA	8.2e-10	PDE2	Guanylate cyclase
Sm_t00018796-RA	1.3e-83	PIGM	Mannosyltransferase
Sm_t00018798-RA	6.4e-09	OCAR_7315	Sel1 domain protein repeat-containing protein
Sm_t00020094-RA	1.7e-19		
Sm_t00006919-RA	5.4e-112		6-phosphofructokinase
Sm_t00006920-RA	9.3e-33		Sodium:solute symporter family
Sm_t00007013-RA	2.7e-102		DUF1295 domain protein
Sm_t00007045-RA	7.3e-14		WW domain
Sm_t00007047-RA	4.3e-16		Tesmin TSO1-like CXC domain containing protein
# 9 queries scanned
# Total time (seconds): 47.9950580597
# Rate: 0.19 q/s
```


**control_2000**  

```
0.933694066305934	Sm_000032F-002-01_112412
0.933205066794933	Sm_000032F-002-01_113777
0.933151066848933	Sm_000032F-002-01_113869
0.933082066917933	Sm_000009F-003-01_42490
0.933620066379934	Sm_000074F_190950
0.933038066961933	Sm_000074F_191578
0.933701066298934	Sm_000001F-008-01_18879
0.932940067059933	Sm_000000F-015-01_698
0.933530066469934	Sm_000010F_548252
0.933540066459934	Sm_000010F_869296
0.933419066580933	Sm_000010F_903693
0.932977067022933	Sm_000010F_903879
0.933313066686933	Sm_000010F_907178
0.933313066686933	Sm_000010F_907182
0.933213066786933	Sm_000010F_907183
0.932964067035933	Sm_000010F_907184
0.933379066620933	Sm_000010F_907188
0.933583066416934	Sm_000010F_913633
0.933255066744933	Sm_000010F_914084
0.933553066446934	Sm_000010F_914302
	
mean_max adjusted	
	
0.933130066869933	Sm_000074F_187741
0.933913066086934	Sm_000010F_863717
0.933551066448934	Sm_000010F_872308
0.933456066543933	Sm_000010F_912741
0.933286066713933	Sm_000010F_912751
0.933011066988933	Sm_000010F_912913
0.933456066543933	Sm_000010F_912921
0.932805067194933	Sm_000010F_912924
0.933353066646933	Sm_000010F_912930
0.933274066725933	Sm_000010F_912948
0.933266066733933	Sm_000010F_913112
0.933461066538933	Sm_000010F_1173785
0.932977067022933	Sm_000054F_146542
0.933078066921933	Sm_000054F_172249
0.933178066821933	Sm_000054F_301651
0.933373066626933	Sm_000111F_9007
0.933269066730933	Sm_000111F_9008
0.933629066370934	Sm_000111F_9017
0.933555066444934	Sm_000111F_9164
0.933637066362934	Sm_000111F_20977
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
9.99999000001e-07	Sm_000010F_1173785
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

```
# emapper version: emapper-a3830d3 emapper DB: 4.5.1
# command: ./emapper.py  -i pep_annotations_control_2000_corr_genenames.fasta --output control_2000_correlation_test_diamond -m diamond --cpu 10
# time: Mon May 13 11:34:45 2019
#query_name	seed_ortholog_evalue	predicted_gene_name	eggNOG annot
Sm_t00018398-RA	2.3e-95		lysine (K)-specific demethylase
Sm_t00018802-RA	1.8e-301		FAD binding domain
Sm_t00018803-RA	3.3e-12		
Sm_t00018502-RA	8.5e-20	UBI4	Ubiquitin exists either covalently attached to another protein, or free (unanchored). When covalently bound, it is conjugated to target proteins via an isopeptide bond either as a monomer (monoubiquitin), a polymer linked via different Lys residues of the ubiquitin (polyubiquitin chains) or a linear polymer linked via the initiator Met of the ubiquitin (linear polyubiquitin chains). Polyubiquitin chains, when attached to a target protein, have different functions depending on the Lys residue of the ubiquitin that is linked Lys-6-linked may be involved in DNA repair
Sm_t00006945-RA	4.2e-08		Protein of unknown function (DUF1336)
Sm_t00007065-RA	2.8e-30	ACTR1A	Actins are highly conserved proteins that are involved in various types of cell motility and are ubiquitously expressed in all eukaryotic cells
Sm_t00007071-RA	1.1e-41	CCNF	cyclin a
Sm_t00007077-RA	2.3e-146	ABCF1	ATP-binding cassette, subfamily F (GCN20), member 1
Sm_t00007083-RA	1.7e-28	BUD32	Component of the EKC KEOPS complex which promotes both telomere uncapping and telomere elongation (By similarity). The complex is required for efficient recruitment of transcriptional coactivators. Important for bud site selection (By similarity)
Sm_t00007082-RA	2.4e-60	JMJD8	jumonji domain containing
# 10 queries scanned
# Total time (seconds): 61.1684319973
# Rate: 0.16 q/s
```


**The warm results:**  

|warm_2000|Permutation result|Contig|Position|GFF (feature and gene name)|eggNOG annotation|non_synonymous|   
|---|---|---|---|---|---|---|   
|**mean_max**|0.997984002015998 |`Sm_000005F`| 1341045 | `CDS` **Sm_t00002637-RA** |`Ankyrin Repeat`| yes |
|**mean_max**|0.997953002046998| `Sm_000005F` | 1341115| `CDS` **Sm_t00002637-RA** | `Ankyrin Repeat` | no| 
|**mean_max**| 0.997906002093998|`Sm_000009F`|440051| `CDS` **Sm_t00006438-RA**|**TRMU** `tRNA 5-methylaminomethyl-2-thiouridylate methyltransferase`| no | 
|**mean_max** | 0.995260004739995 |`Sm_000000F`|2168161| `CDS` **Sm_t00000732-RA** | **UBA1** `Ubiquitin-like modifier activating enzyme` | no | 
|**mean_max** |0.995553004446996|`Sm_000000F`| 2223953 | `3' UTR` **Sm_t00000766-RA** | **DNAJC2**  `Transcription factor`| no | 
|**mean_max** |0.995380004619995 | `Sm_000000F`|  2223956| `3' UTR` **Sm_t00000766-RA**| **DNAJC2**  `Transcription factor` |yes|  
|**mean_mode** |0.998737001262999|`Sm_000006F`|859671|`CDS` **Sm_t00005056-RA**|`GNT1 Glycosyltransferase (GlcNAc)`|no|
|**mean_mode** |0.998016001983998|`Sm_000006F`|859703|`CDS` **Sm_t00005056-RA**|`GNT1 Glycosyltransferase (GlcNAc)`|yes|
|**mean_mode** |0.997981002018998|`Sm_000006F`|859977|`mRNA` **Sm_t00005056-RA**|`GNT1 Glycosyltransferase (GlcNAc)`|no|
|**mean_mode** |0.995986004013996|`Sm_000000F`|222525|`CDS` **Sm_t00000038-RA**|`protein kinase kinase kinase`|no|

Comparing the phenotype genotype correlation test results with the `eq_system.py` FST results showed that 11 of 14 gene names were found in the `eq_system.py` results, however the FST values were of varying results:

`grep -f annotations_warm_2000_corr_uniq_genenames.txt temp_annotation_genenames.txt`  

```
Sm_t00000038-RA
Sm_t00000552-RA
Sm_t00000570-RA
Sm_t00000732-RA
Sm_t00000808-RA
Sm_t00004202-RA
Sm_t00004561-RA
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
The sample sizes were altogether too low for a fair assessment of the environmental impact on the genotypes. 
Nevertheless the analysis provided an example of how much easier the process become when using a customized pipeline in bioinformatics. The tests were performed as if the sample sizes and data were sufficient and the process, programs and models have shown to be successful. 
