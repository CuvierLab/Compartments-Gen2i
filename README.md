# GENERAL INFO

All data analyzes have been performed with  R version 4.4.0 and tested also on 4.0.5.

As for libraries version with R 4.4.0 : 

- crayon : 1.5.2
- GenomicRanges : 1.56.0 
- dplyr : 1.1.4 
- ggplot2 : 3.5.1 
- data.table : 1.15.4 
- ggpubr : 0.6.0 
- karyoploteR : 1.30.0 
- rtracklayer : 1.64.0 
- corrr : 0.4.4 
- GGally : 2.2.1 
- ComplexHeatmap : 2.20.0 
- plyr : 1.8.9 
- CAinterprTools : 1.1.0 
- RColorBrewer : 1.1.3 
- ggVennDiagram : 1.5.2 
- ggdendro: 0.2.0 

and for R 4.0.5 :

- crayon : 1.5.2 
- GenomicRanges : 1.42.0 
- dplyr : 1.1.4 
- ggplot2 : 3.4.2 
- data.table : 1.15.2 
- ggpubr : 0.6.0 
- karyoploteR : 1.16.0 
- rtracklayer : 1.50.0 
- corrr :  0.4.4 
- GGally : 2.2.1 
- ComplexHeatmap : 2.20.0 
- plyr : 1.8.9 
- CAinterprTools : 1.1.0 
- RColorBrewer : 1.1.3 
- ggVennDiagram : 1.2.2 
- ggdendro: 0.2.0 


## PREPROCESSING
### Chromosomes arms

In order to generate chromsome arms bed positions use lift_over_chr_arms.R script.

### Genome bed and Tiled genome bed

Those data can be obtained through tiled_genome_bed.R script

### Union and intersection of peaksets

Union and intersection of peaksets can be retrieved using bedtools.sh script.

### TADscore BM to BigWig

To transform .bm tad separation score files to bigwig please use bm2bw.R.

### Contact matrix in .rds format 
Use g2ih5_To_cm_rds.R to turn g2i.h5 matrices into .rds R format


### ChIP-seq Profiles

After performing union & intersection of peakset, use computeMatrixProfile.sh script to compute ChIP-seq bigwig signal profile around peaks of interest.



## ANALYZES
This repo stores bioinformatic R scripts used to generate manuscript figures.
All scripts start with this line that has to be modified accordingly to  the path used to clone this repo : 

> setwd(dir = "config/src/R/github/")

## DATA PATH
### FOLDER STRUCTURE
Data must be stored in /data folder along with the scripts . It must contain the following folders :

- data/ChIPseq : Contains bed, bigwig and ${ip}_${condition}_newzs_prof.txt profile files.
- data/RNAseq : Contains deseq2_counts_all_conditions.txt and matrecap_genes_ce11.gtf files
- data/HiC : Contains all contact matrices (norm.KR.g2i.h5, obs_exp.g2i.h5, norm.KR.cm.rds,  obs_exp.cm.rds), all bigwig for Eigen vectors (pca1.bw), bigwig for TAD separation scores tad_score.bw and TAD in bed format

### NOTES
The following additional files must be stored in ChIPseq folder :

- ce11_noMT.bed
- ce11_tiled_10kb.bed
- ce11_tiled_1kb.bed 
- ce11_tiled_25kb.bed
- ce11_tiled_5kb.bed
- chr_arms_ce11.bed 
- genes_ce11.bed

If you prefer storing them elsewhere, you will need to adapt scripts for manuscript figures accordingly.
