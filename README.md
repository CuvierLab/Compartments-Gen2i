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

## SEQ DATA ALIGNMENTS

### INIT
Download required files with src/00_init_project.sh
Then use related .sh script to align raw data for each sequencing type


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
Run the scripts from the root of this repository, where `manuscript_lib.R` sits. Each
script starts with a commented `setwd()` line, to be adapted or ignored depending on how
you launch it.

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

## FIGURES OF THE REVISED VERSION (NAR, 2026)

The scripts are named after the panels of the revised manuscript. Each one is standalone:
it rebuilds the genomic object from the tiled genome and the eigenvector bigwigs, sources
`manuscript_lib.R` for the shared functions, and writes its panels and its source tables to
`output/<panel>/`.

| Panel | Script | Analysis |
| --- | --- | --- |
| Fig. 2B, 2C | `figure_2B_2C_S4B_S4C.R` | Meta-profiles of PC1 across B-to-A compartment transitions, and their amplitudes |
| Fig. 3B | `figure_3B.R` | Differential saddle plots, bins ranked by the wild-type PC1 |
| Fig. 3C | `figure_3C.R` | Corner quantification of the Fig. 3B matrices, ratios B-B/A-A and B-B/B-A |
| Fig. 3E | `figure_3E.R` | Cumulative plots and within-map contact preferences, PC1 |
| Supplementary Fig. S4B, S4C | `figure_2B_2C_S4B_S4C.R` | As Fig. 2B and 2C, for PC2 |
| Supplementary Fig. S5B | `figure_S5B.R` | As Fig. 3B, ranked by the wild-type PC2 |
| Supplementary Fig. S5C | `figure_S5C.R` | As Fig. 3C, on the PC2-ranked matrices |
| Supplementary Fig. S5D | `figure_S5D.R` | As Fig. 3E, with compartments called by PC2 |
| Fig. 1B | `figure_1B.R` | Wild-type saddle plot, and density and intensity of the chromatin marks along the PC1 quantiles |
| Fig. 3F | `figure_3F.R` | Ratio of RNA-seq counts between compartments A and B, per genotype |
| Fig. 4A | `figure_4A.R` | Density of the chromatin marks along the quantiles of the mutant-minus-wild-type PC1 |
| Supplementary Fig. S4D | `figure_S4D.R` | PC1 tracks of chromosome I, mutant over wild type |
| Supplementary Fig. S6B | `figure_S6B.R` | TAD separation score along chromosome I |
| Supplementary Fig. S8A | `figure_S8A.R` | As Fig. 4A, for the intensity of the marks |
| Supplementary Fig. S8B | `figure_S8B.R` | Fisher tests between our peak sets and the published ones |
| Supplementary Fig. S8C | `figure_S8C.R` | Venn diagram of the HPL-2, LIN-61, H3K9me2 and H3K9me3 peak sets |
| Supplementary Fig. S8E | `figure_S8E.R` | Aggregated ChIP-seq z-scores in compartment B, mark against HPL-2 |

Every table produced by these scripts was checked against the pipeline that generated the
published figures: the values are identical.

The Hi-C panels use the initial batch (`-old` conditions) and its wild-type eigenvectors, at
25 kb. Bootstrap seeds are fixed in the scripts, so the published values are reproduced
exactly: 2000 resamplings of the corner pixels for the saddle ratios (Fig. 3C, S5C), 2000
resamplings of the transitions for the meta-profiles (Fig. 2B, 2C, S4B, S4C), and an
m-out-of-n bootstrap of 1000 contacts per interaction class and genotype over 2000
iterations for the contact preferences and the group contrast (Fig. 3E, S5D).

The scripts of the remaining panels (Fig. 3D, 4C, 4D, 5C and Supplementary Fig. S2A-C,
S4F, S7D-H, S8D) are being ported from our analysis pipeline and will be added in a further
version of this deposit.

Compartment B is defined, here as in the pipeline, as quantile groups 1 to 3 of the 50
quantiles of the wild-type PC1 of the initial batch (243 tiles of 25 kb); compartment A is
groups 16 to 50. Fig. S8E uses that definition.

The panels are saved as produced by R. The published figures are the same panels relabelled
and laid out for print: genotype names in italics, marks written H3K9me2 rather than
`H3K9me2_GEO_N2_max_signal`, and dendrograms moved. No value differs.

### ADDITIONAL DATA

Besides the files listed above, these scripts read from `data/`:

- `data/ChIPseq/ce11_tiled_25kb.bed`;
- `data/HiC/<condition>_merged.bwa_mem.25kb.pca1.bw` and `...pca2.bw`, the eigenvectors of
  each genotype of the initial batch;
- `data/HiC/<condition>_25kb.obs_exp.cm.rds`, one observed/expected matrix per genotype, for
  the cumulative plots;
- `data/HiC/<condition>_vs_N2-old.bwa_mem.25kb.norm.KR.g2i.h5`, the mutant over wild-type
  ratio matrices produced by `hicCompare.sh`, for the saddle plots;
- `data/HiC/N2-old_merged.bwa_mem.25kb.norm.KR.g2i.h5` for Fig. 1B;
- `data/HiC/<condition>_merged.bwa_mem.25kb_tad_score.bw` for Supplementary Fig. S6B;
- the ChIP-seq peak sets and z-score bigwigs listed above, and
  `data/RNAseq/deseq2_counts_all_conditions.txt` for Fig. 3F.
