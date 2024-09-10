# HEADER ====================================================== =
# GEN2I
# Tue Jun 18 18:22:27 2024
# R version : R version 4.4.0 (2024-04-24)
# System : x86_64, linux-gnu
# ============================================================= =


# WORKING DIRECTORY ----
setwd(dir = "config/src/R/github/")

# SOURCE ----
source("manuscript_lib.R")


# First download chr arms from ce4 annotation : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3046480/ 
# slide 7 of .xls for arms coordinates 

# Download chains
# wget -qO- https://hgdownload.soe.ucsc.edu/goldenPath/ce4/liftOver/ce4ToCe10.over.chain.gz | gunzip -c > ce4ToCe10.over.chain
# wget -qO- https://hgdownload.soe.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz | gunzip -c > ce10ToCe11.over.chain

chr_arms <- read.table("data/ext_data/chr_arms_ce4.gff",sep="\t",h=T)
colnames(chr_arms) <- c("seqnames","chr_anno","start","end")
chr_arms <- GRanges(chr_arms)
seqlevelsStyle(chr_arms) <- "Ensembl"
chr_arms_ce10 <- liftover(resize(chr_arms,1,"end"),chain = import.chain("data/ext_data/ce4ToCe10.over.chain"))
chr_arms_ce11 <- liftover(chr_arms_ce10,import.chain("data/ext_data/ce10ToCe11.over.chain"))
chr_arms_ce11_gap <- gaps(chr_arms_ce11)
chr_arms_ce11_gap$chr_anno <- chr_arms_ce11$chr_anno

