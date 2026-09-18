# HEADER ====================================================== =
# GEN2I
# NAR revision 2026
# R version : R version 4.4.0 (2024-04-24), tested on 4.5.3
# System : x86_64, linux-gnu
# ============================================================= =
#
# Figure 3D - hierarchical clustering of the Hi-C contact matrices of the six
# genotypes, on the B compartment of the left chromosome arms.
#
# Each genotype contributes one vector of contact values: the KR-normalized
# intra-chromosomal contacts between the 25 kb bins that are both in the B
# compartment of wild type (PC1 of N2, initial batch, <= 0) and on a left
# chromosome arm. The three central diagonals are removed, because contacts at
# short genomic distance measure proximity along the polymer rather than
# compartmentalization. The six chromosomes are concatenated into a single
# vector per genotype, the vectors are z-scored, and the genotypes are clustered
# by Euclidean distance with complete linkage. Colours mark the four groups cut
# from the tree; the dashed grey/green branches are the ones above the cut.
#
# This is a clustering of genotypes, not of bins: the leaves are the six maps.
# ============================================================= =

# WORKING DIRECTORY ----
# Run the scripts from the root of this repository (where manuscript_lib.R sits):
# setwd("/path/to/Compartments-Gen2i")

# SOURCE ----
source("manuscript_lib.R")

# PARAMETERS ----
out_d       <- "output/figure_3D"
eigen_name  <- "eigen_pca1_N2.old"                           # PC1 of wild type, initial batch
eigen_bw    <- "data/HiC/N2-old_merged.bwa_mem.25kb.pca1.bw"
arm_bed     <- "data/ChIPseq/leftArm_ce11.bed"               # left arms, ce11 (see lift_over_chr_arms.R)
conditions  <- c("CEC4", "hpl2-lin61-old", "hpl2-old", "lin61-old",
                 "met2-set25-set32-old", "N2-old")           # I158A is excluded from the figures
h5_tpl      <- "data/HiC/%s_merged.bwa_mem.25kb.norm.KR.g2i.h5"
na_offset   <- 20L                                           # NA margin kept in the h5 for the APA
diag_offset <- 1L                                            # main diagonal + 1 off-diagonal on each side
nb_clust    <- 4L                                            # groups cut from the tree, for the colours
dir.create(out_d, showWarnings = FALSE, recursive = TRUE)

# GENOMIC OBJECT ----
## bins of 25 kb, index within chromosome, wild-type eigenvector, left arm flag
Go <- loadranges("data/ChIPseq/ce11_tiled_25kb.bed", genome = "ce11")
mcols(Go)["chr_bin"] <- unlist(lapply(split(Go, seqnames(Go)), function(chr) seq_along(chr)))
Go <- bw_tile(bw = eigen_bw, bed = Go, genome = "ce11", name = eigen_name)

arms_gr <- loadranges(arm_bed, genome = "ce11")
mcols(Go)["leftArm"] <- IRanges::overlapsAny(Go, arms_gr, ignore.strand = TRUE)

## the two filters of the panel: B compartment AND left arm
Go <- subset(Go, !is.na(get(eigen_name)) & get(eigen_name) <= 0 & leftArm)

Go_dt <- go2dt(Go)
chrs  <- unique(as.character(Go_dt$seqnames))

# RUN ----
cm_l <- lapply(setNames(conditions, conditions), function(cond) {
  h5 <- sprintf(h5_tpl, cond)
  cm_info <- rhdf5::h5ls(h5)
  sapply(cm_info$name, function(nm) HDF5Array::HDF5Array(h5, nm))
})

dendro_dt <- NULL
for (chr in chrs) {
  dendro_chr_dt <- NULL
  for (condition in names(cm_l)) {
    cm_c    <- as.matrix(as.matrix(cm_l[[condition]][[chr]]))
    ## the h5 carries an NA margin of 20 bins on each side, so that a 21x21 APA
    ## can be extracted at the very ends of a chromosome; drop it to get back to
    ## the bin indices of the tiled genome
    max_bin <- dim(cm_c)[1] - na_offset
    cm_c    <- cm_c[(na_offset + 1):max_bin, (na_offset + 1):max_bin]
    getKDiag(cm_c, -diag_offset:diag_offset) <- NA
    colRow_sel    <- Go_dt[seqnames == chr, chr_bin]
    dendro_chr_dt <- cbind(dendro_chr_dt,
                           data.table::setnames(
                             data.table(as.vector(cm_c[colRow_sel, colRow_sel])), condition))
  }
  dendro_dt <- rbind(dendro_dt, dendro_chr_dt)
}

## source data: the pairwise correlation of the six contact vectors, and the
## merge order and heights of the tree drawn in the panel
cor_df <- cor(dendro_dt, use = "na.or.complete")
fwrite(data.table(condition = rownames(cor_df), as.data.table(cor_df)),
       file.path(out_d, "correlation_matrix.csv"))

hc <- t(dendro_dt) %>% scale %>% dist(method = "euclidean") %>% hclust
fwrite(data.table(step = seq_len(nrow(hc$merge)), left = hc$merge[, 1],
                  right = hc$merge[, 2], height = hc$height),
       file.path(out_d, "hclust_merge.csv"))
fwrite(data.table(order = seq_along(hc$order), label = hc$labels[hc$order]),
       file.path(out_d, "hclust_leaves.csv"))

gg <- list()
hcdata <- dendro_data_k(hc, nb_clust)
cols   <- RColorBrewer::brewer.pal(nb_clust + 1, "Dark2")
gg[["Dendrogram.All_chrom.Dendogram"]] <- plot_ggdendro(hcdata,
                                                        direction   = "lr",
                                                        scale.color = cols,
                                                        expand.y    = 0.2) + theme_void()

# SAVE ----
for (nm in names(gg))
  ggsave(file.path(out_d, paste0(nm, ".png")), gg[[nm]], width = 8, height = 8, dpi = 150)
