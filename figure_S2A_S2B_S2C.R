# HEADER ====================================================== =
# GEN2I
# NAR revision 2026
# R version : R version 4.4.0 (2024-04-24), tested on 4.5.3
# System : x86_64, linux-gnu
# ============================================================= =
#
# Supplementary Fig. S2A, S2B and S2C - wild-type saddle plots.
#
# Same analysis as Fig. 1B (see figure_1B.R): the 25 kb bins of the initial
# batch wild-type map are ranked by an eigenvector and grouped into 50 quantile
# groups per chromosome; each pixel of the saddle matrix is the mean of the KR
# normalised contacts of the corresponding pair of quantile groups. B-B contacts
# sit in the top left corner, A-A contacts in the bottom right.
#
# The three panels differ only in the ranking eigenvector and in the scope:
#   S2A : PC1 of the initial-batch wild type, chromosome I alone;
#   S2B : the same, chromosome III alone;
#   S2C : PC2 of the merged-batch wild type, all chromosomes merged.
#
# S2C is the pipeline behaviour, kept as is: the ranker is the PC2 of the merged
# batch (eigen_pca2_N2), applied to the contact matrix of the initial batch.
#
# The per-chromosome panels are read off the same run as the merged one: the
# quantile groups are always computed chromosome by chromosome, so a panel of
# chromosome I is the chromosome I sub-matrix, not a separate ranking.
#
# The bars on the two margins give the median eigenvector value of each quantile
# group, over the chromosome shown.
# ============================================================= =

# WORKING DIRECTORY ----
# Run the scripts from the root of this repository (where manuscript_lib.R sits):
# setwd("/path/to/Compartments-Gen2i")

# SOURCE ----
source("manuscript_lib.R")

# PARAMETERS ----
out_d   <- "output/figure_S2A_S2B_S2C"
tiles_p <- "data/ChIPseq/ce11_tiled_25kb.bed"
h5_p    <- "data/HiC/N2-old_merged.bwa_mem.25kb.norm.KR.g2i.h5"  # wild type, initial batch
chr_v   <- c("I", "II", "III", "IV", "V", "X")
N_BINS  <- 50L                     # quantile groups per axis
AGR_SAD <- "mean"                  # how a saddle pixel is summarised
AGR_MCOL <- "median"               # how the margin bars are summarised
AGR_SCM <- "median"                # how the chromosome matrices are merged
NA_OFFSET   <- 20L                 # bins trimmed at both ends of each matrix
DIAG_OFFSET <- 1L                  # off-diagonals set to NA
dir.create(out_d, showWarnings = FALSE, recursive = TRUE)

# The two rankings used by the three panels.
rankers <- list(
  eigen_pca1_N2.old = "data/HiC/N2-old_merged.bwa_mem.25kb.pca1.bw",   # initial batch, PC1
  eigen_pca2_N2     = "data/HiC/N2_merged.bwa_mem.25kb.pca2.bw"        # merged batch, PC2
)

# panel -> ranking + scope ("mergedChr", or a chromosome name)
panels <- list(
  S2A = list(eigen = "eigen_pca1_N2.old", scope = "I"),
  S2B = list(eigen = "eigen_pca1_N2.old", scope = "III"),
  S2C = list(eigen = "eigen_pca2_N2",     scope = "mergedChr")
)

# RUN ----
## SADDLE MATRICES ----
### Contact matrices, one per chromosome ----
cm_info <- rhdf5::h5ls(h5_p)
cm_l    <- sapply(cm_info$name[cm_info$name %in% chr_v],
                  function(nm) HDF5Array::HDF5Array(h5_p, nm))

sadmat_l <- list()
for (eigen_name in names(rankers)) {

  ### MCOLS Chr bin ----
  Go <- loadranges(tiles_p, genome = "ce11")
  mcols(Go)[["chr_bin"]] <- unlist(lapply(split(Go, seqnames(Go)),
                                          function(chr) seq_along(chr)))

  ### MCOLS Eigen values ----
  Go <- bw_tile(bw = rankers[[eigen_name]], bed = Go,
                genome = config.genome, name = eigen_name)
  Go <- subset(Go, seqnames %in% chr_v)

  ### SADDLE ----
  sadmat_l[[eigen_name]] <- getSaddle2D(bin_go = Go,
                                        bait = eigen_name,
                                        anchor = eigen_name,
                                        bait_quant = N_BINS,
                                        anchor_quant = N_BINS,
                                        na_offset = NA_OFFSET,
                                        diag_offset = DIAG_OFFSET,
                                        cm_l = cm_l,
                                        f_score = AGR_MCOL,
                                        f_sad = AGR_SAD,
                                        f_cm = AGR_SCM,
                                        log2 = FALSE)
}

## PLOTS ----
gg <- list()
for (panel in names(panels)) {
  eigen_name <- panels[[panel]]$eigen
  scope      <- panels[[panel]]$scope
  sadmat     <- sadmat_l[[eigen_name]]

  if (scope == "mergedChr") {
    m_raw  <- sadmat[["sadmat_cm_mChr"]]
    stat_b <- sadmat[["sum_stat_mChr_b"]]
    stat_a <- sadmat[["sum_stat_mChr_a"]]
  } else {
    m_raw  <- sadmat[["sadmat_cm_aChr"]][[scope]]
    stat_b <- sadmat[["sum_stat_aChr_b"]][seqnames == scope, ]
    stat_a <- sadmat[["sum_stat_aChr_a"]][seqnames == scope, ]
  }

  # SOURCE TABLE ----
  fwrite(data.table(m_raw, keep.rownames = "bait_quantile"),
         file.path(out_d, sprintf("saddle_matrix.%s.csv", panel)))
  fwrite(merge(stat_b, stat_a, by.x = "quant_bait", by.y = "quant_anchor",
               all = TRUE, sort = TRUE),
         file.path(out_d, sprintf("quantile_medians.%s.csv", panel)))

  m <- matrix2tibble(m_raw, N_BINS, N_BINS)
  gg_center <- ggplot(m, aes(x = Anchor, y = Bait, fill = counts)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.y = element_blank(), axis.title.x = element_blank()) +
    scico::scale_fill_scico(palette = "lajolla", direction = -1) +
    scale_x_discrete(position = "top")

  gg_bait <- ggplot(stat_b, aes(x = factor(quant_bait), y = tot_bait)) +
    geom_bar(stat = "identity", aes(fill = tot_bait)) + theme_bw() +
    theme(legend.position = "none", axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90), axis.text.y = element_blank()) +
    scale_x_discrete(position = "bottom", limits = rev) +
    labs(y = "Tile median values", x = eigen_name) + coord_flip() +
    scico::scale_fill_scico(palette = "bamako", direction = 1)

  gg_anchor <- ggplot(stat_a, aes(x = factor(quant_anchor), y = tot_anchor)) +
    geom_bar(stat = "identity", aes(fill = tot_anchor)) + theme_bw() +
    theme(legend.position = "none", axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    labs(y = "Tile median values", x = eigen_name) +
    scale_x_discrete(position = "top") +
    scico::scale_fill_scico(palette = "bamako", direction = 1)

  tmp_p  <- cowplot::plot_grid(NULL, gg_anchor, gg_bait,
                               gg_center + theme(legend.position = "none"),
                               rel_widths = c(1, 4), rel_heights = c(2, 4),
                               ncol = 2, align = "hv")
  legend <- get_legend(gg_center + theme(legend.box.margin = margin(140, 0, 20)))
  nm <- sprintf("Saddle.%s.N2.%s", panel,
                if (scope == "mergedChr") "mergedChr" else paste0("chr", scope))
  gg[[nm]] <- cowplot::plot_grid(tmp_p, legend, rel_widths = c(5, 1))
}

# SAVE ----
for (nm in names(gg))
  ggsave(file.path(out_d, paste0(nm, ".png")), gg[[nm]],
         width = 10, height = 8, dpi = 150)
