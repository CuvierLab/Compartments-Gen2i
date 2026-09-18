# HEADER ====================================================== =
# GEN2I
# NAR revision 2026
# R version : R version 4.4.0 (2024-04-24), tested on 4.5.3
# System : x86_64, linux-gnu
# ============================================================= =
#
# Supplementary Figure S5B - differential saddle plots, one per mutant.
#
# Bins are ranked by the wild-type PC2 and grouped into 50 quantile groups. Each
# pixel of the saddle matrix is the mean of the mutant/wild-type contact ratios
# of the corresponding pair of quantile groups, and the genome-wide matrix is
# the pixel-wise median of the six chromosome matrices. B-B contacts sit in the
# top left corner, A-A contacts in the bottom right, B-A contacts in the top
# right.
#
# For display each matrix is rescaled to its own range ([-1, 1]): the panels
# show the pattern within a matrix, not the amplitude of the change, which is
# quantified in Supplementary Figure S5C (see figure_S5C.R). The bars on the two margins give
# the median eigenvector value of each quantile group.
#
# Figure 3B is the same analysis ranked by PC1: see figure_3B.R.
# ============================================================= =

# WORKING DIRECTORY ----
# Run the scripts from the root of this repository (where manuscript_lib.R sits):
# setwd("/path/to/Compartments-Gen2i")

# SOURCE ----
source("manuscript_lib.R")

# PARAMETERS ----
out_d      <- "output/figure_S5B"
eigen_name <- "eigen_pca2_N2.old"                            # PC2 of wild type, initial batch
eigen_bw   <- "data/HiC/N2-old_merged.bwa_mem.25kb.pca2.bw"
conditions <- c("CEC4", "hpl2-old", "lin61-old", "hpl2-lin61-old", "met2-set25-set32-old")
h5_tpl     <- "data/HiC/%s_vs_N2-old.bwa_mem.25kb.norm.KR.g2i.h5"   # ratio mutant / wild type
chr_v      <- c("I", "II", "III", "IV", "V", "X")
N_BINS     <- 50L                                             # quantile groups per axis
dir.create(out_d, showWarnings = FALSE, recursive = TRUE)

# GENOMIC OBJECT ----
Go <- loadranges("data/ChIPseq/ce11_tiled_25kb.bed", genome = "ce11")
mcols(Go)["chr_bin"] <- unlist(lapply(split(Go, seqnames(Go)), function(chr) seq_along(chr)))
Go <- bw_tile(bw = eigen_bw, bed = Go, genome = "ce11", name = eigen_name)
Go <- subset(Go, seqnames %in% chr_v)

# RUN ----
gg <- list()

for (cond in conditions) {
  h5      <- sprintf(h5_tpl, cond)
  cm_info <- rhdf5::h5ls(h5)
  cm_l    <- sapply(cm_info$name[cm_info$name %in% chr_v],
                    function(nm) HDF5Array::HDF5Array(h5, nm))

  sadmat <- getSaddle2D(bin_go = Go, bait = eigen_name, anchor = eigen_name,
                        diag_offset = 1, na_offset = 20,
                        bait_quant = N_BINS, anchor_quant = N_BINS,
                        cm_l = cm_l, f_score = "median",
                        f_sad = "mean", f_cm = "median", log2 = FALSE)

  ## source data: the genome-wide matrix, before the display rescaling
  m_raw <- as.data.table(sadmat[["sadmat_cm_mChr"]], keep.rownames = FALSE)
  fwrite(m_raw, file.path(out_d, sprintf("saddle_matrix.%s.csv", cond)))

  ## centre panel, rescaled to [-1, 1] for display
  m <- matrix2tibble(sadmat[["sadmat_cm_mChr"]], N_BINS, N_BINS)
  m$counts <- rangeMinMax(m$counts)
  gg_center <- ggplot(m, aes(x = Anchor, y = Bait, fill = counts)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.y = element_blank(), axis.title.x = element_blank()) +
    scale_x_discrete(position = "top") +
    scale_fill_gradient2(mid = "white", high = "firebrick", low = "dodgerblue")

  ## margins: median eigenvector value of each quantile group
  gg_bait <- ggplot(sadmat[["sum_stat_mChr_b"]], aes(x = factor(quant_bait), y = tot_bait)) +
    geom_bar(stat = "identity", aes(fill = tot_bait)) + theme_bw() +
    theme(legend.position = "none", axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90), axis.text.y = element_blank()) +
    scale_x_discrete(position = "bottom", limits = rev) +
    labs(y = "Tile median values", x = eigen_name) + coord_flip() +
    scico::scale_fill_scico(palette = "bamako", direction = 1)

  gg_anchor <- ggplot(sadmat[["sum_stat_mChr_a"]], aes(x = factor(quant_anchor), y = tot_anchor)) +
    geom_bar(stat = "identity", aes(fill = tot_anchor)) + theme_bw() +
    theme(legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
    labs(y = "Tile median values", x = eigen_name) + scale_x_discrete(position = "top") +
    scico::scale_fill_scico(palette = "bamako", direction = 1)

  tmp_p  <- cowplot::plot_grid(NULL, gg_anchor, gg_bait, gg_center + theme(legend.position = "none"),
                               rel_widths = c(1, 4), rel_heights = c(2, 4), ncol = 2, align = "hv")
  legend <- get_legend(gg_center + theme(legend.box.margin = margin(140, 0, 20)))
  gg[[sprintf("Saddle.%s.mergedChr", cond)]] <- cowplot::plot_grid(tmp_p, legend, rel_widths = c(5, 1))
}

# SAVE ----
for (nm in names(gg))
  ggsave(file.path(out_d, paste0(nm, ".png")), gg[[nm]], width = 8, height = 8, dpi = 150)
