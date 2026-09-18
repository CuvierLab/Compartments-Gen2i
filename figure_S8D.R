# HEADER ====================================================== =
# GEN2I
# NAR revision 2026
# R version : R version 4.4.0 (2024-04-24), tested on 4.5.3
# System : x86_64, linux-gnu
# ============================================================= =
#
# Supplementary Figure S8D - high resolution heatmaps of the net change of
# ChIP-seq signal around H3K9me2/H3K9me3 peaks, mutant against wild type.
#
# Same analysis as Figure 4C (see figure_4C.R for the full description), for the
# two single mutants hpl-2 and lin-61.
#
# Regions : union of the H3K9me2 and H3K9me3 narrow peaks called in N2
#           (6182 peaks), restricted to compartment B, i.e. to the peaks whose
#           nearest 25 kb bin has a negative wild-type PC1 (eigen_pca1_N2.old,
#           initial Hi-C batch) -> 4737 peaks.
# Signal  : deeptools scale-regions profiles of the z-score bigwigs
#           zscore(mutant, N2), 2 kb upstream + 200 bp scaled body + 2 kb
#           downstream in 10 bp bins (420 columns). Positive values mean more
#           of the mark of interest in the mutant.
# Ranking : peaks are sorted by increasing H3K9me2 z-score signal of the SAME
#           mutant over the peak body, so the IPs of a genotype share one row
#           order and can be read side by side.
# Display : the 4737 x 420 matrix is block-averaged to 800 rows x 300 columns
#           (redim_matrix), values are clipped to [-2, 2] and drawn on a fixed
#           viridis scale, identical for every IP and every genotype.
#
# HPL-2 was not profiled in the hpl-2 background, so that panel does not exist.
# The LIN-61 panel of the lin-61 mutant is also printed in Figure 4C, under
# lin-61; hpl-2 and marked with an asterisk, because LIN-61 was not profiled in
# the double mutant: the two panels are the same image.
# ============================================================= =

# WORKING DIRECTORY ----
# Run the scripts from the root of this repository (where manuscript_lib.R sits):
# setwd("/path/to/Compartments-Gen2i")

# SOURCE ----
source("manuscript_lib.R")

# PARAMETERS ----
out_d      <- "output/figure_S8D"
peaks_bed  <- "data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed"  # union H3K9me2 + H3K9me3, N2
tiles_bed  <- "data/ChIPseq/ce11_tiled_25kb.bed"
eigen_name <- "eigen_pca1_N2.old"                                          # PC1 of wild type, initial batch
eigen_bw   <- "data/HiC/N2-old_merged.bwa_mem.25kb.pca1.bw"
prof_tpl   <- "data/ChIPseq/%s_%s_newzs_prof.txt"                          # ip, condition (dots)
zbw_tpl    <- "data/ChIPseq/%s_%s_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw"# ip, condition (dashes)
ranker_ip  <- "H3K9me2"                                                    # IP used to sort the rows
HEIGHT     <- 800L                                                         # rows after block-averaging
WIDTH      <- 300L                                                         # columns after block-averaging
CLIP       <- c(-2, 2)                                                     # value clipping and colour range

# Panels of Supplementary Figure S8D, in the order they appear in the figure.
# prof_cond is kept for symmetry with figure_4C.R: here it is always the
# genotype of the panel, no profile is borrowed from another genotype.
panels <- data.table(
  cond      = c("hpl2", "hpl2", "hpl2",
                "lin61", "lin61", "lin61", "lin61"),
  ip        = c("H3K9me2", "lin61", "LEM2",
                "H3K9me2", "hpl2", "lin61", "LEM2"),
  prof_cond = c("hpl2", "hpl2", "hpl2",
                "lin61", "lin61", "lin61", "lin61")
)
dir.create(out_d, showWarnings = FALSE, recursive = TRUE)

# GENOMIC OBJECT ----
Go <- loadranges(peaks_bed, genome = config.genome)

## ChIP-seq z-score profile matrices -> mcols (one DataFrame of 420 columns each)
for (i in seq_len(nrow(panels))) {
  aname <- paste0(panels$ip[i], "_", panels$prof_cond[i], "_newzs_prof")
  if (aname %in% colnames(mcols(Go))) next()
  prof_dt <- data.table::fread(sprintf(prof_tpl, panels$ip[i], panels$prof_cond[i]),
                               check.names = FALSE)[, name := as.character(name)]
  data.table::setkey(prof_dt, name)                       # retrieve ordering by name
  mcols(Go)[[aname]] <- prof_dt[as.character(Go$name), !"name"]   # match the peak order
}

## H3K9me2 z-score signal over the peak body -> mcols (the ranking key)
for (cond in unique(panels$prof_cond)) {
  sname <- paste0(ranker_ip, "_", cond, "_newzs_signal")
  if (sname %in% colnames(mcols(Go))) next()
  Go <- bw_to_mcols(bw = sprintf(zbw_tpl, ranker_ip, gsub("\\.", "-", cond)),
                    bed = Go, genome = config.genome, anchor = "body",
                    name = sname, upstream = 0, downstream = 0)
}

## Wild-type PC1 of the nearest 25 kb bin -> mcols (the compartment filter)
Go2 <- loadranges(tiles_bed, genome = config.genome)
Go2 <- bw_tile(bw = eigen_bw, bed = Go2, genome = config.genome, name = eigen_name)
idx <- GenomicRanges::nearest(Go, Go2)
mcols(Go)[eigen_name] <- NA
na_idx <- which(is.na(idx))
if (length(na_idx)) {
  mcols(Go[-na_idx])[eigen_name] <- unlist(mcols(Go2[idx[!is.na(idx)]])[eigen_name])
} else {
  mcols(Go)[eigen_name] <- unlist(mcols(Go2[idx])[eigen_name])
}

## Compartment B
Go <- subset(Go, eigen_pca1_N2.old <= 0)

# RUN ----
gg <- list()

for (i in seq_len(nrow(panels))) {
  aname <- paste0(panels$ip[i], "_", panels$prof_cond[i], "_newzs_prof")
  order_n <- paste0(ranker_ip, "_", panels$prof_cond[i], "_newzs_signal")
  nm <- sprintf("Heatmap.Condition_%s.IP_%s", panels$cond[i], panels$ip[i])
  if (panels$prof_cond[i] != panels$cond[i])                    # the borrowed panel
    nm <- sprintf("Heatmap.Condition_%s.IP_%s_from_%s", panels$cond[i], panels$ip[i], panels$prof_cond[i])

  m <- mcols(Go)[[aname]]
  height <- min(HEIGHT, nrow(m))
  width  <- min(WIDTH,  ncol(m))

  ## rows sorted by increasing H3K9me2 z-score, then block-averaged
  f  <- order(mcols(Go)[[order_n]])
  rd <- redim_matrix(as.matrix(m[f, ]), target_height = height, target_width = width, n_core = 1)

  ## source data: the matrix that is drawn, before clipping
  fwrite(as.data.table(rd), file.path(out_d, paste0(nm, ".csv")))

  g <- reshape2::melt(rd)
  min.val <- min(g$value, na.rm = TRUE)
  setDT(g)
  g[is.na(value), value := min.val]
  g[value < min(CLIP), value := min(CLIP)]
  g[value > max(CLIP), value := max(CLIP)]

  gg[[nm]] <- ggplot(g, aes(x = Var2, y = Var1, fill = value)) +
    geom_raster() + theme_minimal() +
    theme(plot.margin = margin(grid::unit(0, "cm")),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = element_blank(),
          plot.caption = element_text(hjust = 0, size = 8, face = "italic"),
          plot.subtitle = element_text(hjust = 0, size = 8),
          plot.title = element_text(hjust = 0, size = 12, face = "bold")) +
    scale_fill_viridis_c(guide = "colorbar", limits = CLIP) +
    labs(x = "Position", y = "Sites", fill = "Diff. score")
}

# SAVE ----
# width/80 x height/80 inches at 300 dpi, the size the published panels were drawn at.
for (nm in names(gg))
  ggsave(file.path(out_d, paste0(nm, ".png")), gg[[nm]],
         width = WIDTH / 80, height = HEIGHT / 80, dpi = 300, limitsize = FALSE)
