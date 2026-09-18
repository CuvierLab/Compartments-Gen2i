# HEADER ====================================================== =
# GEN2I
# NAR revision 2026
# R version : R version 4.4.0 (2024-04-24), tested on 4.5.3
# System : x86_64, linux-gnu
# ============================================================= =
#
# Figure 2D - chromatin marks in compartment B against compartment A.
#
# The 25 kb tiles are ranked by the wild-type PC1 of the initial batch and cut
# into 50 quantile groups; groups 1 to 12 are compartment B, groups 16 to 50 are
# compartment A, and the four groups in between are dropped (this is the
# "onlyBandA" definition of the pipeline, the one used by this panel; the
# stricter "fromSaddle" definition, B = groups 1 to 3, is the one used by
# Supplementary Fig. S8E).
#
# For each mark, the value carried by a tile is the maximum of the wild-type
# ChIP-seq signal over the peaks of that mark that the tile overlaps, i.e. the
# annotation the pipeline calls <mark>_N2_max_signal. The published y axis reads
# "relative peak density"; the quantity actually plotted is that maximum ChIP-seq
# signal per bin. We keep the computation as published and only flag the label.
#
# As in the pipeline, values beyond 1.5 IQR of their own compartment are set to
# NA before the boxplot is drawn (remove_outliers, applied by compartment). The
# boxplot then recomputes its own whiskers, which is why some points are still
# drawn as outliers. The published panel is that version.
#
# The p values printed on the published panel are two-sided Wilcoxon rank sum
# tests between the two compartments; they are written to
# compartment_summary.csv.
# ============================================================= =

# WORKING DIRECTORY ----
# Run the scripts from the root of this repository (where manuscript_lib.R sits):
# setwd("/path/to/Compartments-Gen2i")

# SOURCE ----
source("manuscript_lib.R")

# PARAMETERS ----
out_d      <- "output/figure_2D"
tiles_p    <- "data/ChIPseq/ce11_tiled_25kb.bed"
eigen_name <- "eigen_pca1_N2.old"                            # PC1 of wild type, initial batch
eigen_bw   <- "data/HiC/N2-old_merged.bwa_mem.25kb.pca1.bw"
grp_name   <- "eigen_pca1_N2.old_grps_shifted_A_B"
N_BINS     <- 50L                                            # quantile groups
CUT_AT     <- c(12L, 16L)                                    # B <= 12 < na <= 16 < A
CUT_NAMES  <- c("B", "na", "A")
# mark -> peak set (our own peaks) and wild-type ChIP-seq signal
marks <- c("H3K4me3", "H3K9me2", "LEM2", "hpl2", "lin61")
peak_tpl   <- "data/ChIPseq/%s_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed"
signal_tpl <- "data/ChIPseq/%s_N2_merged.bwa_aln.rmdup.bamCompare.bw"
dir.create(out_d, showWarnings = FALSE, recursive = TRUE)

# RUN ----
## MCOLS Eigen values and compartments ----
Go <- loadranges(tiles_p, genome = "ce11")
Go <- bw_tile(bw = eigen_bw, bed = Go, genome = config.genome, name = eigen_name)

# 50 quantile groups over the whole genome, then the A / na / B cut. NA
# eigenvector values are treated as 0 by the pipeline, hence the same here.
q <- dplyr::ntile(mcols(Go)[[eigen_name]], N_BINS)
q[is.na(q)] <- 0L
mcols(Go)[[grp_name]] <- ifelse(q <= CUT_AT[1], CUT_NAMES[1],
                         ifelse(q <= CUT_AT[2], CUT_NAMES[2], CUT_NAMES[3]))

## MCOLS ChIP-seq signal, maximum over the peaks of each tile ----
for (mark in marks) {
  Go2 <- loadranges(sprintf(peak_tpl, mark), genome = "ce11")
  Go2 <- bw_to_mcols(bw = sprintf(signal_tpl, mark), bed = Go2,
                     genome = config.genome, anchor = "body",
                     name = paste0(mark, "_N2_signal"),
                     upstream = 0, downstream = 0)
  Go <- aggregate_ranges(Go, Go2,
                         out_mc = paste0(mark, "_N2_max_signal"),
                         subject_mc = paste0(mark, "_N2_signal"),
                         fun_agr = "max", fun_opts = "na.rm=T")
}

## PLOTS ----
Go_dt <- data.table(data.frame(mcols(Go)))
Go_dt <- Go_dt[!get(grp_name) %in% c("na", "Bp"), ]   # "onlyBandA", as in the pipeline

gg <- list()
summ <- NULL
for (mark in marks) {
  y <- paste0(mark, "_N2_max_signal")

  # The pipeline draws the boxplot on the values left after remove_outliers,
  # applied compartment by compartment. A tile that overlaps no peak of the mark
  # gets max(numeric(0)) = -Inf, which is what the pipeline stores as well; those
  # tiles carry no information and ggplot drops them when it draws the boxplot.
  # remove_outliers is run before that, on the -Inf included, exactly as in the
  # pipeline: when more than a quarter of a compartment has no peak the first
  # quartile is itself -Inf and nothing is trimmed.
  d  <- Go_dt[, .SD, .SDcols = c(grp_name, y)]
  dn <- data.table::copy(d)
  dn[, (y) := remove_outliers(get(y)), by = get(grp_name)]
  drawn <- is.finite(dn[[y]])

  # SOURCE TABLE ----
  # raw value per tile, and whether it survives to the boxplot
  fwrite(cbind(d, data.table(drawn = drawn)),
         file.path(out_d, sprintf("tile_values.%s.csv", mark)))
  dn <- dn[drawn, ]

  w <- wilcox.test(dn[[y]] ~ dn[[grp_name]])
  summ <- rbind(summ, data.table(
    mark = mark,
    n_A = sum(dn[[grp_name]] == "A"),
    n_B = sum(dn[[grp_name]] == "B"),
    n_no_peak = sum(!is.finite(d[[y]])),
    median_A = median(dn[get(grp_name) == "A"][[y]]),
    median_B = median(dn[get(grp_name) == "B"][[y]]),
    wilcox_p = w$p.value))

  # The boxplot block of the pipeline builds the fill scale from the grouping
  # column alone, so the two compartments get the default ggplot colours; the
  # theme is the one the pipeline sets globally (theme_tufte).
  gg[[sprintf("Boxplot.%s.rmOutliers", mark)]] <-
    ggplot(dn, aes(x = .data[[grp_name]], y = .data[[y]],
                   fill = .data[[grp_name]])) +
    geom_boxplot() + ggthemes::theme_tufte()
}

fwrite(summ, file.path(out_d, "compartment_summary.csv"))

# SAVE ----
for (nm in names(gg))
  ggsave(file.path(out_d, paste0(nm, ".png")), gg[[nm]],
         width = 8, height = 8, dpi = 100)
