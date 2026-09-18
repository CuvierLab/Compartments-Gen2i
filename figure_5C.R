# HEADER ====================================================== =
# GEN2I
# NAR revision 2026
# R version : R version 4.4.0 (2024-04-24), tested on 4.5.3
# System : x86_64, linux-gnu
# ============================================================= =
#
# Figure 5C - average ChIP-seq peak density around the sites bound by both HPL-2
# and LIN-61 in wild type.
#
# The reference set is the intersection of the HPL-2 and LIN-61 narrowPeak sets
# in N2. Each reference site is reduced to one bin at its centre, extended by
# `win_lim` on each side and tiled at `binsize`; sites whose window runs off the
# end of a chromosome are dropped, so that every position of the profile is
# supported by the same number of sites. For each profiled peak set (HPL-2,
# LIN-61, H3K9me2, H3K9me3) the number of overlapping peaks is counted in every
# tile and summed over all reference sites, then divided by the total genomic
# span of that peak set in Mb.
#
# That normalization is what makes the four antibodies comparable on one axis:
# without it a factor with more, or broader, peaks would sit higher everywhere
# for a reason that has nothing to do with the reference sites.
#
# The central bin is a self-overlap: the profiled HPL-2 and LIN-61 peaks contain,
# by construction, the reference set that is their own intersection. It therefore
# rises far above the flanks, and on a linear scale it flattens the flanking
# structure the panel is about. The published panel is the log10 view, which
# shows the spike and the flanks together without cropping anything; the linear
# views are written as well.
# ============================================================= =

# WORKING DIRECTORY ----
# Run the scripts from the root of this repository (where manuscript_lib.R sits):
# setwd("/path/to/Compartments-Gen2i")

# SOURCE ----
source("manuscript_lib.R")

# PARAMETERS ----
out_d    <- "output/figure_5C"
ref_bed  <- "data/ChIPseq/hpl2_lin61_N2_sharp_intersect_peakset.bed"   # HPL-2 and LIN-61 double sites
prof_bed <- c(hpl2    = "data/ChIPseq/hpl2_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
              lin61   = "data/ChIPseq/lin61_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
              H3K9me2 = "data/ChIPseq/H3K9me2_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
              H3K9me3 = "data/ChIPseq/H3K9me3_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed")
win_lim  <- 150e3                                            # half-window, bp
binsize  <- 1e3
fix      <- "center"
zooms_kb <- c(50, 25, 10)                                    # views of the same profile
dir.create(out_d, showWarnings = FALSE, recursive = TRUE)

binnumber <- 2 * win_lim / binsize + 1
x_seq     <- -((binnumber - 1) / 2):((binnumber - 1) / 2)

# GENOMIC OBJECT ----
Go_ref_raw <- loadranges(ref_bed, genome = "ce11")
Go_prof_l  <- lapply(prof_bed, function(p) loadranges(p, genome = "ce11"))

# RUN ----
## one bin at the centre of each reference site, extended and tiled
Go_ref <- GenomicRanges::trim(GenomicRanges::resize(Go_ref_raw, width = binsize, fix = fix))
Go_ref <- Go_ref[width(Go_ref) == binsize]

Go_ref_xt       <- trim(Go_ref + win_lim)
Go_ref_xt       <- Go_ref_xt[width(Go_ref_xt) == binsize + 2 * win_lim]  # complete windows only
Go_ref_xt_tiled <- unlist(GenomicRanges::tile(Go_ref_xt, width = binsize))

all_dt  <- NULL
norm_dt <- NULL                                              # normalize by the total peak span
for (ip in names(Go_prof_l)) {
  Go_prof_s <- Go_prof_l[[ip]]
  norm_dt   <- rbind(norm_dt, data.table(name = ip, norm = sum(width(Go_prof_s)) / 1e6))
  col       <- countOverlaps(Go_ref_xt_tiled, Go_prof_s)
  all_dt    <- rbind(all_dt, data.table(bin = rep(x_seq, length(Go_ref_xt)), value = col, name = ip))
}

val_dt <- all_dt[, .(svalue = sum(value)), by = .(bin, name)]
data.table::setkey(val_dt, name)
data.table::setkey(norm_dt, name)
plot_dt <- norm_dt[val_dt]
plot_dt[, normalized_density := svalue / norm]
plot_dt[, ip := name]
plot_dt[, kb := bin * binsize / 1e3]
setorder(plot_dt, ip, bin)

## source data: the profile itself, and the counts it is built from
fwrite(plot_dt[, .(ip, bin, kb, raw_count = svalue, peak_span_Mb = norm, normalized_density)],
       file.path(out_d, "profile.csv"))
fwrite(data.table(n_reference_sites = length(Go_ref_raw),
                  n_sites_used      = length(Go_ref_xt),
                  win_lim_bp        = win_lim,
                  binsize_bp        = binsize),
       file.path(out_d, "profile_parameters.csv"))

## plots: one view of the profile per half-window. The framing uses
## coord_cartesian and never xlim/ylim, so no point is dropped and the y scale
## is not recomputed: the central spike is never cropped.
profile_plot <- function(dt, half_win_kb, log_y = FALSE) {
  p <- ggplot(dt, aes(x = kb, y = normalized_density)) +
    geom_vline(xintercept = 0, colour = "#BBBBBB", linewidth = 0.4) +
    geom_step(aes(col = ip), linewidth = 0.6) +
    coord_cartesian(xlim = c(-half_win_kb, half_win_kb)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    labs(x = "Distance to peak centre (kb)",
         y = paste0("Normalized ChIP-seq peak density\n(", binsize / 1e3, " kb bins)"),
         colour = NULL) +
    theme_classic() +
    theme(legend.position = "top",
          axis.text  = element_text(size = 9),
          axis.title = element_text(size = 10))
  if (log_y)
    p <- p + scale_y_log10() +
      labs(y = paste0("Normalized ChIP-seq peak density, log10\n(", binsize / 1e3, " kb bins)"))
  p
}

gg     <- list()
win_kb <- win_lim / 1e3
zooms  <- Filter(function(z) z < win_kb, zooms_kb)

gg[["Density_profile.All_features"]] <- profile_plot(plot_dt, win_kb)
for (z in zooms)
  gg[[paste0("Density_profile.Zoom_", z, "kb.All_features")]] <- profile_plot(plot_dt, z)

## log10 views; "Log.Zoom_25kb" is the panel of the figure
gg[["Density_profile.Log.All_features"]] <- profile_plot(plot_dt, win_kb, log_y = TRUE)
for (z in zooms)
  gg[[paste0("Density_profile.Log.Zoom_", z, "kb.All_features")]] <-
    profile_plot(plot_dt, z, log_y = TRUE)

# SAVE ----
for (nm in names(gg))
  ggsave(file.path(out_d, paste0(nm, ".png")), gg[[nm]], width = 6, height = 4.5, dpi = 300)
