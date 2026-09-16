# HEADER ====================================================== =
# GEN2I
# NAR revision 2026
# R version : R version 4.4.0 (2024-04-24), tested on 4.5.3
# System : x86_64, linux-gnu
# ============================================================= =
#
# Supplementary Figure S5C - quantification of the differential saddle plots of Figure 3B.
#
# For every genotype, bins are ranked by the wild-type PC2 and grouped into 50
# quantile groups; each saddle pixel is the mean of the mutant/wild-type ratios
# of the corresponding bin pairs, and the genome-wide matrix is the pixel-wise
# median of the six chromosome matrices. Three corners of 5 x 5 quantile groups
# are then summarized by the median of their 25 pixels:
#   B-B : top left      (rows 1-5,    columns 1-5)
#   A-A : bottom right  (rows 46-50,  columns 46-50)
#   B-A : top right     (rows 1-5,    columns 46-50)
# The figure shows two ratios of these medians, B-B/A-A and B-B/B-A, with a 95%
# bootstrap confidence interval (2000 resamplings of the corner pixels). A ratio
# of 1 means no change relative to wild type.
#
# Corner width: 5 quantile groups. The B-B corner has a steep gradient whereas
# the A-A corner is flat, so a wider corner dilutes the B-B signal and not the
# A-A one; with 10 groups the B-B/A-A comparison flips. The width profile from
# 3 to 15 groups is exported so that the choice can be checked.
#
# The pixels of one corner come from a single aggregated matrix per genotype, so
# they are not independent replicates: the interval is the precision of the
# estimate. The chromosome-level tables exported here treat the chromosome as
# the unit of observation (n = 6), which is the conservative reading.
#
# Figure 3C is the same analysis ranked by PC1: see figure_3C.R.
# ============================================================= =

# WORKING DIRECTORY ----
setwd(dir = "config/src/R/github/")

# SOURCE ----
source("manuscript_lib.R")

# PARAMETERS ----
out_d      <- "output/figure_S5C"
eigen_name <- "eigen_pca2_N2.old"                            # PC2 of wild type, initial batch
eigen_bw   <- "data/HiC/N2-old_merged.bwa_mem.25kb.pca2.bw"
conditions <- c("CEC4", "hpl2-lin61-old", "hpl2-old", "lin61-old", "met2-set25-set32-old")
h5_tpl     <- "data/HiC/%s_vs_N2-old.bwa_mem.25kb.norm.KR.g2i.h5"   # ratio mutant / wild type
chr_v      <- c("I", "II", "III", "IV", "V", "X")
N_BINS     <- 50L                                             # quantile groups per axis
WINDOW_W   <- 5L                                              # corner width, in quantile groups
NULL_VALUE <- 1                                               # ratio matrices: no change = 1
NBOOT      <- 2000
BOOT_SEED  <- 123
colPal     <- c("#666666", "#132CEC", "firebrick", "#D95F02", "#E7298A", "#66A61E", "#E6AB02")
dir.create(out_d, showWarnings = FALSE, recursive = TRUE)

# GENOMIC OBJECT ----
Go <- loadranges("data/ChIPseq/ce11_tiled_25kb.bed", genome = "ce11")
mcols(Go)["chr_bin"] <- unlist(lapply(split(Go, seqnames(Go)), function(chr) seq_along(chr)))
Go <- bw_tile(bw = eigen_bw, bed = Go, genome = "ce11", name = eigen_name)
Go <- subset(Go, seqnames %in% chr_v)

# SADDLE MATRICES ----
l_cm <- list()
for (cond in conditions) {
  h5 <- sprintf(h5_tpl, cond)
  cm_info <- rhdf5::h5ls(h5)
  l_cm[[cond]] <- sapply(cm_info$name[cm_info$name %in% chr_v],
                         function(nm) HDF5Array::HDF5Array(h5, nm))
}

sadmat_l <- list()
for (nm in names(l_cm))
  sadmat_l[[nm]] <- getSaddle2D(bin_go = Go, bait = eigen_name, anchor = eigen_name,
                                diag_offset = 1, na_offset = 20,
                                bait_quant = N_BINS, anchor_quant = N_BINS,
                                cm_l = l_cm[[nm]], f_score = "median",
                                f_sad = "mean", f_cm = "median", log2 = FALSE)

# CORNERS ----
zone_defs <- function(w) {
  b <- 1:w                        # head of the ranking: most negative PC = B compartment
  a <- (N_BINS - w + 1):N_BINS    # tail: most positive PC = A compartment
  list("B-B" = list(rows = b, cols = b),
       "A-A" = list(rows = a, cols = a),
       "B-A" = list(rows = b, cols = a))
}
ZONES      <- zone_defs(WINDOW_W)
ZONE_NAMES <- names(ZONES)
extract_zone <- function(dt, zone) dt[Row %in% zone$rows & Col %in% zone$cols][!is.na(Pixel_Value)]

## saddle matrix -> one row per pixel. Factor levels are converted through
## character: as.numeric() on a factor would return the level index, and the
## corners would be taken at the wrong positions.
saddle_to_dt <- function(sadmat) {
  m <- matrix2tibble(sadmat, N_BINS, N_BINS)
  data.table(Row = as.numeric(as.character(m$Bait)),
             Col = as.numeric(as.character(m$Anchor)),
             Pixel_Value = m$counts)
}

full_l <- list(); chrmed_l <- list()
for (nm in names(sadmat_l)) {
  full_l[[nm]] <- data.table(Condition = nm, saddle_to_dt(sadmat_l[[nm]][["sadmat_cm_mChr"]]))

  for (chr in sadmat_l[[nm]][["chr_v"]]) {                    # one median per chromosome and corner
    dchr <- saddle_to_dt(sadmat_l[[nm]][["sadmat_cm_aChr"]][[chr]])
    for (zn in ZONE_NAMES) {
      z <- extract_zone(dchr, ZONES[[zn]])
      if (nrow(z) > 0) chrmed_l[[length(chrmed_l) + 1]] <-
        data.table(Condition = nm, Chromosome = paste0("chr", chr), Interaction_Type = zn,
                   Median = median(z$Pixel_Value, na.rm = TRUE))
    }
  }
}
pix_matrix  <- rbindlist(full_l)
chr_medians <- rbindlist(chrmed_l)

pix_zones <- rbindlist(lapply(ZONE_NAMES, function(zn) {
  z <- extract_zone(pix_matrix, ZONES[[zn]])
  z[, .(Condition, Interaction_Type = zn, Row, Col, Pixel_Value)]
}))

## sensitivity of the two ratios to the corner width
window_profile <- rbindlist(lapply(3:15, function(w) {
  zs  <- zone_defs(w)
  med <- function(zn) extract_zone(pix_matrix, zs[[zn]])[, .(v = median(Pixel_Value, na.rm = TRUE)), by = Condition]
  bb <- med("B-B"); aa <- med("A-A"); ba <- med("B-A")
  setnames(bb, "v", "BB"); setnames(aa, "v", "AA"); setnames(ba, "v", "BA")
  out <- merge(merge(bb, aa, by = "Condition"), ba, by = "Condition")
  out[, `:=`(Window_W = w, N_pixels = w * w, BB_over_AA = BB / AA, BB_over_BA = BB / BA)][]
}))
setcolorder(window_profile, c("Condition", "Window_W", "N_pixels", "BB", "AA", "BA",
                              "BB_over_AA", "BB_over_BA"))

# RATIOS AND BOOTSTRAP ----
## the estimate is the ratio of the observed medians; the bootstrap only gives
## the interval
boot_ratio_draws <- function(num, den) {
  if (length(num) < 2 || length(den) < 2) return(rep(NA_real_, NBOOT))
  set.seed(BOOT_SEED)
  replicate(NBOOT, median(sample(num, replace = TRUE)) / median(sample(den, replace = TRUE)))
}

stats_ratio <- rbindlist(lapply(unique(pix_zones$Condition), function(cond) {
  bb <- pix_zones[Condition == cond & Interaction_Type == "B-B", Pixel_Value]
  aa <- pix_zones[Condition == cond & Interaction_Type == "A-A", Pixel_Value]
  ba <- pix_zones[Condition == cond & Interaction_Type == "B-A", Pixel_Value]
  rbindlist(lapply(list(list("B-B / A-A", bb, aa), list("B-B / B-A", bb, ba)), function(x) {
    dr <- boot_ratio_draws(x[[2]], x[[3]])
    data.table(Condition = cond, Ratio_Type = x[[1]],
               Ratio   = median(x[[2]], na.rm = TRUE) / median(x[[3]], na.rm = TRUE),
               CI_low  = quantile(dr, 0.025, na.rm = TRUE, names = FALSE),
               CI_high = quantile(dr, 0.975, na.rm = TRUE, names = FALSE))
  }))
}))
stats_ratio[, Excludes_1 := !is.na(CI_low) & (CI_low > NULL_VALUE | CI_high < NULL_VALUE)]

## corner medians, and the same ratios with the chromosome as the unit
stats_vs_wt <- pix_zones[, .(N_Pixels = .N,
                             Median_Ratio = median(Pixel_Value, na.rm = TRUE),
                             Percent_Change = 100 * (median(Pixel_Value, na.rm = TRUE) - NULL_VALUE)),
                         by = .(Condition, Interaction_Type)]

chr_ratios <- dcast(chr_medians, Condition + Chromosome ~ Interaction_Type, value.var = "Median")
chr_ratios[, `:=`(`B-B / A-A` = `B-B` / `A-A`, `B-B / B-A` = `B-B` / `B-A`)]
chr_ratios <- melt(chr_ratios[, .(Condition, Chromosome, `B-B / A-A`, `B-B / B-A`)],
                   id.vars = c("Condition", "Chromosome"),
                   variable.name = "Ratio_Type", value.name = "Ratio")
chr_ratio_stats <- chr_ratios[, {
  v <- Ratio[!is.na(Ratio)]
  set.seed(BOOT_SEED)
  dr <- if (length(v) >= 2) replicate(NBOOT, median(sample(v, replace = TRUE))) else rep(NA_real_, NBOOT)
  .(N_Chromosomes = .N, Median_Ratio = median(v),
    p_val = if (length(v) >= 2) suppressWarnings(wilcox.test(v, mu = NULL_VALUE)$p.value) else NA_real_,
    Min_Attainable_p = if (length(v) >= 2) 2 / (2^length(v)) else NA_real_,
    CI_low = quantile(dr, 0.025, na.rm = TRUE, names = FALSE),
    CI_high = quantile(dr, 0.975, na.rm = TRUE, names = FALSE))
}, by = .(Condition, Ratio_Type)]

# OUTPUT TABLES ----
fwrite(pix_zones,       file.path(out_d, "corner_pixels.csv"))
fwrite(stats_vs_wt,     file.path(out_d, "corner_medians.csv"))
fwrite(stats_ratio,     file.path(out_d, "corner_ratios.csv"))
fwrite(window_profile,  file.path(out_d, "corner_width_profile.csv"))
fwrite(chr_medians,     file.path(out_d, "chromosome_medians.csv"))
fwrite(chr_ratio_stats, file.path(out_d, "chromosome_ratios.csv"))

# PLOT ----
## point and interval, not a boxplot: a ratio of two medians is a single value,
## its spread comes from the resampling and not from a population of observations
gg <- list()
colPal_named <- setNames(colPal[seq_along(conditions)], conditions)
st <- copy(stats_ratio)[, lab := paste0("[", round(CI_low, 2), "-", round(CI_high, 2), "]")]
rng <- diff(range(c(st$CI_low, st$CI_high), na.rm = TRUE)); if (!is.finite(rng) || rng == 0) rng <- 1
st[, y_pos := max(c(CI_high), na.rm = TRUE) + 0.12 * rng]

gg[["SaddlePixel.Ratios"]] <-
  ggplot(st) +
  geom_hline(yintercept = NULL_VALUE, linetype = "dashed", colour = "purple", linewidth = 0.7) +
  geom_errorbar(aes(x = Condition, ymin = CI_low, ymax = CI_high), width = 0.16,
                linewidth = 0.8, colour = "#222222") +
  geom_point(aes(x = Condition, y = Ratio, fill = Condition), shape = 21, size = 3.4,
             stroke = 0.7, colour = "#222222") +
  geom_text(aes(x = Condition, y = y_pos, label = lab), size = 2.9, vjust = 0) +
  facet_wrap(~ Ratio_Type, ncol = 2) +
  scale_fill_manual(values = colPal_named) + guides(fill = "none") +
  labs(title = "Corner ratios",
       subtitle = paste0("all chromosomes | unit = corner pixel (", WINDOW_W^2, " per corner)",
                         " | point = estimate, bar = bootstrap 95% CI (", NBOOT, " draws)"),
       x = "Genotype / Condition", y = "Ratio (mutant / N2)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

gg[["SaddlePixel.WindowProfile"]] <-
  ggplot(window_profile, aes(x = Window_W, y = BB_over_AA, colour = Condition)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = WINDOW_W, linetype = "dotted", colour = "#132CEC") +
  geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
  scale_colour_manual(values = colPal_named) +
  scale_x_continuous(breaks = seq(3, 15, 2)) +
  labs(title = "Sensitivity of the B-B / A-A ratio to corner width",
       subtitle = paste0("dotted line = width used (", WINDOW_W, "); below 1: B-B more depleted than A-A"),
       x = "Corner width (bins)", y = "B-B / A-A", colour = "Genotype") +
  theme_minimal()

# SAVE ----
for (nm in names(gg))
  ggsave(file.path(out_d, paste0(nm, ".png")), gg[[nm]], width = 10, height = 6, dpi = 300)
