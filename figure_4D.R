# HEADER ====================================================== =
# GEN2I
# NAR revision 2026
# R version : R version 4.4.0 (2024-04-24), tested on 4.5.3
# System : x86_64, linux-gnu
# ============================================================= =
#
# Figure 4D - correlation network of the differential ChIP-seq signals over the
# H3K9me2/H3K9me3/LIN-61/HPL-2 union peak set, restricted to the B compartment.
#
# One node per (antibody, mutant) pair: the z-scored mutant/wild-type ChIP-seq
# ratio of that antibody in that mutant, measured on each peak of the union set.
# Two nodes are linked when the Pearson correlation of their per-peak values is
# at least 0.3 in absolute value; the width and the colour of an edge give that
# correlation, blue for positive, red for negative. Node positions come from a
# classical multidimensional scaling of 1 - |r|, so nodes that respond alike sit
# close together: it is a map of co-variation, not a physical layout.
#
# Peaks are assigned to a compartment through the 25 kb bin nearest to them:
# bins are ranked by the wild-type PC1 (merged batch) into 50 quantile groups,
# groups 1-12 are B, groups 17-50 are A, and the groups in between are left
# unassigned. Only B is shown in the figure; the script also writes the A and
# the unfiltered networks.
# ============================================================= =

# WORKING DIRECTORY ----
# Run the scripts from the root of this repository (where manuscript_lib.R sits):
# setwd("/path/to/Compartments-Gen2i")

# SOURCE ----
source("manuscript_lib.R")
library(ggplotify)

# PARAMETERS ----
out_d      <- "output/figure_4D"
peak_bed   <- "data/ChIPseq/H3K9me23_lin61_hpl2_optimal_all_cond_peakset.bed"
zs_tpl     <- "data/ChIPseq/%s_%s_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw"
eigen_name <- "eigen_pca1_N2"                                # PC1 of wild type, merged batch
eigen_bw   <- "data/HiC/N2_merged.bwa_mem.25kb.pca1.bw"
comp_name  <- "eigen_pca1_N2_grps_shifted_A_B"
n_groups   <- 50L                                            # quantile groups of the eigenvector
cut_values <- c(12.0, 16.0)                                  # B = groups 1-12, A = groups 17-50
cut_names  <- c("B", "na", "A")
min_cor    <- 0.3                                            # an edge is drawn above this |r|
repel_seed <- 123                                            # ggrepel only moves the labels
## antibody -> mutants, in the order of the pipeline. I158A is excluded from the
## figures, and CEC4 has no differential ChIP-seq track.
ip_cond_l  <- list(
  hpl2     = c("hpl2.lin61", "lin61", "met2.set25.set32"),   # no hpl2 track in the hpl-2 mutant
  lin61    = c("hpl2", "lin61", "met2.set25.set32"),
  H3K9me2  = c("hpl2", "hpl2.lin61", "lin61", "met2.set25.set32"),
  H3K9me3  = c("hpl2", "hpl2.lin61", "lin61", "met2.set25.set32"))
dir.create(out_d, showWarnings = FALSE, recursive = TRUE)

# GENOMIC OBJECT ----
## union peak set, one mcol per (antibody, mutant) z-score track
Go <- loadranges(peak_bed, genome = "ce11")
for (ip in names(ip_cond_l)) {
  for (cond in ip_cond_l[[ip]]) {
    aname   <- paste0(ip, "_", cond, "_newzs_signal")
    c_nodot <- gsub("\\.", "-", cond)
    Go <- bw_to_mcols(bw = sprintf(zs_tpl, ip, c_nodot), bed = Go, genome = "ce11",
                      anchor = "body", name = aname, upstream = 0, downstream = 0)
  }
}

## tiled genome carrying the compartment call
Go2 <- loadranges("data/ChIPseq/ce11_tiled_25kb.bed", genome = "ce11")
Go2 <- bw_tile(bw = eigen_bw, bed = Go2, genome = "ce11", name = eigen_name)

### 50 quantile groups of the eigenvector, then A / na / B
mcols(Go2)["eigen_50tile"] <- dplyr::ntile(unlist(mcols(Go2)[[eigen_name]]), n_groups)
f <- is.na(mcols(Go2)[["eigen_50tile"]])
mcols(Go2)[["eigen_50tile"]][f] <- 0
mcols(Go2)[comp_name] <- NA
for (i in seq_along(cut_values)) {
  if (i == 1) f <- as.vector(mcols(Go2)[["eigen_50tile"]] <= cut_values[i])
  else f <- as.vector(mcols(Go2)[["eigen_50tile"]] > cut_values[i - 1] &
                      mcols(Go2)[["eigen_50tile"]] <= cut_values[i])
  if (any(f)) mcols(Go2[f])[comp_name] <- cut_names[i]
}
f <- as.vector(mcols(Go2)[["eigen_50tile"]] > cut_values[length(cut_values)])
if (any(f)) mcols(Go2[f])[comp_name] <- cut_names[length(cut_values) + 1]

### each peak inherits the compartment of the nearest 25 kb bin. Peaks with no
### nearest bin (the MtDNA peaks, absent from the tiled genome) keep NA and are
### left out of both compartments.
idx <- GenomicRanges::nearest(Go, Go2)
mcols(Go)[comp_name] <- NA
Go_na_idx <- which(is.na(idx))
if (length(Go_na_idx)) {
  idx <- idx[!is.na(idx)]
  mcols(Go[-Go_na_idx])[comp_name] <- unlist(mcols(Go2[idx])[comp_name])
} else {
  mcols(Go)[comp_name] <- unlist(mcols(Go2[idx])[comp_name])
}
stopifnot(all(c("A", "B") %in% mcols(Go)[[comp_name]]))

# RUN ----
gg    <- list()
Go_dt <- go2dt(Go)

## the seed only drives the label repulsion; the node coordinates come from the
## multidimensional scaling and are deterministic
old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL

for (comp in c("ALL", "A", "B")) {
  Go_s_dt <- if (comp != "ALL") subset(Go_dt, get(comp_name) == comp) else Go_dt
  Go_s_dt <- Go_s_dt[, .SD, .SDcols = colnames(Go_s_dt) %like% "newzs"]

  cor_df <- correlate(Go_s_dt)

  ## source data: the correlation matrix and the edges actually drawn
  fwrite(cor_df, file.path(out_d, sprintf("correlation_matrix.Compartment_%s.csv", comp)))
  m <- as.matrix(cor_df[, -1]); rownames(m) <- cor_df$term
  ed <- as.data.table(which(lower.tri(m) & abs(m) >= min_cor, arr.ind = TRUE))
  fwrite(data.table(from = rownames(m)[ed$row], to = colnames(m)[ed$col],
                    r = m[cbind(ed$row, ed$col)], n_peaks = nrow(Go_s_dt)),
         file.path(out_d, sprintf("network_edges.Compartment_%s.csv", comp)))

  set.seed(repel_seed)
  gg[[paste0("Network_analyzes.Label.Compartment_Filter_", comp)]] <-
    ggplotify::as.ggplot(network_plot(cor_df, min_cor = min_cor, curved = FALSE, repel = TRUE))

  ## unlabelled variant: names replaced by blanks of increasing length
  Go_b_dt <- copy(Go_s_dt)
  for (i in seq_along(names(Go_b_dt)))
    setnames(Go_b_dt, old = names(Go_b_dt)[i], new = paste(rep(" ", i), collapse = ""))
  set.seed(repel_seed)
  gg[[paste0("Network_analyzes.NoLabel.Compartment_Filter_", comp)]] <-
    ggplotify::as.ggplot(Go_b_dt %>% correlate() %>%
                           network_plot(min_cor = min_cor, curved = FALSE, repel = TRUE))
}

## leave the global random stream as it was found
if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
  rm(".Random.seed", envir = .GlobalEnv)

# SAVE ----
for (nm in names(gg))
  ggsave(file.path(out_d, paste0(nm, ".png")), gg[[nm]], width = 8, height = 8, dpi = 150)
