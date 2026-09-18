# HEADER ====================================================== =
# GEN2I
# NAR revision 2026
# R version : R version 4.4.0 (2024-04-24), tested on 4.5.3
# System : x86_64, linux-gnu
# ============================================================= =
#
# Supplementary Figure S5D - cumulative plots of intrachromosomal contacts by interaction
# class, and quantification of the contact preferences within each map.
#
# Left  : empirical cumulative distributions of z-scored log2(observed/expected)
#         contacts, one curve per genotype, for B-B, A-B and A-A contacts of the
#         autosomes (chromosomes I-III) and B-B contacts of the X chromosome.
# Right : changes relative to wild type of quantities measured WITHIN each map
#         on raw log2(observed/expected) values:
#           B-B preference  = median(B-B) - median(A-B)
#           A-A preference  = median(A-A) - median(A-B)
#           A-B segregation = median(A-B) - mean(median(A-A), median(B-B))
#         Differences are assessed by an m-out-of-n bootstrap (m = 1000 contacts
#         per class and genotype, B = 2000 iterations). The bracket compares the
#         mean of the two mutants with impaired compartmentalization with the
#         mean of wild type and the three single mutants, on the same draws.
#
# The z-score of the left panels is computed over all contacts of a map, with
# the three classes pooled: it says where a class sits relative to the others.
# A loss of B-B contacts therefore shifts A-A and A-B upwards together. The
# right panels compare classes within one map, where such a common shift
# cancels out; this is why A-A can move on the left and not on the right.
#
# Figure 3E is the same analysis ranked by PC1: see figure_3E.R.
# ============================================================= =

# WORKING DIRECTORY ----
# Run the scripts from the root of this repository (where manuscript_lib.R sits):
# setwd("/path/to/Compartments-Gen2i")

# SOURCE ----
source("manuscript_lib.R")
library(matrixStats)

# PARAMETERS ----
out_d      <- "output/figure_S5D"
res        <- "25kb"
norm       <- "obs_exp"
eigen_name <- "eigen_pca2_N2.old"                            # PC2 of wild type, initial batch
eigen_bw   <- "data/HiC/N2-old_merged.bwa_mem.25kb.pca2.bw"
gp_xp      <- c("N2-old", "CEC4", "hpl2-lin61-old", "hpl2-old", "lin61-old", "met2-set25-set32-old")
ref_xp     <- "N2-old"
gp_chr_l   <- list(chrom_I_II_III = c("I", "II", "III"),     # "autosomes" of the figure
                   chrom_X        = "X")
n_boot     <- 2000                                            # bootstrap iterations
n_sub      <- 1000                                            # contacts drawn per class and genotype
boot_seed  <- 20260911
dir.create(out_d, showWarnings = FALSE, recursive = TRUE)

# GENOMIC OBJECT ----
## bins of 25 kb, index within chromosome, wild-type eigenvector on each bin
Go <- loadranges("data/ChIPseq/ce11_tiled_25kb.bed", genome = "ce11")
mcols(Go)["chr_bin"] <- unlist(lapply(split(Go, seqnames(Go)), function(chr) seq_along(chr)))
Go <- bw_tile(bw = eigen_bw, bed = Go, genome = "ce11", name = eigen_name)

go_dt <- go2dt(Go)
go_dt[, class := ifelse(get(eigen_name) < 0, "B", "A")]       # B compartment = negative eigenvector

# RUN ----
gg <- list()

for (chrom_idx in names(gp_chr_l)) {

  ## ---- contacts of each class, genotype by genotype ------------------------
  dat <- NULL
  for (xp in gp_xp) {
    cm_l <- readRDS(sprintf("data/HiC/%s_%s.%s.cm.rds", xp, res, norm))

    for (chr in gp_chr_l[[chrom_idx]]) {
      cm <- as.matrix(cm_l[[chr]])

      ## bins whose coverage is an outlier are masked, rows and columns alike
      Q    <- quantile(rowSums(cm, na.rm = TRUE), c(0.05, 0.95), na.rm = TRUE)
      iqr  <- IQR(rowSums(cm, na.rm = TRUE), na.rm = TRUE)
      keep <- rowSums(cm) > (Q[1] - 1.5 * iqr) & rowSums(cm) < (Q[2] + 1.5 * iqr)
      cm[!keep, ] <- NA; cm[, !keep] <- NA

      ## the main diagonal and the first two off-diagonals carry proximity, not
      ## compartmentalization
      getKDiag(cm, -2:2) <- NA

      ## extreme contacts are masked as well, then values are log2-transformed
      Q   <- quantile(cm, c(0.05, 0.95), na.rm = TRUE)
      iqr <- IQR(cm, na.rm = TRUE)
      cm[cm > (Q[2] + 1.5 * iqr) | cm < (Q[1] - 1.5 * iqr)] <- NA
      cm[cm == 0] <- NA
      cm <- log2(cm)

      idx.A <- go_dt[seqnames == chr & class == "A", chr_bin]
      idx.B <- go_dt[seqnames == chr & class == "B", chr_bin]

      AA <- cm[idx.A, idx.A]; AA[lower.tri(AA)] <- NA        # upper triangle only
      BB <- cm[idx.B, idx.B]; BB[lower.tri(BB)] <- NA
      AB <- cm[idx.A, idx.B]

      dat <- rbind(dat,
                   data.table(value = as.vector(AA), class = "A to A", xp = xp, chr = chr),
                   data.table(value = as.vector(BB), class = "B to B", xp = xp, chr = chr),
                   data.table(value = as.vector(AB), class = "A to B", xp = xp, chr = chr))
    }
  }

  ## ---- left panels: cumulative distributions of the z-scores ---------------
  dat[, zscore := scale(value), by = xp]                     # one z-score per genotype, classes pooled
  dat <- dat[is.finite(zscore)]

  for (clss in unique(dat$class)) {
    gg[[sprintf("Lineplot.%s.%s", clss, chrom_idx)]] <-
      ggplot(dat[class == clss], aes(x = zscore, colour = xp, group = xp)) +
      stat_ecdf(geom = "line", n = 400) +
      coord_cartesian(xlim = c(-4, 4)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
      scale_colour_brewer(palette = "Set1") +
      labs(x = "Z-score", y = "ECDF", title = sprintf("%s (%s)", clss, chrom_idx)) +
      theme_classic2()
  }

  ## ---- right panels: preferences measured within each map ------------------
  dat_f <- dat[is.finite(value)]
  three <- c("A to A", "B to B", "A to B")
  if (!(all(three %in% dat_f$class) && ref_xp %in% dat_f$xp)) next()

  cls3    <- c("B to B", "A to B", "A to A")
  dcol_of <- c(`B to B` = "delta_BB", `A to A` = "delta_AA", `A to B` = "delta_AB")
  ylab_of <- c(`B to B` = "Δ preference B-B vs A-B (mutant − N2)",
               `A to B` = "Δ A-B vs homotypic (mutant − N2)",
               `A to A` = "Δ preference A-A vs A-B (mutant − N2)")
  qty_of  <- function(m_aa, m_bb, m_ab)
    list(`B to B` = m_bb - m_ab, `A to A` = m_aa - m_ab, `A to B` = m_ab - (m_aa + m_bb) / 2)

  ## source data: medians per chromosome and pooled, and the three quantities
  med_dt <- rbind(dat_f[, .(med = median(value), n = .N), by = .(xp, chr, class)],
                  dat_f[, .(med = median(value), n = .N), by = .(xp, class)][, chr := "all"])
  med_dt[, cl := c(`A to A` = "AA", `B to B` = "BB", `A to B` = "AB")[class]]
  w <- dcast(med_dt, xp + chr ~ cl, value.var = c("med", "n"))
  w[, `:=`(pref_BB = med_BB - med_AB, pref_AA = med_AA - med_AB,
           segr_AB = med_AB - (med_AA + med_BB) / 2)]
  ref_w <- w[xp == ref_xp, .(chr, r_BB = pref_BB, r_AA = pref_AA, r_AB = segr_AB)]
  w <- ref_w[w, on = "chr"][, `:=`(delta_BB = pref_BB - r_BB, delta_AA = pref_AA - r_AA,
                                   delta_AB = segr_AB - r_AB)][, c("r_BB", "r_AA", "r_AB") := NULL]
  xp_lev  <- intersect(gp_xp, unique(dat_f$xp))
  mut_lev <- setdiff(xp_lev, ref_xp)
  w[, xp := factor(xp, levels = xp_lev)]; setorder(w, xp, chr)
  fwrite(w, file.path(out_d, sprintf("Preference.%s.csv", chrom_idx)))

  ## m-out-of-n bootstrap of the within-map quantities
  vals  <- lapply(setNames(xp_lev, xp_lev), function(x)
             lapply(setNames(three, three), function(cl) dat_f[xp == x & class == cl, value]))
  m_eff <- min(n_sub, min(unlist(lapply(vals, lengths))))
  set.seed(boot_seed)
  boot_q <- lapply(vals, function(vl) {
    med <- lapply(vl, function(v)
      colMedians(matrix(sample(v, m_eff * n_boot, replace = TRUE), nrow = m_eff)))
    qty_of(med[["A to A"]], med[["B to B"]], med[["A to B"]])
  })

  obs     <- w[chr == "all"]
  dist_dt <- rbindlist(lapply(mut_lev, function(x) rbindlist(lapply(cls3, function(k)
               data.table(xp = x, class = k, delta = boot_q[[x]][[k]] - boot_q[[ref_xp]][[k]])))))
  bs <- dist_dt[, .(lo = quantile(delta, .025, names = FALSE),
                    hi = quantile(delta, .975, names = FALSE),
                    p  = max(2 * min(mean(delta <= 0), mean(delta >= 0)), 1 / n_boot)),
                by = .(xp, class)]
  bs[, delta_obs := mapply(function(x, k) obs[xp == x][[dcol_of[[k]]]], xp, class)]
  bs[, label := as.character(cut(p, breaks = c(-Inf, 1e-3, 1e-2, 5e-2, Inf),
                                 labels = c("***", "**", "*", "ns")))]
  bs[, `:=`(m = m_eff, B = n_boot)]
  fwrite(bs[, .(xp, class, delta_obs, lo, hi, p, label, m, B)],
         file.path(out_d, sprintf("Preference_Bootstrap.%s.csv", chrom_idx)))

  ## group contrast: the two mutants with impaired compartmentalization against
  ## wild type and the single mutants, on the same bootstrap draws
  imp_lev  <- grep("hpl2-lin61|met2-set25-set32", xp_lev, value = TRUE)
  comp_lev <- setdiff(xp_lev, imp_lev)
  obs_col  <- c(`B to B` = "pref_BB", `A to A` = "pref_AA", `A to B` = "segr_AB")
  grp <- rbindlist(lapply(cls3, function(k) {
    d <- Reduce(`+`, lapply(imp_lev,  function(x) boot_q[[x]][[k]])) / length(imp_lev) -
         Reduce(`+`, lapply(comp_lev, function(x) boot_q[[x]][[k]])) / length(comp_lev)
    o <- mean(sapply(imp_lev,  function(x) obs[xp == x][[obs_col[[k]]]])) -
         mean(sapply(comp_lev, function(x) obs[xp == x][[obs_col[[k]]]]))
    data.table(class = k, group_1 = paste(imp_lev, collapse = ";"),
               group_2 = paste(comp_lev, collapse = ";"), delta_obs = o,
               lo = quantile(d, .025, names = FALSE), hi = quantile(d, .975, names = FALSE),
               p  = max(2 * min(mean(d <= 0), mean(d >= 0)), 1 / n_boot), m = m_eff, B = n_boot)
  }))
  grp[, label := as.character(cut(p, breaks = c(-Inf, 1e-3, 1e-2, 5e-2, Inf),
                                  labels = c("***", "**", "*", "ns")))]
  fwrite(grp, file.path(out_d, sprintf("Preference_Group.%s.csv", chrom_idx)))

  ## panels: one point per mutant with its bootstrap interval, genotypes in the
  ## order of the figure, and the group contrast as a bracket over the two
  ## mutants with impaired compartmentalization
  pal <- setNames(RColorBrewer::brewer.pal(max(3, length(mut_lev)), "Set1")[seq_along(mut_lev)], mut_lev)
  bs[, xp := factor(xp, levels = mut_lev)]
  yr_pt   <- range(c(bs$lo, bs$hi, 0))
  span    <- diff(yr_pt)
  fig_lev <- c(setdiff(mut_lev, imp_lev), intersect(imp_lev[order(!grepl("hpl2-lin61", imp_lev))], mut_lev))

  for (k in cls3) {
    sub <- copy(bs[class == k])[, xp := factor(as.character(xp), levels = fig_lev)]
    sub[, ystar := yr_pt[2] + 0.08 * span]
    g  <- grp[class == k]
    xs <- sort(c(match(imp_lev[1], fig_lev), match(imp_lev[2], fig_lev)))
    yb <- yr_pt[2] + 0.20 * span

    gg[[sprintf("Boxplot.PreferenceGroup.%s.%s", k, chrom_idx)]] <-
      ggplot(sub, aes(x = xp, y = delta_obs, colour = xp)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = .18, linewidth = .6) +
      geom_point(size = 2.8) +
      geom_text(aes(y = ystar, label = label), colour = "black", size = 3.2) +
      annotate("segment", x = xs[1], xend = xs[2], y = yb, yend = yb, linewidth = .4) +
      annotate("segment", x = xs, xend = xs, y = yb, yend = yb - 0.03 * span, linewidth = .4) +
      annotate("text", x = mean(xs), y = yb + 0.10 * span, size = 2.9, lineheight = .9,
               label = sprintf("vs N2 and single mutants\n%+.2f [%+.2f, %+.2f] %s",
                               g$delta_obs, g$lo, g$hi, g$label)) +
      scale_colour_manual(values = pal, guide = "none") +
      scale_x_discrete(drop = FALSE) +
      coord_cartesian(ylim = c(yr_pt[1], yr_pt[2] + 0.36 * span), clip = "off") +
      labs(x = NULL, y = ylab_of[[k]],
           title = sprintf("%s (%s) - within-map quantity, m-out-of-n bootstrap (m=%d, B=%d)",
                           k, chrom_idx, m_eff, n_boot),
           subtitle = "stars: mutant vs N2; bracket: mean of the two mutants vs mean of the other genotypes") +
      theme_classic2() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
}

# SAVE ----
for (nm in names(gg))
  ggsave(file.path(out_d, paste0(nm, ".png")), gg[[nm]], width = 6, height = 5, dpi = 300)
