# HEADER ====================================================== =
# GEN2I
# Supplementary Figure S7H
# Where the repeat families of Supplementary Fig. S7F sit in the wild-type genome
# ============================================================= =
#
# Answers the objection that the derepression of MULE-MuDR and CMC-Mirage in the
# triple mutant could CAUSE decompartmentalization: if these families acted on
# compartmentalization, their copies should be located in the regions concerned.
#
# EVERYTHING IS MEASURED IN THE WILD TYPE, and that is the methodological point.
# Comparing the copies to "the regions that lose B most" in the double and
# triple mutants would be meaningless, since those mutants no longer have
# compartments. The two features tested are therefore properties of N2: the
# wild-type B compartment and the triple peaks (bins overlapping lin-61, hpl-2
# and H3K9me2 peaks at once).
#
# THE COMPARATOR IS THE OTHER REPEAT COPIES, NOT THE GENOME. Repeats are
# concentrated on the chromosome arms, which are B-rich; a genome-wide
# comparison would call any family enriched.
#
# The published panel is the first scope, "B compartment (N2)".
#
# As for Supplementary Fig. S7F/S7G, the repeat annotation is read from the
# tables produced upstream; nothing is re-mapped or re-counted here.
# ============================================================= =


# WORKING DIRECTORY ----
# Run the scripts from the root of this repository (where manuscript_lib.R sits):
# setwd("/path/to/Compartments-Gen2i")

# SOURCE ----
source("manuscript_lib.R")

# PARAMETERS ----
out_d <- "output/figure_S7H"
dir.create(out_d, showWarnings = FALSE, recursive = TRUE)

bed_p   <- "data/RNAseq/rmsk_ce11.bed"
anno_p  <- "data/RNAseq/rmsk_ce11_annotations.txt"
genic_p <- "data/RNAseq/rmsk_ce11_genic.txt"

# The three families of Supplementary Fig. S7F.
fam_v <- c("MULE-MuDR","CMC-Mirage","Pao")

comp_value <- "B"
drop_class <- c("rRNA","ARTEFACT")
drop_tentative <- TRUE

# SETSEED ----
set.seed(123)

# RUN ----
## 1. BINS OF THE WILD-TYPE GENOME, 10 kb ----
### THE 25 kb COMPARTMENT CALL ----
# B is the strict definition (quantile groups 1-12 of the 50 quantiles of the
# wild-type PC1), not the sign of the eigenvector: the ambiguous middle band is
# set aside. It is the definition already used by the enrichment analyses of the
# project. This panel uses the merged-batch PC1, as the pipeline does.
Go25 <- loadranges("data/ChIPseq/ce11_tiled_25kb.bed", genome = "ce11")

params.aname = "eigen_pca1_N2"
params = list(
  bw_p = 'data/HiC/N2_merged.bwa_mem.25kb.pca1.bw',
  genome = 'ce11',
  desc = 'Bigwig signal of eigen_pca1 N2'
)
Go25 <- bw_tile(bw = params$bw_p,
                bed = Go25,
                genome = config.genome,
                name = params.aname)

params.aname = "eigen_pca1_N2_50tile"
params.dependencies = "eigen_pca1_N2"
params = list(ngroups = 50)
params$ngroups <- as.integer(params$ngroups)
mcols(Go25)[params.aname] <- NA
mcols(Go25)[params.aname] <- dplyr::ntile(unlist(mcols(Go25)[[params.dependencies]]),
                                          params$ngroups)

params.aname = "eigen_pca1_N2_grps_shifted_A_B"
params.dependencies = "eigen_pca1_N2_50tile"
params = list(values = c(12, 16), names = c("B", "na", "A"))
mcols(Go25)[params.aname] <- NA
# convert NA to 0
f <- is.na(mcols(Go25)[[params.dependencies]])
mcols(Go25)[[params.dependencies]][f] <- 0
for(i in 1:length(params$values)){
  if(i == 1){
    f <- as.vector(mcols(Go25)[[params.dependencies]] <= params$values[i])
    if(any(f))
      mcols(Go25[f])[params.aname] <- params$names[i]
  }
  else{
    f <- as.vector(mcols(Go25)[[params.dependencies]]>params$values[i-1] & mcols(Go25)[[params.dependencies]]<=params$values[i])
    if(any(f))
      mcols(Go25[f])[params.aname] <- params$names[i]
  }
}
f <- as.vector(mcols(Go25)[[params.dependencies]]>params$values[i])
if(any(f))
  mcols(Go25[f])[params.aname] <- params$names[i+1]

### THE 10 kb BINS, WHICH CARRY THE PEAK-OVERLAP ANNOTATION ----
Go <- loadranges("data/ChIPseq/ce11_tiled_10kb.bed", genome = "ce11")

params.aname = "eigen_pca1_N2_grps_shifted_A_B"
params = list(aname = 'eigen_pca1_N2_grps_shifted_A_B')
idx <- GenomicRanges::nearest(Go, Go25)
mcols(Go)[params.aname] <- NA
Go_na_idx <- which(is.na(idx))
if (length(Go_na_idx)) {
  idx <- idx[!is.na(idx)]
  mcols(Go[-Go_na_idx])[params.aname] <- mcols(Go25[idx])[params$aname]
} else {
  mcols(Go)[params.aname] <- mcols(Go25[idx])[params$aname]
}

### TRIPLE PEAKS : lin-61 x hpl-2 x H3K9me2 ----
# The landmark is NOT "the H3K9me2 peaks": it is the intersection of the lin-61,
# hpl-2 and H3K9me2 peaks, the triple peaks of the manuscript.
for (an in c("lin61_hpl2_intersect_WF_narrowPeaks", "H3K9me2_WF_narrowPeak")) {
  bed_i <- switch(an,
                  lin61_hpl2_intersect_WF_narrowPeaks = "data/ChIPseq/hpl2_lin61_N2_sharp_intersect_peakset.bed",
                  H3K9me2_WF_narrowPeak               = "data/ChIPseq/H3K9me2_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed")
  Go2 <- loadranges(bed_i, genome = "ce11")
  f <- GenomicRanges::findOverlaps(Go, Go2)@from
  mcols(Go)[an] <- F
  mcols(Go[f])[an] <- T
}

params.aname = "lin61_hpl2_H3K9me2_WF_bin_ovlp"
params.dependencies = c("lin61_hpl2_intersect_WF_narrowPeaks", "H3K9me2_WF_narrowPeak")
mcols(Go)[params.aname] <- NA
mcols(Go)[params.aname] <- Gooperation(unlist(mcols(Go)[[params.dependencies[1]]]),
                                       unlist(mcols(Go)[[params.dependencies[2]]]),
                                       method = "and")

comp_col <- "eigen_pca1_N2_grps_shifted_A_B"
peak_col <- "lin61_hpl2_H3K9me2_WF_bin_ovlp"

bin_dt <- go2dt(Go)
bin_dt <- bin_dt[, .(chr = as.character(seqnames), bs = as.integer(start), be = as.integer(end),
                     comp = as.character(get(comp_col)),
                     peak = as.logical(get(peak_col)))]
bin_dt <- bin_dt[!is.na(peak)]
bin_dt[, is_B := comp %in% comp_value]

cat(sprintf("  %d bins | strict B : %d (%.1f %%) | triple peaks : %d (%.1f %%)\n",
            nrow(bin_dt), sum(bin_dt$is_B), 100 * mean(bin_dt$is_B),
            sum(bin_dt$peak), 100 * mean(bin_dt$peak)))

## 2. REPEAT COPIES - same scope as Supplementary Fig. S7F/S7G ----
bed <- fread(bed_p, col.names = c("chr","start","end","name","score","strand"))
ann <- merge(fread(anno_p, select = c("name","repClass","repFamily")),
             fread(genic_p), by = "name")
ann[, genic_ovl := as.logical(genic_ovl)]

cp <- merge(bed[, .(name, chr, start, end)], ann, by = "name")
cp <- cp[!repClass %in% drop_class & !genic_ovl]
if (drop_tentative)
  cp <- cp[!grepl("\\?$", repFamily) & !grepl("\\?$", repClass)]
# A copy is attached to ONE bin, by its midpoint: a repeat copy is much shorter
# than a bin, splitting it would bring nothing.
cp[, mid := (start + end) %/% 2L]

setkey(bin_dt, chr, bs, be)
# The position columns carry explicit names: `s` and `e` would collide with the
# local variables of the test block, and data.table would silently resolve the
# COLUMN instead of the variable.
cpb <- foverlaps(cp[, .(name, repFamily, chr, pos_s = mid, pos_e = mid)], bin_dt,
                 by.x = c("chr","pos_s","pos_e"), by.y = c("chr","bs","be"),
                 type = "within", nomatch = NULL)
cat(sprintf("  %d intergenic copies attached to a bin\n", nrow(cpb)))

miss <- setdiff(fam_v, unique(cpb$repFamily))
if (length(miss)) warning("Families absent from the copies kept: ",
                          paste(miss, collapse = ", "), call. = FALSE, immediate. = TRUE)
fam_v <- intersect(fam_v, unique(cpb$repFamily))

## 3. ENRICHMENT TESTS ----
# Fisher exact on a 2x2: (family vs other copies) x (inside the scope or not).
# The odds ratio reads directly: > 1 = the family is over-represented in the
# scope RELATIVE TO THE OTHER REPEATS.
#
# The third scope is the most informative: restricted to the B compartment, it
# asks whether the peaks add anything beyond the compartment that carries them.
SCOPES <- list(
  `B compartment (N2)`                     = list(sub = function(d) d,               hit = "is_B"),
  `Triple peaks: lin-61 x hpl-2 x H3K9me2` = list(sub = function(d) d,               hit = "peak"),
  `Triple peaks, within B`                 = list(sub = function(d) d[is_B == TRUE], hit = "peak"))

res <- rbindlist(lapply(names(SCOPES), function(sc) {
  sco <- SCOPES[[sc]]
  d <- copy(sco$sub(cpb))
  if (!nrow(d)) return(NULL)
  d[, hit := as.logical(d[[sco$hit]])]
  rbindlist(lapply(fam_v, function(f) {
    a <- d[repFamily == f, hit]; b <- d[repFamily != f, hit]
    if (length(a) < 5) return(NULL)
    ft <- fisher.test(matrix(c(sum(a), length(a) - sum(a),
                               sum(b), length(b) - sum(b)), nrow = 2))
    data.table(scope = sc, family = f,
               n_copies = length(a), n_hit = sum(a),
               pct_family = 100 * mean(a), pct_others = 100 * mean(b),
               odds_ratio = unname(ft$estimate),
               ci_low = ft$conf.int[1], ci_high = ft$conf.int[2],
               p_value = ft$p.value)
  }))
}))

# BH over all the tests of the panel: 3 families x 3 scopes form one family of
# tests.
res[, padj := p.adjust(p_value, method = "BH")]
res[, stars := fifelse(is.na(padj), "",
             fifelse(padj < 0.001, "***",
             fifelse(padj < 0.01,  "**",
             fifelse(padj < 0.05,  "*", "ns"))))]
print(res[, .(scope, family, n_copies, pct_family = round(pct_family,1),
              pct_others = round(pct_others,1), OR = round(odds_ratio,2),
              padj = signif(padj,2), stars)])

fwrite(res, file.path(out_d, "Positional_enrichment.csv"))

## 4. FIGURES ----
theme_g2i <- theme_minimal(base_size = 10) +
  theme(plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_blank(),
        axis.text  = element_text(face = "bold", size = 9),
        strip.text = element_text(face = "bold", size = 9),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, colour = "#555555"),
        plot.caption  = element_text(hjust = 0, size = 7.5, colour = "#555555"))

forest_of <- function(pl, ttl) {
  # Fisher confidence intervals reach 0 or infinity when a cell is empty: bound
  # them for display rather than let the bar disappear.
  lo <- min(c(pl$ci_low[pl$ci_low > 0], pl$odds_ratio), na.rm = TRUE) / 1.6
  hi <- max(c(pl$ci_high[is.finite(pl$ci_high)], pl$odds_ratio), na.rm = TRUE) * 1.6
  pl <- copy(pl)[, `:=`(ci_lo_d = pmax(ci_low, lo), ci_hi_d = pmin(ci_high, hi))]
  ggplot(pl, aes(x = odds_ratio, y = family, colour = odds_ratio > 1)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "#CC3333", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = ci_lo_d, xmax = ci_hi_d), height = 0.16, linewidth = 0.7) +
    geom_point(size = 3) +
    # The label is set to the right of the bar rather than above the point:
    # above, it overlaps the confidence interval from the first family on.
    geom_text(aes(x = ci_hi_d, label = paste0("  ", sprintf("%.0f", pct_family), " %  ", stars)),
              hjust = 0, vjust = 0.5, size = 3, show.legend = FALSE) +
    scale_x_continuous(trans = "log2", limits = c(lo, hi * 2.2)) +
    scale_colour_manual(values = c(`TRUE` = "#1B7837", `FALSE` = "#762A83"), guide = "none") +
    facet_wrap(~ scope, ncol = 3) +
    labs(title = ttl,
         subtitle = paste0("Odds ratio vs all other intergenic repeat copies | Fisher exact, ",
                           "BH-adjusted | label: % of the family's copies, then significance\n",
                           "Triple peaks = bins overlapping lin-61, hpl-2 and H3K9me2 peaks at once"),
         x = "Odds ratio (log2 scale)", y = NULL,
         caption = paste0("Everything is measured in N2. A ratio of 1 means the family is ",
                          "distributed like the other repeats.\nThe comparator is the other repeat ",
                          "copies, not the genome: repeats are concentrated on the arms, which are ",
                          "B-rich,\nso a genome-wide comparison would call any family enriched.")) +
    theme_g2i +
    theme(aspect.ratio = 0.42, panel.grid.major.y = element_blank())
}

pl <- copy(res)
pl[, scope  := factor(scope, levels = names(SCOPES))]
pl[, family := factor(family, levels = rev(fam_v))]

gg <- list()
# Forest plot: the test itself, read directly. A confidence interval crossing 1
# is a family distributed like the other repeats.
gg[["Enrichment.Forest"]] <- forest_of(pl, "Where the repeat families sit, in the wild-type genome")
# The published panel is the first scope alone.
gg[["Enrichment.Forest.B_compartment"]] <-
  forest_of(pl[scope == "B compartment (N2)"],
            "Positional enrichment in the wild-type B compartment")

# Bars: the same result as raw proportions, for whoever prefers shares to odds
# ratios.
bar_dt <- rbindlist(list(
  pl[, .(scope, family = as.character(family), grp = "this family", pct = pct_family)],
  unique(pl[, .(scope, grp = "all other repeats", pct = pct_others)])[
    , .(scope, family = "all other repeats", grp, pct)]))
bar_dt[, family := factor(family, levels = c(fam_v, "all other repeats"))]

gg[["Enrichment.Proportions"]] <-
  ggplot(bar_dt, aes(x = family, y = pct, fill = grp)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.0f %%", pct)), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c(`this family` = "#1B9E77",
                               `all other repeats` = "#999999"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  facet_wrap(~ scope, ncol = 3) +
  labs(title = "Same result as proportions",
       subtitle = paste0("Share of each family's intergenic copies falling in the given wild-type ",
                         "feature\nTriple peaks = bins overlapping lin-61, hpl-2 and H3K9me2 peaks at once"),
       x = NULL, y = "% of copies") +
  theme_g2i +
  theme(aspect.ratio = 0.9,
        axis.text.x = element_text(angle = 30, hjust = 1))

# SAVE ----
ggsave(file.path(out_d, "Enrichment.Forest.png"),               gg[["Enrichment.Forest"]],               width = 11, height = 4, dpi = 150)
ggsave(file.path(out_d, "Enrichment.Forest.B_compartment.png"), gg[["Enrichment.Forest.B_compartment"]], width = 5.5, height = 4, dpi = 150)
ggsave(file.path(out_d, "Enrichment.Proportions.png"),          gg[["Enrichment.Proportions"]],          width = 11, height = 5, dpi = 150)
