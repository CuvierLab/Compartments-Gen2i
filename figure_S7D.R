# HEADER ====================================================== =
# GEN2I
# Supplementary Figure S7D
# RNA-seq levels per genotype, in compartment A and in compartment B
# ============================================================= =


# WORKING DIRECTORY ----
# Run the scripts from the root of this repository (where manuscript_lib.R sits):
# setwd("/path/to/Compartments-Gen2i")

# SOURCE ----
source("manuscript_lib.R")

# PARAMETERS ----
out_d <- "output/figure_S7D"
dir.create(out_d, showWarnings = FALSE, recursive = TRUE)

# Conditions, in the order used for the published panel. hpl_2I158A is present in
# the count table but is not shown in any figure of the manuscript; it is kept in
# the source table below and left out of the plot only.
cond_v <- c("N2_WT", "hpl2", "lin61", "lin61_hpl2", "met2_set25_set32")
cond_all_v <- c("N2_WT", "hpl2", "hpl_2I158A", "lin61", "lin61_hpl2", "met2_set25_set32")
# The published panel is the rep1 one; the pipeline also produces rep2 and the
# average of the two replicates, so all three are kept here.
rep_v <- c("rep1", "rep2", "mean")
rep_pub <- "rep1"
# Low counts are discarded, as in the pipeline.
min_cnt <- 5

# SETSEED ----
set.seed(123)

# RUN ----
## COMPARTMENTS ON THE 25 kb TILED GENOME ----
### MCOLS Chr bin ----
Go <- loadranges("data/ChIPseq/ce11_tiled_25kb.bed", genome = "ce11")
params.aname = "chr_bin"
params = list(bychr= T)
mcol <- NULL
if(params$bychr) {
  mcol <- unlist(lapply(split(Go,seqnames(Go)),function(chr){1:length(chr)}))
} else mcol <- 1:length(Go)
mcols(Go)[params.aname] <- NA
mcols(Go)[params.aname] <- mcol

### MCOLS Eigen values ----
params.aname = "eigen_pca1_N2.old"
params = list(
  bw_p = 'data/HiC/N2-old_merged.bwa_mem.25kb.pca1.bw',
  genome = 'ce11',
  desc = 'Bigwig signal of eigen_pca1 N2-old'
)

Go <- bw_tile(bw = params$bw_p,
              bed = Go,
              genome = config.genome,
              name = params.aname)

### TILE GENOME ----
params.aname = "eigen_pca1_N2.old_50tile"
params.process =
  params = list(
    ngroups = 50
  )
params.dependencies = "eigen_pca1_N2.old"
params$ngroups <- as.integer(params$ngroups)

mcols(Go)[params.aname] <- NA
mcols(Go)[params.aname] <- dplyr::ntile(unlist(mcols(Go)[[params.dependencies]]),
                                        params$ngroups)

### CUT INTO A / B GROUPS ----
# "The eigen groups correspond to quantiles : B = 1-12, A = 16-50", the wide
# definition used by this panel (annotation eigen_pca1_N2.old_grps_shifted_A_B).
# Groups 13-15 form the ambiguous middle band and are labelled "na".
params.aname = "eigen_pca1_N2.old_grps_shifted_A_B"
params.dependencies = "eigen_pca1_N2.old_50tile"
params = list(
  values = c(12, 16),
  names  = c("B", "na", "A")
)

mcols(Go)[params.aname] <- NA
# convert NA to 0
f <- is.na(mcols(Go)[[params.dependencies]])
mcols(Go)[[params.dependencies]][f] <- 0

for(i in 1:length(params$values)){
  if(i == 1){
    f <- as.vector(mcols(Go)[[params.dependencies]] <= params$values[i])
    if(any(f))
      mcols(Go[f])[params.aname] <- params$names[i]
  }
  else{
    f <- as.vector(mcols(Go)[[params.dependencies]]>params$values[i-1] & mcols(Go)[[params.dependencies]]<=params$values[i])
    if(any(f))
      mcols(Go[f])[params.aname] <- params$names[i]
  }
}
f <- as.vector(mcols(Go)[[params.dependencies]]>params$values[i])
if(any(f))
  mcols(Go[f])[params.aname] <- params$names[i+1]

## GENES ----
### GENES COUNTS ----
wb_genes <- loadranges("data/ChIPseq/genes_ce11.bed",genome = "ce11")

for (condition in cond_all_v) {
  for (rep in "rep" %+% 1:2) {
    params.aname = condition %+% "." %+% rep %+% ".cnts"
    params = list(
      path = 'data/RNAseq/deseq2_counts_all_conditions.txt',
      key = 'rn',
      col = condition %+% '.' %+% rep
    )

    dt <- fread(params$path)
    if(!is.null(params$key)){
      if(is.numeric(params$key))
        params$key <- colnames(dt)[params$key]
      setkeyv(dt, params$key)
    }else{
      if("V1" %in% colnames(dt))
        setkey(dt, V1)
      if("ID" %in% colnames(dt))
        setkey(dt, ID)
      if("name" %in% colnames(dt))
        setkey(dt, name)
    }
    mcols(wb_genes)[params.aname] <- NA
    mcols(wb_genes)[params.aname] <-  dt[wb_genes$name, .SD, .SDcols = params$col]
  }
}

#### DO THE MEAN ----
for (condition in cond_all_v) {

  params.aname = condition %+% ".mean.cnts"
  params.dependencies = c(condition %+% ".rep1.cnts",condition %+% ".rep2.cnts")
  params = list(
    operation = 'mean'
  )

  mcols(wb_genes)[params.aname] <- NA
  method <- ifelse(!is.null(params$method),params$method, params$operation)
  mcols(wb_genes)[params.aname] <- Gooperation(unlist(mcols(wb_genes)[[params.dependencies[1]]]),
                                               unlist(mcols(wb_genes)[[params.dependencies[2]]]),
                                               method = method)
}

### INHERIT THE COMPARTMENT CALL ON GENES ----
params.aname = "eigen_pca1_N2.old_grps_shifted_A_B"
params = list(
  aname = 'eigen_pca1_N2.old_grps_shifted_A_B'
)
idx <- GenomicRanges::nearest(wb_genes, Go)
# Initlalize query Mcol
mcols(wb_genes)[params.aname] <- NA
# if there is NA somewhere (mostly cause when using seqnames in one of the ranges and not the other)
# you should remove those indexes from both Go ranges
wb_genes_na_idx <- which(is.na(idx))
if (length(wb_genes_na_idx)) {
  idx <- idx[!is.na(idx)]
  mcols(wb_genes[-wb_genes_na_idx])[params.aname] <- mcols(Go[idx])[params$aname]
} else {
  mcols(wb_genes)[params.aname] <- mcols(Go[idx])[params$aname]
}

# If a distance threshold is set to retrieve nearest:
tokeep <- rep(T, length(wb_genes))
if (!is.null(params$distThreshold)) tokeep <- distanceToNearest(wb_genes, Go)@elementMetadata$distance <= params$distThreshold
mcols(wb_genes)[params.aname][!tokeep,] <- NA


## PLOT ----
# The published panel keeps the six-genotype colour scale of the pipeline and
# simply drops the hpl-2(I158A) box, so the five remaining boxes keep the fill
# they had with six levels. The levels are set explicitly rather than left to
# factor(), whose ordering depends on the collation locale.
pal_6 <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")

gg <- list()
stat_dt <- NULL
box_dt <- NULL

for (rep_nm in rep_v) {

  y_all <- cond_all_v %+% "." %+% rep_nm %+% ".cnts"
  y_pl  <- cond_v     %+% "." %+% rep_nm %+% ".cnts"
  pal_pl <- pal_6[match(y_pl, y_all)]
  # Control of every pairwise comparison.
  ctrl_nm <- "N2_WT." %+% rep_nm %+% ".cnts"

  for (comp in c("A", "B")) {

    # filter : eigen_pca1_N2.old_grps_shifted_A_B == comp
    Go_f <- subset(wb_genes,
                   unlist(mcols(wb_genes)[["eigen_pca1_N2.old_grps_shifted_A_B"]]) == comp)

    Go_dt <- go2dt(Go_f)
    Go_dt <- reshape(Go_dt,
                     varying = y_all,
                     v.names = "score",
                     timevar = "sig_mcol",
                     times = y_all,
                     direction = "long")

    # Filter low counts
    Go_dt <- Go_dt[eval(expression(eval(parse(text = paste0("score>", min_cnt))))),]

    # +/- OUTLIERS : outliers are set to NA within each genotype
    Gon_dt <- data.table::copy(Go_dt)
    Gon_dt[, score := remove_outliers(score), by = sig_mcol]

    # Source tables : every genotype, hpl_2I158A included
    box_dt <- rbind(box_dt, Gon_dt[!is.na(score),
                                   .(compartment = comp, replicate = rep_nm, n = .N,
                                     ymin = min(score), lower = quantile(score, .25),
                                     middle = median(score), upper = quantile(score, .75),
                                     ymax = max(score), mean = mean(score)),
                                   by = sig_mcol])

    for (nm in setdiff(y_all, ctrl_nm)) {
      # Same one-sided test as the panel: ggpubr passes the first member of the
      # comparison as x, and N2 comes first in the comparison list, so the
      # alternative "less" reads "N2 lower than the mutant".
      wt <- wilcox.test(Gon_dt[sig_mcol == ctrl_nm, score],
                        Gon_dt[sig_mcol == nm, score],
                        paired = FALSE, alternative = "less")
      stat_dt <- rbind(stat_dt, data.table(compartment = comp, replicate = rep_nm,
                                           group1 = ctrl_nm, group2 = nm,
                                           method = "wilcox.test",
                                           alternative = "less",
                                           p = wt$p.value,
                                           p.signif = symnum(wt$p.value,
                                                             cutpoints = c(0, 1e-4, 1e-3, 1e-2, 5e-2, 1),
                                                             symbols = c("****", "***", "**", "*", "ns"))))
    }

    # Plotted data : hpl_2I158A dropped
    pl_dt <- Gon_dt[sig_mcol %in% y_pl]
    pl_dt[, sig_mcol := factor(sig_mcol, levels = y_pl)]

    myComp <- combn(x = y_pl, m = 2, simplify = F)
    myComp <- myComp[unlist(lapply(myComp, function(x) ctrl_nm %in% x))]

    nm_base <- "Boxplot.rmOutliers.Compartment_" %+% comp %+% ".Replicate_" %+% rep_nm

    gg[[nm_base %+% ".with_pval"]] <-
      ggboxplot(pl_dt, x = "sig_mcol", y = "score", fill = "sig_mcol",
                palette = pal_pl) +
      stat_compare_means(aes(label = ..p.signif..),
                         comparisons = myComp,
                         paired = FALSE,
                         method = "wilcox.test",
                         method.args = list(alternative = "less")) +
      labs(x = comp %+% " compartment", y = "Score [arb. units]") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    gg[[nm_base %+% ".no_pval"]] <-
      ggboxplot(pl_dt, x = "sig_mcol", y = "score", fill = "sig_mcol",
                palette = pal_pl) +
      labs(x = comp %+% " compartment", y = "Score [arb. units]") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
}

# SAVE ----
fwrite(box_dt,  file.path(out_d, "Boxplot.Compartment_boxstats.csv"))
fwrite(stat_dt, file.path(out_d, "Boxplot.Compartment_wilcoxon.csv"))

for (nm in names(gg))
  ggsave(file.path(out_d, paste0(nm, ".png")), gg[[nm]], width = 6, height = 6, dpi = 150)
