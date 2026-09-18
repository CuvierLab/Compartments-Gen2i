# HEADER ====================================================== =
# GEN2I
# Supplementary Figure S7E
# RNA-seq levels per genotype, in active and inactive TAD-like structures
# ============================================================= =


# WORKING DIRECTORY ----
# Run the scripts from the root of this repository (where manuscript_lib.R sits):
# setwd("/path/to/Compartments-Gen2i")

# SOURCE ----
source("manuscript_lib.R")

# PARAMETERS ----
out_d <- "output/figure_S7E"
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

# The seven wild-type peak sets whose density over the TADs is clustered.
# The order is the one of the pipeline, and it sets the row order of the
# correspondence analysis.
ip_p <- c(H3K4me3_WF_narrowPeak_cov  = "data/ChIPseq/H3K4me3_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
          H3K9me3_WF_narrowPeak_cov  = "data/ChIPseq/H3K9me3_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
          hpl2_WF_narrowPeak_cov     = "data/ChIPseq/hpl2_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
          H3K9me2_WF_narrowPeak_cov  = "data/ChIPseq/H3K9me2_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
          H3K27me3_WF_narrowPeak_cov = "data/ChIPseq/H3K27me3_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
          LEM2_WF_narrowPeak_cov     = "data/ChIPseq/LEM2_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
          lin61_WF_narrowPeak_cov    = "data/ChIPseq/lin61_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed")

# SETSEED ----
set.seed(123)

# RUN ----
## TAD CLASSES FROM THE DENSITY OF THE WILD-TYPE CHIP PEAKS ----
### TADs OF THE INITIAL BATCH, 25 kb ----
Go <- loadranges("data/HiC/N2-old_merged.bwa_mem.25kb_domains.bed", genome = "ce11")

### PEAK DENSITY OVER EACH TAD, CHROMOSOME BY CHROMOSOME ----
for (params.aname in names(ip_p)) {
  Go2 <- loadranges(ip_p[[params.aname]], genome = "ce11")
  params = list(bychr = T)
  if (is.null(params$bychr)) {
    params$bychr <- T
  }
  mcols(Go)[[params.aname]] <- getCoverageDensity(Go, Go2, bychr = params$bychr)
}

### CORRESPONDENCE ANALYSIS CLUSTERING ----
# Verbatim from the pipeline annotation script
# resources/scripts/db/anno/set_tad_class_ca_clustering.R, which is what fills the
# tad_classes_chipdensity annotation of tad_N2.old_25kb.
params.aname <- "tad_classes_chipdensity"
params.rname <- "tad_N2.old_25kb"

# chromosomes
chr_v <- c("I","II","III","IV","V","X")

go_chr <- subset(Go,seqnames %in% chr_v)
go_df <- as.data.frame(go_chr)
toclust_bin_df <- subset(t(go_df),grepl("cov",colnames(go_df)))
toclust_bin_df <- t(apply(toclust_bin_df,1,as.numeric))
colnames(toclust_bin_df) <- 1:ncol(toclust_bin_df)

# caCluster draws on the active device; keep those diagnostics in output/.
pdf(file = file.path(out_d, "TAD_clustering.CAcluster_diagnostics.pdf"))
res <- caCluster(toclust_bin_df, which="both", dim=2, opt.part=TRUE)
resC <- caCluster(toclust_bin_df, which="cols", dim=2, opt.part=TRUE)
resCfix <- caCluster(toclust_bin_df,part = 4, which="cols", dim=2, opt.part=TRUE)
resR <- caCluster(toclust_bin_df, which="rows", dim=2, opt.part=TRUE)
dev.off()

names(resR) <- paste0("class_",1:length(resR))
names(resC) <- paste0("class_",1:length(resC))
names(resCfix) <- paste0("class_",1:length(resCfix))
col_ord <- plyr::ldply(resC,.id = "tad_class",function(X)data.frame(tad_idx_chr=names(X)))
col_ord_fix <- plyr::ldply(resCfix,.id = "tad_class_fix",function(X)data.frame(tad_idx_chr=names(X)))
row_ord <- plyr::ldply(resR,.id = "chip_class",function(X)data.frame(chip_ip=names(X)))

toclust_bin_df <- toclust_bin_df[match(row_ord$chip_ip,rownames(toclust_bin_df)),]
toclust_bin_df <- toclust_bin_df[,match(col_ord$tad_idx_chr,colnames(toclust_bin_df))]

pdf(file = file.path(out_d, "TAD_clustering.Caclustering_heatmap.pdf"))
print(Heatmap(scale(toclust_bin_df), name = "mat",
              cluster_rows = F,cluster_columns = F,
              border="black",column_title = "CA clustering",
              row_split = row_ord$chip_class,
              column_split = col_ord$tad_class))
dev.off()

mcols(Go)[params.aname] <- NA
if(grepl("old",params.rname)){
  GenomicRanges::mcols(Go[as.numeric(col_ord$tad_idx_chr)])[params.aname] <- ifelse(col_ord$tad_class=="class_1","Active","Inactive")
} else {
  GenomicRanges::mcols(Go[as.numeric(col_ord$tad_idx_chr)])[params.aname] <- ifelse(col_ord$tad_class=="class_2","Active","Inactive")
}

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

### INHERIT THE TAD CLASS ON GENES ----
params.aname = "N2-old_tad_classes_chipdensity"
params = list(
  aname = 'tad_classes_chipdensity'
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
# Same construction as Supplementary Fig. S7D: the published panel keeps the
# six-genotype colour scale of the pipeline and simply drops the hpl-2(I158A)
# box, so the five remaining boxes keep the fill they had with six levels.
pal_6 <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")

gg <- list()
stat_dt <- NULL
box_dt <- NULL
tad_cls <- c("Active", "Inactive")

for (rep_nm in rep_v) {

  y_all <- cond_all_v %+% "." %+% rep_nm %+% ".cnts"
  y_pl  <- cond_v     %+% "." %+% rep_nm %+% ".cnts"
  pal_pl <- pal_6[match(y_pl, y_all)]
  ctrl_nm <- "N2_WT." %+% rep_nm %+% ".cnts"

  for (cls in tad_cls) {

    # filter : N2-old_tad_classes_chipdensity == cls
    Go_f <- subset(wb_genes,
                   unlist(mcols(wb_genes)[["N2-old_tad_classes_chipdensity"]]) == cls)

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
                                   .(tad_class = cls, replicate = rep_nm, n = .N,
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
      stat_dt <- rbind(stat_dt, data.table(tad_class = cls, replicate = rep_nm,
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

    nm_base <- "Boxplot.rmOutliers.TAD_class_" %+% cls %+% ".Replicate_" %+% rep_nm

    gg[[nm_base %+% ".with_pval"]] <-
      ggboxplot(pl_dt, x = "sig_mcol", y = "score", fill = "sig_mcol",
                palette = pal_pl) +
      stat_compare_means(aes(label = ..p.signif..),
                         comparisons = myComp,
                         paired = FALSE,
                         method = "wilcox.test",
                         method.args = list(alternative = "less")) +
      labs(x = cls %+% " TAD-like structures", y = "Score [arb. units]") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    gg[[nm_base %+% ".no_pval"]] <-
      ggboxplot(pl_dt, x = "sig_mcol", y = "score", fill = "sig_mcol",
                palette = pal_pl) +
      labs(x = cls %+% " TAD-like structures", y = "Score [arb. units]") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
}

# SAVE ----
tad_dt <- go2dt(Go)
fwrite(tad_dt,  file.path(out_d, "Boxplot.TAD_classes_chipdensity.csv"))
fwrite(box_dt,  file.path(out_d, "Boxplot.TAD_class_boxstats.csv"))
fwrite(stat_dt, file.path(out_d, "Boxplot.TAD_class_wilcoxon.csv"))

for (nm in names(gg))
  ggsave(file.path(out_d, paste0(nm, ".png")), gg[[nm]], width = 6, height = 6, dpi = 150)
