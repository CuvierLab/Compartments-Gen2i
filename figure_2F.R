# HEADER ====================================================== =
# GEN2I
# Tue Jun 18 18:22:27 2024
# R version : R version 4.4.0 (2024-04-24)
# System : x86_64, linux-gnu
# ============================================================= =


# WORKING DIRECTORY ----
setwd(dir = "config/src/R/github/")

# SOURCE ----
source("manuscript_lib.R")

# RUN ----
## RNA SEQ COUNTS RATIO ----
### MCOLS Chr bin ----
Go <- loadranges("data/ChIPseq/ce11_tiled_25kb.bed",genome = "ce11")
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

### GENES COUNTS ----
wb_genes <- loadranges("data/ChIPseq/genes_ce11.bed",genome = "ce11")

for (condition in c("hpl_2I158A","lin61","N2_WT","lin61_hpl2","met2_set25_set32","hpl2")) {
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
    # save.image("join.RData")
    mcols(wb_genes)[params.aname] <- NA
    mcols(wb_genes)[params.aname] <-  dt[wb_genes$name, .SD, .SDcols = params$col]
  }
}

#### DO THE MEAN ----
for (condition in c("hpl_2I158A","lin61","N2_WT","lin61_hpl2","met2_set25_set32","hpl2")) {
  
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

### INHERIT EIGEN PCA1 N2 ON GENES ----
params.aname = "eigen_pca1_N2.old"
params = list(
  aname = 'eigen_pca1_N2.old'
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


### PLOT ----
params = list(
  y = c('N2_WT.mean.cnts', 'hpl2.mean.cnts', 'hpl_2I158A.mean.cnts', 'lin61.mean.cnts', 'lin61_hpl2.mean.cnts', 'met2_set25_set32.mean.cnts'),
  eigenMcol = 'eigen_pca1_N2.old',
  legend = TRUE,
  height = 800,
  width = 800
)

y <- params$y
y_nm <- "log2 Ratio( A counts / B counts )"
eigenMcol <- params$eigenMcol

wb_genes_dt <- go2dt(wb_genes)
wb_genes_dt[,comp:=ifelse(get(eigenMcol) > 0, "A", "B")]

plot_dt <- NULL
for (i_smpl in 1:50) {
  wb_genes_spl_dt <- wb_genes_dt[sample(.N, 10000)]
  ratio_dt <- wb_genes_spl_dt[,lapply(.SD,sum),  .SDcols = y, by = comp]
  plot_dt <- rbind(plot_dt,setnames(melt(log2(ratio_dt[comp=="A",.SD,.SDcols = y]/ratio_dt[comp=="B",.SD,.SDcols = y])),c("condition","ratio")))
}

plot_dt$condition <- as.character(plot_dt$condition)

myComp <- combn(x = unique(plot_dt$condition), m = 2, simplify = F)
f0 <- unlist(lapply(myComp, function(X) { any(grepl("N2", X))}))
f1 <- unlist(lapply(myComp, function(X) { any(grepl("^lin61\\.", X)) & any(grepl("lin61_hpl2",X))}))
f2 <- unlist(lapply(myComp, function(X) { any(grepl("^hpl2\\.", X)) & any(grepl("lin61_hpl2",X))}))
myComp <- myComp[f0 | f1 | f2]

old_options <- options(ggpubr.parse_aes = FALSE)
gg <- list()
gg[["Boxplot.Ratio_all_conditions"]] <- ggboxplot(plot_dt, x = "condition", y = "ratio", fill = "condition") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  stat_compare_means(aes(label = ..p.signif..),
                     comparisons = myComp,
                     paired = T) + ylab(y_nm)
