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

### LOAD EIGEN----
## Get nearest TAD IDs at 5kb with TAD of N2 from WF
Go <- loadranges("data/ChIPseq/ce11_tiled_10kb.bed",genome = "ce11")
params.aname = "chr_bin"
params = list(bychr= T)
mcol <- NULL
if(params$bychr) {
  mcol <- unlist(lapply(split(Go,seqnames(Go)),function(chr){1:length(chr)}))
} else mcol <- 1:length(Go)
mcols(Go)[params.aname] <- NA
mcols(Go)[params.aname] <- mcol

Go2 <- loadranges("data/ChIPseq/ce11_tiled_25kb.bed",genome = "ce11")
params.aname = "eigen_pca1_N2.old"
params = list(
  bw_p = 'data/HiC/N2-old_merged.bwa_mem.25kb.pca1.bw',
  genome = 'ce11',
  desc = 'Bigwig signal of eigen_pca1 N2-old'
)
# convert to Granges

cat("bigwig_signal_tile:... ")
# save.image("dev_l.RData")

Go2 <- bw_tile(bw = params$bw_p,
               bed = Go2,
               genome = config.genome,
               name = params.aname)

params.aname = "eigen_pca1_N2.old"
params.process = 
  params = list(
    aname = 'eigen_pca1_N2.old'
  )

idx <- GenomicRanges::nearest(Go, Go2)
# Initlalize query Mcol
mcols(Go)[params.aname] <- NA
# if there is NA somewhere (mostly cause when using seqnames in one of the ranges and not the other)
# you should remove those indexes from both Go ranges
Go_na_idx <- which(is.na(idx))
if (length(Go_na_idx)) {
  idx <- idx[!is.na(idx)]
  mcols(Go[-Go_na_idx])[params.aname] <- unlist(mcols(Go2[idx])[params$aname])
} else {
  mcols(Go)[params.aname] <- unlist(mcols(Go2[idx])[params$aname])
}

# If a distance threshold is set to retrieve nearest:
tokeep <- rep(T, length(Go))
if (!is.null(params$distThreshold)) tokeep <- distanceToNearest(Go, Go2)@elementMetadata$distance <= params$distThreshold
mcols(Go)[params.aname][!tokeep,] <- NA




### GET TAD CLASS ON TAD N2 OLD ----
Go3 <- loadranges("data/HiC/N2-old_merged.bwa_mem.25kb_domains.bed",genome = "ce11")
params.aname = "chr_bin"
params = list(bychr= T)
mcol <- NULL
if(params$bychr) {
  mcol <- unlist(lapply(split(Go3,seqnames(Go3)),function(chr){1:length(chr)}))
} else mcol <- 1:length(Go3)
mcols(Go3)[params.aname] <- NA
mcols(Go3)[params.aname] <- mcol

### MCOLS COVERAGE ----
peaks_lst <- c("H3K4me3_GSE49739.bed",
               "H3K4me1_GSE50262.bed",
               "H3K9me3_GSE49732.bed",
               "hpl2_GSE100829.bed",
               "H3K9me2_GSE113841.bed",
               "H3K27ac_GSE49734.bed",
               "H3K27me3_GSE49738.bed",
               "LEM2_GSE24241.bed",
               "lin61_GSE49209.bed")

for (item in c(peaks_lst)) {
  
  params.aname = paste0(sub("(.*)_GSE.+","\\1",item),"_peaks_cov")
  params = list(
    bychr = TRUE,
    bed = paste0("data/ChIPseq/",item)
  )
  
  Go2 <- loadranges(params$bed, genome = "ce11")
  
  
  if (is.null(params$bychr)) {
    params$bychr <- T
  }
  
  # Go2 <- rangesDB2GO(filterRangesDB(rname = list(in_=list(value=params$rname, strict=TRUE))))
  mcols(Go3)[[params.aname]] <- getCoverageDensity(Go3, Go2, bychr=params$bychr)
}

params.aname = c('H3K4me3_peaks_cov', 'H3K4me1_peaks_cov', 'H3K9me3_peaks_cov', 'hpl2_peaks_cov', 'H3K9me2_peaks_cov', 'H3K27ac_peaks_cov', 'H3K27me3_peaks_cov', 'LEM2_peaks_cov', 'lin61_peaks_cov')
params = list(
  resolution = '25kb',
  height = 800,
  width = 800
)

output <- "results"
gp_res <- params$resolution
# chromosomes
chr_v <- c("I","II","III","IV","V","X")
# Graphical parameters
col_fun_sad <-  circlize::colorRamp2(seq(-1,1,length.out = 10), rev(brewer.pal(n = 10, name = "RdBu")))
col_fun_cov <-  circlize::colorRamp2(seq(0,1,length.out = 9), rev(brewer.pal(n = 9, name = "PuOr")))

go_chr <- subset(Go3,seqnames %in% chr_v)
go_df <- adf(go_chr)
toclust_bin_df <- subset(t(go_df),grepl("cov",colnames(go_df)))
toclust_bin_df <- t(apply(toclust_bin_df,1,as.numeric))
colnames(toclust_bin_df) <- 1:ncol(toclust_bin_df)
# Normalize min max for plotting scale
resC <- caCluster(toclust_bin_df, which="cols", dim=2, opt.part=TRUE)
resCfix <- caCluster(toclust_bin_df,part = 4, which="cols", dim=2, opt.part=TRUE)
names(resC) <- paste0("class_",1:length(resC))
names(resCfix) <- paste0("class_",1:length(resCfix))
col_ord <- ldply(resC,.id = "tad_class",function(X)data.frame(tad_idx_chr=names(X)))
col_ord_fix <- ldply(resCfix,.id = "tad_class_fix",function(X)data.frame(tad_idx_chr=names(X)))


mcols(Go3)["tad_classes_chipdensity_4classes_GEO"] <- NA
GenomicRanges::mcols(Go3[as.numeric(col_ord$tad_idx_chr)])["tad_classes_chipdensity_4classes_GEO"] <- as.character(col_ord_fix$tad_class)
go_l <- split(Go3,Go3$tad_classes_chipdensity_4classes_GEO)







### INHERIT TAD CLASS AT TILED 10kb ----
for (class in unique(Go3$tad_classes_chipdensity_4classes_GEO)) {
  
  tad <- subset(Go3, Go3$tad_classes_chipdensity_4classes_GEO==class)
  params.aname = paste0("tad_N2.old_25kb_GEO_",class,"_cov")
  params = list(
    bychr = TRUE
  )
  
  if (is.null(params$bychr)) {
    params$bychr <- T
  }
  
  # Go2 <- rangesDB2GO(filterRangesDB(rname = list(in_=list(value=params$rname, strict=TRUE))))
  mcols(Go)[[params.aname]] <- getCoverageDensity(Go, tad, bychr=params$bychr)
  
}



### PLOT ----
params.aname = c('tad_N2.old_25kb_GEO_class_1_cov', 'tad_N2.old_25kb_GEO_class_2_cov', 'tad_N2.old_25kb_GEO_class_3_cov', 'tad_N2.old_25kb_GEO_class_4_cov', 'eigen_pca1_N2.old')
params = list(
  signal = c('tad_N2.old_25kb_GEO_class_1_cov', 'tad_N2.old_25kb_GEO_class_2_cov', 'tad_N2.old_25kb_GEO_class_3_cov', 'tad_N2.old_25kb_GEO_class_4_cov'),
  ranker = 'eigen_pca1_N2.old',
  binsize = 50,
  binFun = 'mean',
  chr = c('I', 'II', 'III', 'IV', 'V', 'X'),
  chrAll = TRUE,
  legend = TRUE,
  height = 800,
  width = 2000
)

signal <- params[["signal"]]
ranker <- params[["ranker"]]
binsize <- params[["binsize"]]
binFun <- params[["binFun"]]
chr_v <- params[["chr"]]
chrAll <- params[["chrAll"]]
gg <- NULL

# Compute coverage density per peaks and per chromosomes
Go$chr <- as.character(GenomicRanges::seqnames(Go))
Go_dt <- data.table(data.frame(mcols(Go)))
cols <- colnames(Go_dt)[colnames(Go_dt) %in% signal]
Go_dt[,(cols) := lapply(.SD, function(x) ifelse(is.finite(x),x,NA)),.SDcols = cols]
Go_dt <- Go_dt[,.SD,.SDcols = colnames(Go_dt) %in% c(signal, ranker,"chr")]

# Digitalize by chromsoome and for all chromosomes
Go_dt[,digitalizedRankingByChrom:=dplyr::ntile(get(ranker),binsize),by = chr]
Go_dt[,digitalizedRankingAll:=dplyr::ntile(get(ranker),binsize)]
# Summarize peak density by digitalized groups from ranking value + Scaling
sum_All_dt <- data.table(Go_dt %>% dplyr::group_by(digitalizedRankingAll) %>% dplyr::summarise_all(binFun,na.rm=T))
sum_All_dt[,(signal) := lapply(.SD, function(x) as.vector(rangeMinMax(x, 0, 1))),.SDcols = colnames(sum_All_dt) %in% signal]

htmp_mat <-as.matrix(t(sum_All_dt %>% select(all_of(signal))))
anno_df <- data.frame(eigen = sum_All_dt %>% select(all_of(ranker)))
colnames(htmp_mat) <- 1:ncol(htmp_mat)
rownames(anno_df) <- colnames(htmp_mat)
gg[[paste0("HeatmapDensity.All_chromosomes.All")]] <- ggplotify::as.ggplot(pheatmap::pheatmap(htmp_mat,cluster_cols = F,color = hcl.colors(50,fixup = T, "RdBu",rev = T),border_color = NA,annotation_col = anno_df,cellwidth = 15,cellheight = 15))
