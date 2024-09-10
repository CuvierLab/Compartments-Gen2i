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


## HEATMAP DENSITY ----
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
  genome = 'ce11'
)

Go <- bw_tile(bw = params$bw_p,
              bed = Go,
              genome = config.genome,
              name = params.aname)


for (condition in paste0(c("hpl2","lin61","hpl2-lin61","met2-set25-set32","I158A"),"-old")) {
  
  dot_cond = sub("-","\\.",condition)
  params.aname = paste0("eigen_pca1_",dot_cond)
  params = list(
    bw_p = paste0('data/HiC/',condition,'_merged.bwa_mem.25kb.pca1.bw'),
    genome = 'ce11'
  )
  
  Go <- bw_tile(bw = params$bw_p,
                bed = Go,
                genome = config.genome,
                name = params.aname)
  
  #### DELTA EIGEN ----
  params.aname = paste0("delta_eigen_pca1_",dot_cond,"_N2.old")
  params.process = 
    params = list(
      operation = 'subtract'
    )
  params.dependencies = c('eigen_pca1_' %+% dot_cond, 'eigen_pca1_N2.old')
  mcols(Go)[params.aname] <- NA
  method <- ifelse(!is.null(params$method),params$method, params$operation)
  mcols(Go)[params.aname] <- Gooperation(unlist(mcols(Go)[[params.dependencies[1]]]),
                                         unlist(mcols(Go)[[params.dependencies[2]]]),
                                         method = method)
}


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
  mcols(Go)[[params.aname]] <- getCoverageDensity(Go, Go2, bychr=params$bychr)
}


### PLOT ----
params.aname = c('H3K4me3_peaks_cov', 'H3K4me1_peaks_cov', 'H3K9me3_peaks_cov', 'hpl2_peaks_cov', 'H3K9me2_peaks_cov', 'H3K27ac_peaks_cov', 'H3K27me3_peaks_cov', 'LEM2_peaks_cov', 'lin61_peaks_cov', 'delta_eigen_pca1_I158A.old_N2.old')
for (condition in paste0(c("hpl2","lin61","hpl2.lin61","met2.set25.set32","I158A"),".old")) {
  params = list(
    signal = c('H3K4me3_peaks_cov', 'H3K4me1_peaks_cov', 'H3K9me3_peaks_cov', 'hpl2_peaks_cov', 'H3K9me2_peaks_cov', 'H3K27ac_peaks_cov', 'H3K27me3_peaks_cov', 'LEM2_peaks_cov', 'lin61_peaks_cov'),
    ranker = paste0('delta_eigen_pca1_',condition,'_N2.old'),
    binsize = 30,
    binFun = 'mean',
    chr = c('I', 'II', 'III', 'IV', 'V', 'X'),
    chrAll = TRUE,
    legend = TRUE,
    height = 500,
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
  sum_meanChr_dt <- data.table(Go_dt %>% dplyr::group_by(digitalizedRankingByChrom) %>% dplyr::summarise_all(binFun,na.rm=T))
  sum_meanChr_dt[,(signal) := lapply(.SD, function(x) as.vector(rangeMinMax(x, 0, 1))),.SDcols = colnames(sum_meanChr_dt) %in% signal]
  
  htmp_mat <-as.matrix(t(sum_meanChr_dt %>% select(all_of(signal))))
  anno_df <- data.frame(eigen = sum_meanChr_dt %>% select(all_of(ranker)))
  colnames(htmp_mat) <- 1:ncol(htmp_mat)
  rownames(anno_df) <- colnames(htmp_mat)
  gg[[paste0("HeatmapDensity.Mean_chromosomes.All.",condition)]] <- ggplotify::as.ggplot(pheatmap::pheatmap(htmp_mat,cluster_cols = F,color = hcl.colors(50,fixup = T, "RdBu",rev = T),annotation_col = anno_df,cellwidth = 15,cellheight = 15))
  
}
