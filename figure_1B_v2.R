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
## SADDLEPLOT ----
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
# convert to Granges

cat("bigwig_signal_tile:... ")
# save.image("dev_l.RData")

Go <- bw_tile(bw = params$bw_p,
              bed = Go,
              genome = config.genome,
              name = params.aname)
### PLOT ----
params = list(
  anchor = 'eigen_pca1_N2.old',
  bait = 'eigen_pca1_N2.old',
  bait_binsize = 50,
  anchor_binsize = 50,
  h5path = list(N2 = 'data/HiC/N2-old_merged.bwa_mem.25kb.norm.KR.g2i.h5'),
  scaled = FALSE,
  agrMcol = 'median',
  agrSad = 'mean',
  agrSCM = 'median',
  title = '',
  legend = TRUE,
  fig_type = 'png',
  threads = 1,
  height = 800,
  width = 800
)


chr_v <- params$chrs # vector of chromosomes to be used in saddle, NULL if all chromosomes
log2_b <- ifelse(is.null(params$log),F,params$log)
bait_n <- params$bait_binsize # Size in bin of the bait (rows) saddle matrix. The greater the longer it will take to compute
anchor_n <- params$anchor_binsize # Size in bin of the a anchor(cols) saddle matrix. The greater the longer it will take to compute
genome <- params$genome # Reference genome
na_offset <- ifelse(is.null(params$na_offset), 20,params$na_offset)
diag_offset <- ifelse(is.null(params$diag_offset), 1,params$diag_offset)
####### HIC
h5path <- params$h5path # Path to H5 ContactMatrix HiC file
norm <- params$norm # Normalization of ContactMatrix. Used to dump contact matrix if H5 path is not provided
cond <- params$cond # Condition of ContactMatrix. Used to dump contact matrix if H5 path is not provided
rep <- params$rep # Replicate of ContactMatrix. Used to dump contact matrix if H5 path is not provided
####### AGGREG PARAMETERS
agr_mcols <- ifelse(!is.null(params$agrMcol),params$agrMcol, "median") # How should be summarized mcols ranker statistics ?
agr_sad <- ifelse(!is.null(params$agrSad),params$agrSad, "median") # How should be coomputed each saddle pixel ?
agr_scm <- ifelse(!is.null(params$agrSCM),params$agrSCM, "median") # how should be merged sub saddle matrices ( by chrom )

if(params$bait_binsize > params$anchor_binsize){
  bait <- params$anchor
  anchor <- params$bait
  bait_bs <- params$anchor_binsize
  anchor_bs <- params$bait_binsize
  ratio_ba <- bait_bs/anchor_bs
} else {
  bait <- params$bait
  anchor <- params$anchor
  bait_bs <- params$bait_binsize
  anchor_bs <- params$anchor_binsize
  ratio_ba <- bait_bs/anchor_bs
}

l_cm <- list()

# Filter Go with chr_v if not null, keep all chromp otherwise
if (!is.null(chr_v))
  Go <- subset(Go,seqnames%in%chr_v)
# Load h5
if (!is.null(h5path)) {
  for(h5.nm in names(h5path)){
    # get h5 name info by h5path
    cm_info <- rhdf5::h5ls(h5path[[h5.nm]])
    if (!is.null(chr_v))
      l_cm[[h5.nm]] <- sapply(cm_info$name[cm_info$name %in% chr_v], function(nm) { HDF5Array::HDF5Array(h5path[[h5.nm]], nm)})
    else # Keep all chrom
      l_cm[[h5.nm]] <- sapply(cm_info$name, function(nm) { HDF5Array::HDF5Array(h5path[[h5.nm]], nm)})
  }
} else { # unserialize
  # TODO : put loop here to load cm in l_cm
  mat_nm <- paste0(cond,"_",rep,"_",norm,"_",as.integer(res))
  filt_cm <- filterCmDB(cmname=list(in_=list(value=mat_nm,strict=T)))
  l_cm[["to_dev"]] <- dumpCM(filt_cm,name = chr_v) # Will return all chr if chr_v is NULL
}

sadmat_l <- list()
gg <- list()
# Create saddle matrix
for (nm in names(l_cm)){
  sadmat_l[[nm]] <- getSaddle2D(bin_go=Go,
                                bait=bait, # mcols name containing score for ranking genomic tiles
                                diag_offset=diag_offset, # mcols name containing score for ranking genomic tiles
                                na_offset=na_offset, # mcols name containing score for ranking genomic tiles
                                anchor=anchor, # mcols name containing score for ranking genomic tiles
                                bait_quant=bait_bs, # 30 digit (size of the saddle) by default
                                anchor_quant=anchor_bs, # 30 digit (size of the saddle) by default
                                cm_l=l_cm[[nm]], # list of H5 matrix (by chromosomes)
                                f_score=agr_mcols, # How should be summarized mcols ranker statistics ?
                                f_sad=agr_sad, # How should be coomputed each saddle pixel ?
                                f_cm=agr_scm, # How should we merge all sub-saddle matrix between chromosomes ?
                                log2=log2_b)
  
  
  # Do the plot
  m <- matrix2tibble(sadmat_l[[nm]][["sadmat_cm_mChr"]], bait_bs, anchor_bs)
  gg_center <- ggplot2::ggplot(m, ggplot2::aes(x=Anchor, y=Bait, fill=counts)) +
    ggplot2::geom_tile() + theme(axis.text.x = element_text(angle=90),axis.title.y = element_blank(),axis.title.x = element_blank())+
    scico::scale_fill_scico(palette = "lajolla",direction = -1)+ scale_x_discrete(position="top")
  gg_bait <- ggplot2::ggplot(data = sadmat_l[[nm]][["sum_stat_mChr_b"]], ggplot2::aes_string(x = "factor(quant_bait)", y = "tot_bait")) +
    ggplot2::geom_bar(stat = "identity", aes_string(fill = "tot_bait")) + ggplot2::theme_bw() +
    theme(legend.position="none", axis.ticks.y = element_blank(),axis.text.x = element_text(angle=90),
          axis.text.y = element_blank())+ scale_x_discrete(position="bottom",limits=rev) +
    ggplot2::labs(y = "Tile median values", x = bait) +  coord_flip() +scico::scale_fill_scico(palette = "bamako",direction = 1)
  
  gg_anchor <- ggplot2::ggplot(data = sadmat_l[[nm]][["sum_stat_mChr_a"]], ggplot2::aes_string(x = "factor(quant_anchor)", y = "tot_anchor")) +
    ggplot2::geom_bar(stat = "identity", aes_string(fill = "tot_anchor")) + ggplot2::theme_bw()+
    theme(legend.position="none", axis.ticks.x = element_blank(),axis.text.x = element_blank()) +
    ggplot2::labs(y = "Tile median values", x = anchor) + scale_x_discrete(position="top")+scico::scale_fill_scico(palette = "bamako",direction = 1)
  
  
  # Remove legend entry
  tmp_p <- cowplot::plot_grid(NULL,
                              gg_anchor ,
                              gg_bait ,
                              gg_center + theme(legend.position="none"),
                              rel_widths = c(1,4),
                              rel_heights = c(2,4),
                              ncol = 2, align = "hv")
  
  legend <- get_legend(gg_center + theme(legend.box.margin = margin(140,0, 20)))
  gg[[paste0("Saddle.",nm,".mergedChr")]] <- cowplot::plot_grid(tmp_p, legend, rel_widths = c(5,1))
}





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
  genome = 'ce11',
  desc = 'Bigwig signal of eigen_pca1 N2-old'
)
# convert to Granges

cat("bigwig_signal_tile:... ")
# save.image("dev_l.RData")

Go <- bw_tile(bw = params$bw_p,
              bed = Go,
              genome = config.genome,
              name = params.aname)
### MCOLS COVERAGE ----
peaks_lst <- c("H3K4me3_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "H3K4me1_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "H3K9me3_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "hpl2_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "H3K9me2_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "H3K27ac_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "H3K27me3_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "LEM2_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "LIN61_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed")

for (item in c(peaks_lst)) {
  
  params.aname = paste0(sub("(.*)_GEO.+","\\1",item),"_peaks_cov")
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
params = list(
  signal = c('H3K4me3_peaks_cov', 'H3K4me1_peaks_cov', 'H3K9me3_peaks_cov', 'hpl2_peaks_cov', 'H3K9me2_peaks_cov', 'H3K27ac_peaks_cov', 'H3K27me3_peaks_cov', 'LEM2_peaks_cov', 'LIN61_peaks_cov'),
  ranker = 'eigen_pca1_N2.old',
  binsize = 50,
  binFun = 'mean',
  chr = c('I', 'II', 'III', 'IV', 'V', 'X'),
  chrAll = TRUE,
  title = '',
  legend = TRUE
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
sum_byChr_dt <- data.table(Go_dt %>% dplyr::group_by(digitalizedRankingByChrom,chr) %>% dplyr::summarise_all(binFun,na.rm=T))
sum_byChr_dt[,(signal) := lapply(.SD, function(x) as.vector(rangeMinMax(x, 0, 1))), by = chr,.SDcols = colnames(sum_byChr_dt) %in% signal]
sum_meanChr_dt <- data.table(Go_dt %>% dplyr::group_by(digitalizedRankingByChrom) %>% dplyr::summarise_all(binFun,na.rm=T))
sum_meanChr_dt[,(signal) := lapply(.SD, function(x) as.vector(rangeMinMax(x, 0, 1))),.SDcols = colnames(sum_meanChr_dt) %in% signal]
sum_All_dt <- data.table(Go_dt %>% dplyr::group_by(digitalizedRankingAll) %>% dplyr::summarise_all(binFun,na.rm=T))
sum_All_dt[,(signal) := lapply(.SD, function(x) as.vector(rangeMinMax(x, 0, 1))),.SDcols = colnames(sum_All_dt) %in% signal]

htmp_mat <-as.matrix(t(sum_All_dt %>% select(all_of(signal))))
anno_df <- data.frame(eigen = sum_All_dt %>% select(all_of(ranker)))
colnames(htmp_mat) <- 1:ncol(htmp_mat)
rownames(anno_df) <- colnames(htmp_mat)
gg[[paste0("HeatmapDensity.All_chromosomes.All")]] <- ggplotify::as.ggplot(pheatmap::pheatmap(htmp_mat,cluster_cols = F,color = hcl.colors(50,fixup = T, "RdBu",rev = T),annotation_col = anno_df,cellwidth = 15,cellheight = 15))







## HEATMAP INTENSITY ----
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
### MCOLS INTENSITY ----
peaks_lst <- c("H3K4me3_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "H3K4me1_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "H3K9me3_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "hpl2_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "H3K9me2_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "H3K27ac_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "H3K27me3_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "LEM2_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "LIN61_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed")

for (item in c(peaks_lst)) {
  
  ip = sub("(.*)_GEO.+","\\1",item)
  params.aname = paste0(ip,"_N2_signal")
  
  params = list(
    bw_p = paste0("data/ChIPseq/",ip,"_GEO_merged.bwa_aln.rmdup.bamCompare.bw"),
    downstream = 0,
    upstream = 0,
    anchor = 'body'
  )
  Go2 <- loadranges(paste0("data/ChIPseq/",item),genome = "ce11")
  
  Go2 <- bw_to_mcols(bw = params$bw_p,
                     bed = Go2,
                     genome = config.genome,
                     anchor = params$anchor,
                     name = params.aname,
                     upstream = params$upstream,
                     downstream = params$downstream)
  params.aname <- paste0(ip,"_GEO_N2_max_signal")
  params = list(
    aname = paste0(ip,'_N2_signal'),
    fun = 'max',
    fun_opts = 'na.rm=T'
  )
  
  Go <- aggregate_ranges(Go, Go2, out_mc=params.aname,  subject_mc = params$aname, fun_agr = params$fun,fol_opts = params$fol_opts, fun_opts=params$fun_opts)
}



### PLOT ----
params = list(
  signal = c('H3K4me3_GEO_N2_max_signal', 'H3K4me1_GEO_N2_max_signal', 'H3K9me3_GEO_N2_max_signal', 'hpl2_GEO_N2_max_signal', 'H3K9me2_GEO_N2_max_signal', 'H3K27ac_GEO_N2_max_signal', 'H3K27me3_GEO_N2_max_signal', 'LEM2_GEO_N2_max_signal', 'LIN61_GEO_N2_max_signal'),
  ranker = 'eigen_pca1_N2.old',
  binsize = 50,
  binFun = 'mean',
  chr = c('I', 'II', 'III', 'IV', 'V', 'X'),
  chrAll = TRUE,
  title = '',
  legend = TRUE
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
sum_byChr_dt <- data.table(Go_dt %>% dplyr::group_by(digitalizedRankingByChrom,chr) %>% dplyr::summarise_all(binFun,na.rm=T))
sum_byChr_dt[,(signal) := lapply(.SD, function(x) as.vector(rangeMinMax(x, 0, 1))), by = chr,.SDcols = colnames(sum_byChr_dt) %in% signal]
sum_meanChr_dt <- data.table(Go_dt %>% dplyr::group_by(digitalizedRankingByChrom) %>% dplyr::summarise_all(binFun,na.rm=T))
sum_meanChr_dt[,(signal) := lapply(.SD, function(x) as.vector(rangeMinMax(x, 0, 1))),.SDcols = colnames(sum_meanChr_dt) %in% signal]
sum_All_dt <- data.table(Go_dt %>% dplyr::group_by(digitalizedRankingAll) %>% dplyr::summarise_all(binFun,na.rm=T))
sum_All_dt[,(signal) := lapply(.SD, function(x) as.vector(rangeMinMax(x, 0, 1))),.SDcols = colnames(sum_All_dt) %in% signal]

htmp_mat <-as.matrix(t(sum_All_dt %>% select(all_of(signal))))
anno_df <- data.frame(eigen = sum_All_dt %>% select(all_of(ranker)))
colnames(htmp_mat) <- 1:ncol(htmp_mat)
rownames(anno_df) <- colnames(htmp_mat)
gg[[paste0("HeatmapDensity.All_chromosomes.All")]] <- ggplotify::as.ggplot(pheatmap::pheatmap(htmp_mat,cluster_cols = F,color = hcl.colors(50,fixup = T, "RdBu",rev = T),annotation_col = anno_df,cellwidth = 15,cellheight = 15))

