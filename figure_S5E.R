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
## GENOMIC DISTANCE ----
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

## BEDFILE----
Go2 <- loadranges("data/ChIPseq/hpl2_lin61_N2_narrow_union_peakset.bed",genome = "ce11")

### LOAD ZSCORE BIGWIG----
for (ip in c("H3K9me3","H3K9me2","lin61","hpl2")) {
  for (condition in c("hpl2.lin61","I158A","met2.set25.set32","hpl2","lin61")) {
    if (condition == "hpl2" & ip == "hpl2") next()
    if (condition == "hpl2.lin61" & ip == "lin61") next()
    
    params.aname = paste0(ip,"_",condition,"_newzs_signal")
    c_nodot <- gsub("\\.","-",condition)
    params = list(
      bw_p = paste0("data/ChIPseq/",ip,"_",c_nodot,"_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw"),
      downstream = 0,
      upstream = 0,
      anchor = 'body'
    )
    cat("bigwig_signal:... ")
    Go2 <- bw_to_mcols(bw = params$bw_p,
                       bed = Go2,
                       genome = config.genome,
                       anchor = params$anchor,
                       name = params.aname,
                       upstream = params$upstream,
                       downstream = params$downstream)
  }
  
}


### INHERITS BIGWIG SIGNAL FROM DOUBLE PEAKS ----

for (ip in c("H3K9me3","H3K9me2","lin61","hpl2")) {
  for (condition in c("hpl2.lin61","I158A","met2.set25.set32","hpl2","lin61")) {
    if (condition == "hpl2" & ip == "hpl2") next()
    if (condition == "hpl2.lin61" & ip == "lin61") next()
    
    params.aname = "hpl2_lin61_u_nr_WF_N2_" %+% paste0(ip,"_",condition,"_newzs_signal")
    
    params = list(
      aname = paste0(ip,"_",condition,"_newzs_signal"),
      fun = 'sum',
      fun_opts = 'na.rm=T'
    )
    Go <- aggregate_ranges(Go, Go2, out_mc=params.aname,  subject_mc = params$aname, fun_agr = params$fun,fol_opts = params$fol_opts, fun_opts=params$fun_opts)
  }
}




### PLOT ----

go <- go2dt(subset(Go, Go$eigen_pca1_N2.old_50tile <= 12 ))
gg <- list()
for (ip in c("H3K9me3","H3K9me2","lin61")) {
  for (condition in c("hpl2.lin61","I158A","met2.set25.set32","hpl2","lin61")) {
    if (condition == "hpl2") next()
    if (condition == "hpl2.lin61" & ip == "lin61") next()
    
    params.aname = "hpl2_lin61_u_nr_WF_N2_" %+% paste0(ip,"_",condition,"_newzs_signal")
    
    
    params.aname = c('hpl2_lin61_u_nr_WF_N2_hpl2_' %+% condition %+% '_newzs_signal', "hpl2_lin61_u_nr_WF_N2_" %+% paste0(ip,"_",condition,"_newzs_signal"))
    params = list(
      x = params.aname[1],
      y = params.aname[2],
      lab_x = 'Aggregated ChIP-seq Zscore'%+% paste0("hpl2 in cond ",condition),
      lab_y = 'Aggregated ChIP-seq Zscore'%+% paste0(ip," in cond ",condition)
    )
    
      sgg <- ggplot(go,
                    aes_string(x=params$x,
                               y=params$y,
                               col=params$col))
    
    
    sgg <- sgg + geom_point()
    gg[[ip%+%condition]] <- sgg + stat_regline_equation() + stat_smooth(method = "lm") + stat_cor(label.x.npc = "center")
    
  }
}
