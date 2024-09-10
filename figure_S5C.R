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



## SIMPLE PEAK ANALYZES ----
### TILE 1KB ----
Go <- loadranges("data/ChIPseq/ce11_tiled_1kb.bed","ce11")
### WF PEAKS ----

peaks_lst <- c("H3K9me3",
               "hpl2",
               "H3K9me2",
               "lin61")

for (ip in peaks_lst) {
  
  
  pk  <-  ip%+%"_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed"
  
  params.aname = ip %+% "_WF_narrowPeak"
  params = list(
    genome = 'ce11',
    downstream = 0,
    upstream = 0,
    anchor = 'body'
  )
  Go2 <- g2i:::loadranges("data/ChIPseq/"%+%pk, genome = params$genome)
  
  params$ignore.strand <- ifelse(is.null(params$ignore.strand),T,params$ignore.strand)
  
  if(is.null(params$anchor_from))
    params$anchor_from <- "body"
  if(is.null(params$upstream_from))
    params$upstream_from <- 0
  if(is.null(params$downstream_from))
    params$downstream_from <- 0
  if(is.null(params$anchor_to))
    params$anchor_to <- "body"
  if(is.null(params$upstream_to))
    params$upstream_to <- 0
  if(is.null(params$downstream_to))
    params$downstream_to <- 0
  
  # For easier code handling, In case where downstream and upstream opts are 0 apply dummy reframe
  Go <- reframeGR(Go, params$anchor_from, params$upstream_from, params$downstream_from)
  Go2 <- reframeGR(Go2, params$anchor_to, params$upstream_to, params$downstream_to)
  
  if(params$ignore.strand) f <- GenomicRanges::findOverlaps(Go, Go2,ignore.strand=T)@from
  else  f <- GenomicRanges::findOverlaps(Go, Go2)@from
  
  mcols(Go)[params.aname] <- F
  mcols(Go[f])[params.aname] <- T
}




### PLOT ----
params.aname = c('hpl2_WF_narrowPeak', 'lin61_WF_narrowPeak', 'H3K9me2_WF_narrowPeak', 'H3K9me3_WF_narrowPeak')
Go <- go2dt(Go)
l <- list()
for(i in params.aname){
  setkeyv(Go, i)
  l[[i]] <- Go[Go[[i]] == T, name]
}

ggVennDiagram(l,set_color = RColorBrewer::brewer.pal(4,"Set1"))


