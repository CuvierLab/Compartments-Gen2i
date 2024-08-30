# WORKING DIRECTORY ----
setwd(dir = "config/src/R/github/")

# SOURCE ----
source("manuscript_lib.R")



## SIMPLE PEAK ANALYZES ----
### TILE 1KB ----
Go <- loadranges("data/ChIPseq/ce11_tiled_1kb.bed","ce11")
### GEO PEAKS ----

peaks_lst <- c("H3K9me3_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "hpl2_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "H3K9me2_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "LIN61_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed")

for (pk in peaks_lst) {
  
  
  ip  <-  sub("(.*)_GEO.+","\\1",pk)
  
  params.aname = ip %+% "_peaks"
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
params.aname = c('hpl2_peaks', 'LIN61_peaks', 'H3K9me2_peaks', 'H3K9me3_peaks', 'hpl2_WF_narrowPeak', 'lin61_WF_narrowPeak', 'H3K9me2_WF_narrowPeak', 'H3K9me3_WF_narrowPeak')
params = list(
  row = c('hpl2_peaks', 'LIN61_peaks', 'H3K9me2_peaks', 'H3K9me3_peaks'),
  col = c('hpl2_WF_narrowPeak', 'lin61_WF_narrowPeak', 'H3K9me2_WF_narrowPeak', 'H3K9me3_WF_narrowPeak'),
  height = 800,
  width = 800
)

Go <- go2dt(Go)
skip_go <- T
df <- fisher.bool.bool(df = Go, grad_row = params$row, grad_col = params$col,identity = F)
df$ROW <- "Geo_peaks"
setDT(df)
df[, pv := as.numeric(pv)]
mn.pv <- df[pv>0, min(pv)]
df[, pv := pv+mn.pv]
ok <- fisher_plot_asym_pdf(DF = df, main = params$title, norm = F, range = c(-2,2))



## DOUBLE PEAK ANALYZES ----
## TILE 1KB ----
Go <- loadranges("data/ChIPseq/ce11_tiled_1kb.bed","ce11")
### GEO PEAKS ----

peaks_lst <- c("H3K9me3_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "hpl2_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "H3K9me2_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
               "LIN61_GEO.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed")

for (pk in peaks_lst) {
  
  
  ip  <-  sub("(.*)_GEO.+","\\1",pk)
  
  params.aname = ip %+% "_peaks"
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

#### DOUBLE PEAKS ----
for (ip1 in c("H3K9me3","H3K9me2")) {
  for (ip2 in c("hpl2","LIN61")) {
    params.aname <- ip1%+%"_"%+%ip2%+%"_peaks"
    
    params.dependencies = c(ip1%+%"_peaks", ip2%+%"_peaks")
    
    mcols(Go)[params.aname] <- NA
    method <- "and"
    mcols(Go)[params.aname] <- Gooperation(unlist(mcols(Go)[[params.dependencies[1]]]),
                                           unlist(mcols(Go)[[params.dependencies[2]]]),
                                           method = method)
  }
}




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

#### DOUBLE PEAKS ----
for (ip1 in c("H3K9me3","H3K9me2")) {
  for (ip2 in c("hpl2","lin61")) {
    params.aname <- ip1%+%"_"%+%ip2%+%"_WF_narrowPeak"
    
    params.dependencies = c(ip1%+%"_WF_narrowPeak", ip2%+%"_WF_narrowPeak")
    
    mcols(Go)[params.aname] <- NA
    method <- "and"
    mcols(Go)[params.aname] <- Gooperation(unlist(mcols(Go)[[params.dependencies[1]]]),
                                           unlist(mcols(Go)[[params.dependencies[2]]]),
                                           method = method)
  }
}





### PLOT ----
params.aname = c('H3K9me2_hpl2_peaks', 'H3K9me2_LIN61_peaks', 'H3K9me3_hpl2_peaks', 'H3K9me3_LIN61_peaks', 'H3K9me2_hpl2_WF_narrowPeak', 'H3K9me2_lin61_WF_narrowPeak', 'H3K9me3_hpl2_WF_narrowPeak', 'H3K9me3_lin61_WF_narrowPeak')
params = list(
  row = c('H3K9me2_hpl2_peaks', 'H3K9me2_LIN61_peaks', 'H3K9me3_hpl2_peaks', 'H3K9me3_LIN61_peaks'),
  col = c('H3K9me2_hpl2_WF_narrowPeak', 'H3K9me2_lin61_WF_narrowPeak', 'H3K9me3_hpl2_WF_narrowPeak', 'H3K9me3_lin61_WF_narrowPeak'),
  height = 800,
  width = 800
)

Go <- go2dt(Go)
skip_go <- T
df <- fisher.bool.bool(df = Go, grad_row = params$row, grad_col = params$col,identity = F)
df$ROW <- "Geo_peaks"
setDT(df)
df[, pv := as.numeric(pv)]
mn.pv <- df[pv>0, min(pv)]
df[, pv := pv+mn.pv]
ok <- fisher_plot_asym_pdf(DF = df, main = params$title, norm = F, range = c(-2,2))

