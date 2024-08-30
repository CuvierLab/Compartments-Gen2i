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

## BEDFILE----
Go <- loadranges("data/ChIPseq/H3K9me23_lin61_hpl2_optimal_N2_peakset.bed",genome = "ce11")

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
    Go <- bw_to_mcols(bw = params$bw_p,
                      bed = Go,
                      genome = config.genome,
                      anchor = params$anchor,
                      name = params.aname,
                      upstream = params$upstream,
                      downstream = params$downstream)
  }
  
}


## TILED GENOME ----
Go2 <- loadranges("data/ChIPseq/ce11_tiled_25kb.bed",genome = "ce11")
### LOAD EIGEN ----
params.aname = "eigen_pca1_N2"
params = list(
  bw_p = 'data/HiC/N2_merged.bwa_mem.25kb.pca1.bw',
  genome = 'ce11'
)
Go2 <- bw_tile(bw = params$bw_p,
               bed = Go2,
               genome = config.genome,
               name = params.aname)
### EIGEN CLASSES ----
params.aname = "eigen_pca1_N2_50tile"
params = list(
  ngroups = 50
)
params.dependencies = "eigen_pca1_N2"
params$ngroups <- as.integer(params$ngroups)

mcols(Go2)[params.aname] <- NA
mcols(Go2)[params.aname] <- dplyr::ntile(unlist(mcols(Go2)[[params.dependencies]]),params$ngroups)

params.aname = "eigen_pca1_N2_grps_shifted_A_B"
params = list(
  values = c(12.0, 16.0),
  names = c('B', 'na', 'A')
)
params.dependencies = "eigen_pca1_N2_50tile"
# convert NA to 0
f <- is.na(mcols(Go2)[[params.dependencies]])
mcols(Go2)[[params.dependencies]][f] <- 0
mcols(Go2)[params.aname] <- NA

for(i in 1:length(params$values)){
  if(i == 1){
    f <- as.vector(mcols(Go2)[[params.dependencies]] <= params$values[i])
    if(any(f))
      mcols(Go2[f])[params.aname] <- params$names[i]
  }
  else{
    f <- as.vector(mcols(Go2)[[params.dependencies]]>params$values[i-1] & mcols(Go2)[[params.dependencies]]<=params$values[i])
    if(any(f))
      mcols(Go2[f])[params.aname] <- params$names[i]
  }
}
f <- as.vector(mcols(Go2)[[params.dependencies]]>params$values[i])
if(any(f))
  mcols(Go2[f])[params.aname] <- params$names[i+1]



# INHERITS EIGEN CLASSES ----
params.aname = "eigen_pca1_N2_grps_shifted_A_B"
params = list(
  aname = 'eigen_pca1_N2_grps_shifted_A_B'
)
idx <- GenomicRanges::nearest(Go, Go2)
# Initlalize query Mcol
mcols(Go)[params.aname] <- NA
# if there is NA somewhere (mostly cause when using seqnames in one of the ranges and not the other)
# you should remove those indexes from both Go ranges
Go_na_idx <- which(is.na(idx))
if (length(Go_na_idx)) {
  idx <- idx[!is.na(idx)]
  mcols(Go[-Go_na_idx])[params.aname] <- unlist(mcols(Go2[idx])[params$aname])[[1]]
} else {
  mcols(Go)[params.aname] <- unlist(mcols(Go2[idx])[params$aname])
}

# If a distance threshold is set to retrieve nearest:
tokeep <- rep(T, length(Go))
if (!is.null(params$distThreshold)) tokeep <- distanceToNearest(Go, Go2)@elementMetadata$distance <= params$distThreshold
mcols(Go)[params.aname][!tokeep,] <- NA


# PLOT ----
params.aname = c('hpl2_I158A_newzs_signal', 'hpl2_hpl2.lin61_newzs_signal', 'hpl2_lin61_newzs_signal', 'hpl2_met2.set25.set32_newzs_signal', 'lin61_I158A_newzs_signal', 'lin61_hpl2_newzs_signal', 'lin61_lin61_newzs_signal', 'lin61_met2.set25.set32_newzs_signal', 'H3K9me2_I158A_newzs_signal', 'H3K9me2_hpl2_newzs_signal', 'H3K9me2_hpl2.lin61_newzs_signal', 'H3K9me2_lin61_newzs_signal', 'H3K9me2_met2.set25.set32_newzs_signal', 'H3K9me3_I158A_newzs_signal', 'H3K9me3_hpl2_newzs_signal', 'H3K9me3_hpl2.lin61_newzs_signal', 'H3K9me3_lin61_newzs_signal', 'H3K9me3_met2.set25.set32_newzs_signal', 'eigen_pca1_N2_grps_shifted_A_B')
params = list(
  comp = 'eigen_pca1_N2_grps_shifted_A_B',
  legend = TRUE,
  height = 800,
  width = 800
)

gg <- NULL
Go_dt <- go2dt(Go)

# Loop on eigen compartments
comp = "B" 
Go_s_dt <- subset(Go_dt, get(params$comp) == comp)
# keep only zscore columns
Go_s_dt <- Go_s_dt[,.SD,.SDcols = colnames(Go_s_dt) %like% "newzs"]
gg[[paste0("Network_analyzes.Label.Compartment_Filter_",comp)]] <- ggplotify::as.ggplot(Go_s_dt %>% correlate()  %>% network_plot(curved = F,repel=T))

for (i in seq_along(names(Go_s_dt))) {
  setnames(Go_s_dt, old = names(Go_s_dt)[i], new = paste(rep(" ", i), collapse = ""))
}
gg[[paste0("Network_analyzes.NoLabel.Compartment_Filter_",comp)]] <- ggplotify::as.ggplot(Go_s_dt %>% correlate()  %>% network_plot(curved = F,repel=T))

# =========================================================== =
