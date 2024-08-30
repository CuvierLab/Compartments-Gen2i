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

## Deeptools ChIP-seq Zscore Matrix profiles ----
### NOT RUN ###
'computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 \
--downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 \
-R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed \
-S data/ChIPseq/${ip}_${condition}_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw \
-o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/{ip}_${condition}_newzs_prof.matrix.gz \
--outFileNameMatrix data/ChIPseq/{ip}_${condition}_newzs_prof.txt) > matrix.log 2>&1'

## BEDFILE----
Go <- loadranges("data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed",genome = "ce11")

## ChIP-seq Zscore Matrix profile to mcols ----
for (condition in c("hpl2.lin61","I158A","met2.set25.set32","hpl2","lin61")) {
  for (ip in c("H3K9me2","lin61","hpl2","LEM2")) {
    if (condition == "met2.set25.set32" & ip == "LEM2") next()
    if (condition == "hpl2" & ip == "hpl2") next()
    if (condition == "hpl2.lin61" & ip == "lin61") next()
    params.aname = paste0(ip,"_",condition,"_newzs_prof")
    annotxt_dt <-  data.table::fread(paste0("data/ChIPseq/",ip,"_",condition,"_newzs_prof.txt"),check.names = F)[, name:=as.character(name)]
    data.table::setkey(annotxt_dt, name) # Set keys to retrieve ordering by name column
    ao_dt <- annotxt_dt[as.character(Go$name), !"name"] # Sort to match GO ordering
    mcols(Go)[[params.aname]] <- ao_dt
  }
}

# ChIP-seq zscore signal to mcols ----
for (condition in c("hpl2.lin61","I158A","met2.set25.set32","hpl2","lin61")) {
  ip <- "H3K9me2"
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


## Inherit eigen values ----
### MCOLS Eigen values ----
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

## Plot ----
Go <- subset(Go, eigen_pca1_N2.old<=0)
gg <- list()
for (condition in c("hpl2.lin61","I158A","met2.set25.set32","hpl2","lin61")) {
  for (ip in c("H3K9me2","lin61","hpl2","LEM2")) {
    if (condition == "met2.set25.set32" & ip == "LEM2") next()
    if (condition == "hpl2" & ip == "hpl2") next()
    if (condition == "hpl2.lin61" & ip == "lin61") next()
  
    
    cat(paste0("computeMatrix scale-regions --sortRegions keep --binSize 10 --upstream 2000 --downstream 2000 --regionBodyLength 200 --unscaled5prime 0 --unscaled3prime 0 --missingDataAsZero --numberOfProcessors 1 -R data/ChIPseq/H3K9me2_H3K9me3_N2_narrow_union_peakset.bed -S data/ChIPseq/",ip,"_",condition,"_vs_N2_zscore.bwa_aln.rmdup.bamCompare.bw -o data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/",ip,"_",condition,"_newzs_prof.matrix.gz --outFileNameMatrix data/ChIPseq/H3K9me2_H3K9me3_u_nr_WF_N2/",ip,"_",condition,"_newzs_prof.matrix.tab 
"))
      
  }
}
    height = 800
    width = 300
    params.aname = c(paste0(ip,'_',condition,'_newzs_prof'), paste0('H3K9me2_',condition,'_newzs_signal'))
    params = list(
      order = params.aname[2],
      decreasing = FALSE,
      yscale = c(-2.0, 2.0),
      log = FALSE,
      colorscale = c(-2.0, 2.0),
      filter = list(lst = list(f1 = list(aname = 'eigen_pca1_N2.old', operator = '<=', value = 0))),
      legend = TRUE
      )
    
    m <- mcols(Go)[[params.aname[[1]]]]
    if(height > nrow(m))
      height = nrow(m)
    
    if(width > ncol((m)))
      width = ncol((m))
    
    if(!is.null(params$order)){
      f <- order(mcols(Go)[[params$order]])
      rd=redim_matrix(as.matrix(m[f,]), target_height=height, target_width=width, n_core=1)
    }
    else
      rd=redim_matrix(as.matrix(m[order(rowSums(m))]), target_height=height, target_width=width, n_core=1)
    
    ratio_p <- 160/min(width, height)
    if (ratio_p > 1){
      height <- height * ratio_p
      width <- width * ratio_p
    }
    
    g <- reshape2::melt(rd)
    min.val <- min(g$value, na.rm = T)
    setDT(g)
    g[is.na(value), value := min.val]
    
    if (!is.null(params$yscale)){
      min <- min(params$yscale)
      max <- max(params$yscale)
      g[value<min, value := min]
      g[value>max, value := max]
    }
    
    if(params$log)
      sgg <- ggplot(g, aes(x=Var2, y=Var1, fill=log(value)))
    else
      sgg <- ggplot(g, aes(x=Var2, y=Var1, fill=value))
    
    sgg <- sgg + geom_raster( ) + theme_minimal() +
      theme(
        plot.margin=margin(grid::unit(0, "cm")),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = element_blank(),
        plot.caption = element_text(hjust=0, size=8, face = "italic"),
        plot.subtitle = element_text(hjust=0, size=8),
        plot.title   = element_text(hjust=0, size=12, face="bold"))
    
    if (!is.null(params$colorscale))
      sgg <- sgg + scale_fill_viridis_c(guide="colorbar", limits = c(min(params$colorscale),max(params$colorscale)))
    else
      sgg <- sgg + scale_fill_viridis_c( guide="colorbar")
    
    gg[[paste0(ip,"_",condition,"_heatmap")]] <- sgg
  }
}
# ----
