# HEADER ====================================================== =
# GEN2I
# Tue Jun 18 18:22:27 2024
# R version : R version 4.4.0 (2024-04-24)
# System : x86_64, linux-gnu
# ============================================================= =


# # WORKING DIRECTORY ----
setwd(dir = "config/src/R/github/")

# SOURCE ----
source("manuscript_lib.R")

# RUN ----

## Get nearest TAD IDs at 5kb with TAD of N2 from WF
Go <- loadranges("data/ChIPseq/ce11_tiled_5kb.bed",genome = "ce11")
params.aname = "chr_bin"
params = list(bychr= T)
mcol <- NULL
if(params$bychr) {
  mcol <- unlist(lapply(split(Go,seqnames(Go)),function(chr){1:length(chr)}))
} else mcol <- 1:length(Go)
mcols(Go)[params.aname] <- NA
mcols(Go)[params.aname] <- mcol


Go2 <- loadranges("data/HiC/N2_merged.bwa_mem.5kb_domains.bed",genome = "ce11")
params.aname = "bin"
params = list(bychr= F)
mcol <- NULL
if(params$bychr) {
  mcol <- unlist(lapply(split(Go2,seqnames(Go2)),function(chr){1:length(chr)}))
} else mcol <- 1:length(Go2)
mcols(Go2)[params.aname] <- NA
mcols(Go2)[params.aname] <- mcol



params = list(aname = "bin")
params.aname = "tad_id_N2_5kb"

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


## CONTACT MATRIX PLOT ----
params.rname <- "tads"
params = list(
  upperTAD_id = 'tad_id_N2_5kb',
  lowerTAD_id = 'tad_id_N2_5kb',
  upperMat = list(N2 = 'data/HiC/N2_merged.bwa_mem.5kb.norm.KR.g2i.h5'),
  lowerMat = list(N2 = 'data/HiC/N2_merged.bwa_mem.5kb.norm.KR.g2i.h5'),
  na_offset = 20,
  diag_offset = 0,
  palette = 'turbo',
  quant = 70,
  transform = 'log2',
  title = '',
  legend = TRUE,
  fig_type = 'png',
  threads = 1,
  height = 800,
  width = 800
)
upperMat <- params[["upperMat"]]
lowerMat <- params[["lowerMat"]]
vMin <- ifelse(is.null(params[["vMin"]]), -Inf, params[["vMin"]])
vMax <- ifelse(is.null(params[["vMax"]]), Inf, params[["vMax"]])
na_offset <- ifelse(is.null(params[["na_offset"]]), 0, params[["na_offset"]])
diag_offset <- ifelse(is.null(params[["diag_offset"]]), 1, params[["diag_offset"]])
## GRAPHICAL PARAMETERS
palette <- params[["palette"]]
quant <- params[["quant"]]
transform <- params[["transform"]]
redim_mat <- params[["redim_mat"]]

# Create dt for better handling of each comparisons
a <- tibble::enframe(params[["upperMat"]]) %>% tidyr::unnest(cols = value) %>% data.table::data.table()
b <- tibble::enframe(params[["lowerMat"]]) %>% tidyr::unnest(cols = value) %>% data.table::data.table()
dt <- cbind(a,b)
data.table::setnames(dt, c("nameUp", "pathUp", "nameLow","pathLow"))

Go <- go2dt(Go)

gg <- NULL
for (i in 1:nrow(dt)) {
  nm_u <- as.character(dt[i, "nameUp"])
  nm_l <- as.character(dt[i, "nameLow"])
  
  h5u_l <- list(as.vector(unlist(dt[i, "pathUp"])))
  h5l_l <- list(as.vector(unlist(dt[i, "pathLow"])))
  
  names(h5u_l) <- nm_u
  names(h5l_l) <- nm_l
  
  cm_u_lbc <- load_h5_matrix(h5path = h5u_l) # lbc = list by chromosomes
  cm_l_lbc <- load_h5_matrix(h5path = h5l_l)
  
  chr_v <- intersect(names(cm_u_lbc[[nm_u]]), names(cm_l_lbc[[nm_l]]))
  
  chr_v <- chr_v[!grepl("mt|Mt|MT|chrM|M",chr_v)]
  # load("pl")
  for(chr in chr_v){
    
    cm_u <- as.matrix(as.matrix(cm_u_lbc[[nm_u]][[chr]]))
    cm_l <- as.matrix(as.matrix(cm_l_lbc[[nm_l]][[chr]]))
    
    max_bin <- dim(cm_u)[1] - na_offset # Remove NA's extend use in APA
    cm_u <- cm_u[(na_offset+1):max_bin, (na_offset+1):max_bin] # Remove NA's extend use in APA
    cm_l <- cm_l[(na_offset+1):max_bin, (na_offset+1):max_bin] # Remove NA's extend use in APA
    
    # Get TAD info to draw on ggplot object
    if (params.rname != ""){
      u_tad_dt <- Go[!is.na(get(params$upperTAD_id)) & seqnames==chr, .(chr_bin, u_tad_id=get(params$upperTAD_id))]
      u_tad_left_dt <- copy(u_tad_dt)
      u_tad_left_dt[,x:=chr_bin[1],by = u_tad_id]
      u_tad_left_dt[,y:=chr_bin]
      u_tad_right_dt <- copy(u_tad_dt)
      u_tad_right_dt[,y:=chr_bin[.N],by = u_tad_id]
      u_tad_right_dt[,x:=chr_bin]
      
      l_tad_dt <- Go[!is.na(get(params$lowerTAD_id)) & seqnames==chr, .(chr_bin, l_tad_id=get(params$lowerTAD_id))]
      l_tad_left_dt <- copy(l_tad_dt)
      l_tad_left_dt[,y:=chr_bin[1],by = l_tad_id]
      l_tad_left_dt[,x:=chr_bin]
      l_tad_right_dt <- copy(l_tad_dt)
      l_tad_right_dt[,x:=chr_bin[.N],by = l_tad_id]
      l_tad_right_dt[,y:=chr_bin]
    } 
    
    cm_u[cm_u == 0] <- NA
    cm_l[cm_l == 0] <- NA
    # 
    getKDiag(cm_u, -diag_offset:diag_offset) <- NA
    getKDiag(cm_l, -diag_offset:diag_offset) <- NA
    # remove lower and upper accordingly
    
    cm_u[lower.tri(cm_u)] <- 0
    cm_l[upper.tri(cm_l)] <- 0
    getKDiag(cm_l, 0) <- 0
    cm <- if(!is.null(redim_mat)) redim_matrix(cm_u+cm_l, target_height=redim_mat, target_width=redim_mat, n_core=12) else  cm_u+cm_l
    cm <- if(!is.null(transform)) get(transform)(cm) else cm
    
    # Melt data to go ggplot
    cm <- data.table::data.table(melt(cm))
    colnames(cm) <- c("x","y","value")
    
    cm[,value:=ifelse(value>vMin,value,vMin)] 
    cm[,value:=ifelse(value<vMax,value,vMax)] 
    
    if (!is.null(quant)) {
      cmq_u <- matrix(dplyr::ntile(as.vector(cm_u), quant), nrow = nrow(cm_u), ncol = ncol(cm_u))
      cmq_l <- matrix(dplyr::ntile(as.vector(cm_l), quant), nrow = nrow(cm_l), ncol = ncol(cm_l))
      # remove lower and upper accordingly
      cmq_u[lower.tri(cmq_u)] <- 0
      cmq_l[upper.tri(cmq_l)] <- 0
      cmq <- if(!is.null(redim_mat)) redim_matrix(cmq_u+cmq_l, target_height=redim_mat, target_width=redim_mat, n_core=12) else cmq_u+cmq_l
      
      cmq <- melt(cmq)
      colnames(cmq) <- c("x","y","value")
      
      cmqp <- ggplot(cmq, aes(x = x, y = y, fill = value)) +
        geom_raster() +
        scale_fill_gradientn(colors = viridis::viridis_pal(option =palette)(quant)) +
        theme_void() + ggtitle(paste0("Upper : ",nm_u,"\n", "Lower : ",nm_l,"\n","chr: ", chr))
      if (params.rname != "") {
        cmqp <- cmqp + geom_rect(data = u_tad_left_dt,
                                 aes(xmin = x - 0.5, xmax = x - 0.495, ymin = y - 0.5, ymax = y + 0.5),
                                 fill = NA, color = "white") +
          geom_rect(data = u_tad_right_dt,
                    aes(xmin = x - 0.5, xmax = x - 0.495, ymin = y - 0.5, ymax = y + 0.5),
                    fill = NA, color = "white") + 
          geom_rect(data = l_tad_left_dt,
                    aes(xmin = x - 0.5, xmax = x - 0.495, ymin = y - 0.5, ymax = y + 0.5),
                    fill = NA, color = "white") +
          geom_rect(data = l_tad_right_dt,
                    aes(xmin = x - 0.5, xmax = x - 0.495, ymin = y - 0.5, ymax = y + 0.5),
                    fill = NA, color = "white")
      }
      gg[[paste0("ContactMatrices.Quantilized.",nm_u,"_vs_",nm_l,".Chromosome_", chr)]] <- cmqp
    }
  }
}


## PEAK DENSITY ----
Go <- loadranges("data/ChIPseq/ce11_tiled_5kb.bed",genome = "ce11")

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
  
  params.aname = paste0(sub("(.*)_GSE.+","\\1",item),"_GEO")
  params.rname = "tile_5kb"
  params = list(
    checkFile = FALSE,
    genome = 'ce11',
    downstream = 0,
    bed = paste0("data/ChIPseq/",item),
    upstream = 0,
    anchor = 'body',
    desc = 'Intersect Tile 5kb with Peakset of H3K4me3 from GEO - ',
    rname = 'H3K4me3_peaks'
  )  
  
  params$ignore.strand <- ifelse(is.null(params$ignore.strand),T,params$ignore.strand)
  
  
  if(!is.null(params$bed)){
    if(is.null(params$genome)) cat("params genome should be specified when params bed is used for intersect_ranges method")
    Go2 <- loadranges(params$bed, genome = "ce11")
  }  
  
  
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
  
  if(params$ignore.strand) {
    f <- GenomicRanges::findOverlaps(Go, Go2,ignore.strand=T)@from
  } else  f <- GenomicRanges::findOverlaps(Go, Go2)@from
  
  mcols(Go)[params.aname] <- F
  mcols(Go[f])[params.aname] <- T
}

### PLOT ----
params.aname = c('H3K4me3_GEO', 'H3K4me1_GEO', 'H3K9me3_GEO', 'hpl2_GEO', 'H3K9me2_GEO', 'H3K27ac_GEO', 'H3K27me3_GEO', 'LEM2_GEO', 'lin61_GEO')
params = list(
  res_kb = 5,
  chr_l = 'All',
  legend = TRUE,
  colors = c('#666666', '#132CEC', 'firebrick', '#D95F02', '#E7298A', '#66A61E', '#E6AB02', '#A6761D', '#1B9E77'),
  height = 500,
  width = 1500
)

res_kb <- params$res_kb
chr_l <- params$chr_l
my_mcols <- params.aname
Go_l <- split(Go, seqnames(Go))
ratio <- 1000/res_kb

if (chr_l == "All"){
  loop <- names(Go_l)
} else {
  loop <- chr_l
}

gg <- list()
for (chr in loop){
  dt <- as.data.table(mcols(Go_l[[chr]]))[,.SD, .SDcols = my_mcols]
  dt$chr_bin <- 1:nrow(dt)
  
  dt_p <- dt %>% tidyr::pivot_longer(cols=my_mcols, values_to = "is_here")
  dt_p$peaks <- ifelse(dt_p$is_here == T, dt_p$name, "EmptyBins")
  
  map_color <- NULL
  for (i in seq_along(my_mcols)) {
    if (!is.null(params[["colors"]]))
      map_color <- c(map_color,  params[["colors"]][i])
    else
      map_color <- c(map_color,  brewer.pal(length(my_mcols), "Set1")[i])
  }
  map_color <- c(map_color, "#ffffff")
  names(map_color) <- c(my_mcols, "EmptyBins")
  p <- ggplot(dt_p, aes(x = chr_bin, y = name, fill=peaks)) +  geom_tile(col=NA) +
    theme_minimal()+
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top",expand=c(-1,0)) +
    xlab(paste0('Peaks distribution over chromosome ',chr)) +
    scale_fill_manual(values = map_color, breaks = setdiff(unique(dt_p$peaks), "EmptyBins"), 
                      labels = setdiff(unique(dt_p$peaks), "EmptyBins")) +
    theme(title = element_text(size=12,face = "bold"),
          axis.text.x=element_text(angle = 90,size=10, vjust=0, hjust=0.6,face = "bold",margin=margin(20,0,0,0)),
          axis.text.y=element_text(size=10, margin=margin(0,0,0,20),face = "bold"),
          axis.title.x =element_text(size=18, vjust = 20, margin=margin(0,200,0,0),face = "bold"),
          axis.title.y =element_text(size=18,vjust = 0.8, margin=margin(0,0,0,0),face = "bold"),
          plot.margin = unit(c(2,1,2,1), "cm"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) 
  gg[[paste0("Peak_distribution.chromosome_",chr)]] <- p
}
ggsave(filename = "tst.png",plot = gg$Peak_distribution.chromosome_I,width = params$width/80,height = params$height/80)




## TAD CLUSTERING ----
### MCOLS Chr bin ----
Go <- loadranges("data/HiC/N2-old_merged.bwa_mem.25kb_domains.bed",genome = "ce11")
params.aname = "chr_bin"
params = list(bychr= T)
mcol <- NULL
if(params$bychr) {
  mcol <- unlist(lapply(split(Go,seqnames(Go)),function(chr){1:length(chr)}))
} else mcol <- 1:length(Go)
mcols(Go)[params.aname] <- NA
mcols(Go)[params.aname] <- mcol

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
params.aname = c('H3K4me3_peaks_cov', 'H3K4me1_peaks_cov', 'H3K9me3_peaks_cov', 'hpl2_peaks_cov', 'H3K9me2_peaks_cov', 'H3K27ac_peaks_cov', 'H3K27me3_peaks_cov', 'LEM2_peaks_cov', 'lin61_peaks_cov')
params = list(
  resolution = '25kb',
  height = 800,
  width = 800
)

output <- "results"
gp_res <- params$resolution
# chromosomes
gp_chr_l <- list(ALL=c("I","II","III","IV","V","X"))
# Graphical parameters
col_fun_sad <-  circlize::colorRamp2(seq(-1,1,length.out = 10), rev(brewer.pal(n = 10, name = "RdBu")))
col_fun_cov <-  circlize::colorRamp2(seq(0,1,length.out = 9), rev(brewer.pal(n = 9, name = "PuOr")))

for(chr_n in names(gp_chr_l)){
  chr_v <- gp_chr_l[[chr_n]]
  go_chr <- subset(Go,seqnames %in% chr_v)
  go_df <- adf(go_chr)
  toclust_bin_df <- subset(t(go_df),grepl("cov",colnames(go_df)))
  toclust_bin_df <- t(apply(toclust_bin_df,1,as.numeric))
  colnames(toclust_bin_df) <- 1:ncol(toclust_bin_df)
  # Normalize min max for plotting scale
  pdf(paste0(output,".TAD_clustering.Chromosome_",chr_n,".pdf"),width=10)
  res <- caCluster(toclust_bin_df, which="both", dim=2, opt.part=TRUE)
  resC <- caCluster(toclust_bin_df, which="cols", dim=2, opt.part=TRUE)
  resCfix <- caCluster(toclust_bin_df,part = 4, which="cols", dim=2, opt.part=TRUE)
  resR <- caCluster(toclust_bin_df, which="rows", dim=2, opt.part=TRUE)
  
  names(resR) <- paste0("class_",1:length(resR))
  names(resC) <- paste0("class_",1:length(resC))
  names(resCfix) <- paste0("class_",1:length(resCfix))
  col_ord <- ldply(resC,.id = "tad_class",function(X)data.frame(tad_idx_chr=names(X)))
  col_ord_fix <- ldply(resCfix,.id = "tad_class_fix",function(X)data.frame(tad_idx_chr=names(X)))
  row_ord <- ldply(resR,.id = "chip_class",function(X)data.frame(chip_ip=names(X)))
  
  toclust_bin_df <- toclust_bin_df[match(row_ord$chip_ip,rownames(toclust_bin_df)),]
  toclust_bin_df <- toclust_bin_df[,match(col_ord$tad_idx_chr,colnames(toclust_bin_df))]
  
  print(Heatmap(scale(toclust_bin_df), name = "mat",
                column_title = "CA clustering - TAD cluster fix : 4",
                cluster_rows = F,cluster_columns = F,border="black",
                row_split = row_ord$chip_class, column_split = col_ord_fix$tad_class_fix))
  dev.off()
}
