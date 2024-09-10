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
# convert to Granges

cat("bigwig_signal_tile:... ")
# save.image("dev_l.RData")

Go <- bw_tile(bw = params$bw_p,
              bed = Go,
              genome = config.genome,
              name = params.aname)

### PLOT ----
gg <- list()

for (condition in c("CEC4",paste0(c("hpl2","hpl2-lin61","I158A","lin61","met2-set25-set32"),"-old"))) {
  
  h5p <- list()
  h5p[[condition]] <- paste0("data/HiC/",condition,"_vs_N2-old.bwa_mem.25kb.norm.KR.g2i.h5")
  params = list(
    anchor = 'eigen_pca1_N2.old',
    bait = 'eigen_pca1_N2.old',
    bait_binsize = 50,
    anchor_binsize = 50,
    h5path = h5p,
    scaled = TRUE,
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
  for(h5.nm in names(h5path)){
    # get h5 name info by h5path
    cm_info <- rhdf5::h5ls(h5path[[h5.nm]])
    if (!is.null(chr_v))
      l_cm[[h5.nm]] <- sapply(cm_info$name[cm_info$name %in% chr_v], function(nm) { HDF5Array::HDF5Array(h5path[[h5.nm]], nm)})
    else # Keep all chrom
      l_cm[[h5.nm]] <- sapply(cm_info$name, function(nm) { HDF5Array::HDF5Array(h5path[[h5.nm]], nm)})
  }
  
  sadmat_l <- list()
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
    m$counts <- g2ismk::rangeMinMax(m$counts)
    gg_center <- ggplot2::ggplot(m, ggplot2::aes(x=Anchor, y=Bait, fill=counts)) +
      ggplot2::geom_tile() + theme(axis.text.x = element_text(angle=90),axis.title.y = element_blank(),axis.title.x = element_blank())+
      scale_x_discrete(position="top") +   scale_fill_gradient2(mid = "white",high = "firebrick", low = "dodgerblue") 
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
}
