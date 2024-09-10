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


Go <- loadranges("data/ChIPseq/ce11_noMT.bed",genome = "ce11")
params.aname = ''
params = list(
  resolution = '25kb',
  legend = TRUE,
  height = 800,
  width = 800
)

res <- params[["resolution"]]
outdir <- gsub(".done","",snakemake@output[[1]])

gp_xp_rnaseq <- c("hpl2","hpl2_lin61","met2Set25Set32","lin61","hpl2I158A")
gp_gtf_path <- "data/RNAseq/matrecap_genes_ce11.gtf"
genes_mrc_df <- data.table::fread(file = gp_gtf_path)

pp <- NULL

outdir <- "Figure_S4B"

for (xp in gp_xp_rnaseq) {
  
  # Exctract LFC column information
  sub_col <- colnames(genes_mrc_df)[grepl(paste0("WT_vs_",xp,"_(FDR|pval|log2FoldChange)"),colnames(genes_mrc_df))]
  sel_col <- c("V1","V2","V3",sub_col)
  deg_mrc <- genes_mrc_df[,.SD, .SDcols = colnames(genes_mrc_df) %in% sel_col]
  colnames(deg_mrc) <- c("seqnames","start","end",gsub(paste0("WT_vs_",xp,"_"),"\\1",sub_col))
  deg_gr <- GRanges(deg_mrc) %>%
    GenomeInfoDb::dropSeqlevels("MtDNA",pruning.mode = "coarse") 
  
  urg_gr <- GRanges(data.frame(deg_gr) %>% filter(log2FoldChange < 0 & pval < 5e-2 & FDR < 0.1))
  drg_gr <- GRanges(data.frame(deg_gr) %>% filter(log2FoldChange >= 0 & pval < 5e-2 & FDR < 0.1))
  ce11_tiled_gr <- tileGenome(seqlengths = seqlengths(Go),tilewidth = as.numeric(sub("kb","",res))*1000,cut.last.tile.in.chrom = T)
  
  chr <- "I"
    # browser()
    ce11_tiled_chr_gr <- subset(ce11_tiled_gr,seqnames==chr)
    ur_tiled_gr <- ce11_tiled_chr_gr
    dr_tiled_gr <- ce11_tiled_chr_gr
    ur_tiled_gr$dens <- as.numeric(cut(countOverlaps(ur_tiled_gr,urg_gr),5))
    ur_tiled_gr$col <- brewer.pal(5,"Blues")[ur_tiled_gr$dens]
    dr_tiled_gr$dens <- as.numeric(cut(countOverlaps(dr_tiled_gr,drg_gr),5))
    dr_tiled_gr$col <- brewer.pal(5,"Reds")[dr_tiled_gr$dens]
    
    png(paste0(outdir,".KaryoPlot.Condition_",xp,".Chromosome.",chr,".png"))
    # KARYPLOT
    kp <- plotKaryotype(genome=Go,chromosomes = chr,plot.type = 1,plot.param=pp,
                        main = paste0("Chromosome ",chr),labels.plotter = NULL) # force removal of mtdna)
    kpAddBaseNumbers(kp,tick.len = 2,minor.tick.dist = 1e6,minor.tick.len = 1,tick.dist = 2e6, add.units = TRUE)
    
    at <- autotrack(1,1,margin=0.1)
    # Add gene density
    r0n <- at$r0+(at$r1 - at$r0)/8
    r00n <- at$r0+(r0n - at$r0)/2
    kpPlotRegions(kp, ur_tiled_gr, r0=r00n, r1=r0n,col=ur_tiled_gr$col,num.layers = 1)
    kpPlotRegions(kp, dr_tiled_gr,  r0=r00n, r1=at$r0,col=dr_tiled_gr$col,num.layers = 1)
    kpAddLabels(kp, labels = paste0("UPREG"),r0=r00n, r1=r0n, label.margin = 0.06,cex=0.4, srt=0, pos=3, col=rgb(51/255, 133/255, 1))
    kpAddLabels(kp, labels = paste0("DOWNREG"),r0=r00n, r1=at$r0, label.margin = 0.06,cex=0.4, srt=0, pos=3,col=rgb(1, 0, 0))
    dev.off()
}