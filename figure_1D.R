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
## CONTACT MATRICES HEATMAP ----

Go <- loadranges("data/ChIPseq/ce11_noMT.bed",genome = "ce11")

for (condition in c("CEC4",paste0(c("hpl2","hpl2-lin61","I158A","lin61","met2-set25-set32"),"-old"))) {
  eigen_bw_l <- list(N2 = 'data/HiC/N2-old_merged.bwa_mem.25kb.pca1.bw')
  eigen_bw_l[[condition]] = paste0('data/HiC/',condition,'_merged.bwa_mem.25kb.pca1.bw')
  
  params = list(
    ymin = -0.2,
    ymax = 0.2,
    bw_l = eigen_bw_l,
    chr_v = c('I', 'II', 'III', 'IV', 'V', 'X'),
    title = '',
    legend = TRUE,
    height = 500,
    width = 1500
  )
  
  pp <- karyoploteR::getDefaultPlotParams(plot.type=1)
  pp$ideogramheight <- 0.1
  pp$leftmargin <- 0.1
  pp$rightmargin <- 0.1
  pp$topmargin <- 20
  pp$bottommargin <- 10
  pp$data1inmargin <- 10
  
  width <- params$width
  height <- params$height
  
  ymin <- params$ymin
  ymax <- params$ymax
  
  bw_l <- params$bw_l
  chr_v <- params$chr
  
  mycol <- RColorBrewer::brewer.pal(length(bw_l),"Set1")
  
  chr <- "I"
  output <- paste0("./",condition)
  fn <- paste0("karyo.",chr)
  png(paste0(output,".",fn,".png"), width = width, height = height)
  kp <- karyoploteR::plotKaryotype(genome=subset(Go,seqnames==chr), plot.param=pp,
                                   main = paste0("BigWig signal snapshot \n Zoom on region : " ,chr),
                                   labels.plotter = NULL)
  karyoploteR::kpAddBaseNumbers(kp,tick.len = 2,cex = 1.3,minor.tick.dist = 1e6,minor.tick.len = 1,tick.dist = 2e6, add.units = TRUE)
  lab <- glue::glue_collapse(names(bw_l),"\n")
  # karyoploteR::kpAddLabels(kp, labels = lab,label.margin = 0.03)
  toEval <- NULL
  for (i in 1:length(bw_l)) {
    bw_p <- bw_l[i]
    toEval <- paste0(toEval,
                     "\n",
                     paste0('karyoploteR::kpPlotBigWig(karyoplot = kp,data = "',bw_p,'", ymin = ',ymin,', ymax = ',ymax,', col = NA, border = "',mycol[i],'",lwd=3)'))
    
  }
  eval(parse(text=toEval))
  karyoploteR::kpAxis(kp,ymin = ymin,ymax = ymax)
  karyoploteR::kpAbline(kp, h=0, ymin=ymin,ymax=ymax, lty=2,lwd=3, col="#666666")
  legend(x = "topright", fill = mycol[1:length(bw_l)], legend = names(bw_l))
  dev.off()
}