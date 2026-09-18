# HEADER ====================================================== =
# GEN2I
# NAR revision 2026
# R version : R version 4.4.0 (2024-04-24), tested on 4.5.3
# System : x86_64, linux-gnu
# ============================================================= =
#
# Supplementary Fig. S4F - PC2 tracks of the right arm of chromosome I,
# mutant over wild type.
#
# Same karyoplot as Supplementary Fig. S4D (see figure_S4D.R), with two
# differences: the eigenvector is the second principal component instead of the
# first, and the plot is zoomed on the right arm of the chromosome
# (I:10,120,013-15,072,433 in ce11) instead of the whole chromosome. The
# published panel is a further crop of that arm, from 11 to 15 Mbp.
#
# Wild type is drawn in red, the mutant in blue, both from the initial batch of
# Hi-C at 25 kb.
# ============================================================= =

# WORKING DIRECTORY ----
# Run the scripts from the root of this repository (where manuscript_lib.R sits):
# setwd("/path/to/Compartments-Gen2i")

# SOURCE ----
source("manuscript_lib.R")

# PARAMETERS ----
out_d      <- "output/figure_S4F"
conditions <- c("CEC4", "hpl2-old", "lin61-old", "hpl2-lin61-old", "met2-set25-set32-old")
zoom_bed   <- "data/ChIPseq/rightArm_ce11.bed"   # chromosome arms, see lift_over_chr_arms.R
chr        <- "I"
dir.create(out_d, showWarnings = FALSE, recursive = TRUE)

# RUN ----
## PC2 TRACKS ALONG THE RIGHT ARM OF CHROMOSOME I ----

Go      <- loadranges("data/ChIPseq/ce11_noMT.bed", genome = "ce11")
go_zoom <- loadranges(zoom_bed, genome = "ce11")

for (condition in conditions) {
  eigen_bw_l <- list(N2 = 'data/HiC/N2-old_merged.bwa_mem.25kb.pca2.bw')
  eigen_bw_l[[condition]] = paste0('data/HiC/', condition, '_merged.bwa_mem.25kb.pca2.bw')

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

  # SAVE ----
  png(file.path(out_d, sprintf("PC2.%s.karyo.%s.png", condition, chr)),
      width = width, height = height)
  kp <- karyoploteR::plotKaryotype(genome=subset(Go,seqnames==chr), plot.param=pp,
                                   zoom = subset(go_zoom,seqnames==chr),
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

  # SOURCE TABLE ----
  # The values behind the two tracks, over the plotted region, at the 25 kb
  # resolution of the eigenvectors.
  reg <- subset(go_zoom, seqnames == chr)
  sig <- lapply(names(bw_l), function(nm) {
    gr <- rtracklayer::import(bw_l[[nm]], which = reg)
    data.table(condition = nm, seqnames = as.vector(seqnames(gr)),
               start = start(gr), end = end(gr), score = gr$score)
  })
  fwrite(rbindlist(sig),
         file.path(out_d, sprintf("PC2.%s.karyo.%s.csv", condition, chr)))
}
