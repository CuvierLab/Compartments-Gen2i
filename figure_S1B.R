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
## SADDLEPLOT ----
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
  genome = 'ce11'
)
# convert to Granges

cat("bigwig_signal_tile:... ")
# save.image("dev_l.RData")

Go <- bw_tile(bw = params$bw_p,
              bed = Go,
              genome = config.genome,
              name = params.aname)


params.aname = "eigen_pca1_N2"
params = list(
  bw_p = 'data/HiC/N2_merged.bwa_mem.25kb.pca1.bw',
  genome = 'ce11'
)
# convert to Granges

cat("bigwig_signal_tile:... ")
# save.image("dev_l.RData")

Go <- bw_tile(bw = params$bw_p,
              bed = Go,
              genome = config.genome,
              name = params.aname)


### PLOT ----
Go <- go2dt(subset(Go,seqnames(Go)=="I"))
params = list(
  x = 'eigen_pca1_N2',
  y = 'eigen_pca1_N2.old',
  legend = TRUE,
  height = 800,
  width = 800
)

sgg <- ggplot(Go,
              aes_string(x=params$x,
                         y=params$y))
sgg <- sgg + geom_point()
gg2 <- sgg + stat_regline_equation() + stat_smooth(method = "lm") + stat_cor(label.x.npc = "center")

sgg <- sgg + geom_abline(slope = 1, intercept = 0)

gg <- list()
gg[["Scatterplot.Classic"]] <- sgg 
gg[["Scatterplot.Regression_line"]] <- gg2
