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

h5p <- list()
for (condition in c("CEC4",paste0(c("N2","hpl2","hpl2-lin61","I158A","lin61","met2-set25-set32"),"-old"))) {
  h5p[[condition]] <- paste0("data/HiC/",condition,"_merged.bwa_mem.25kb.obs_exp.g2i.h5")
}

params = list(
  h5path = h5p,
  na_offset = 20,
  diag_offset = 1,
  nb_clust = 4,
  legend = TRUE,
  height = 800,
  width = 800
)

# Only B compartments
Go <- Go[Go$eigen_pca1_N2.old<0]


Go_dt <- data.table(data.frame(Go))

## GRAPHICAL PARAMETERS
h5path <- params[["h5path"]]
na_offset <- ifelse(is.null(params[["na_offset"]]), 0, params[["na_offset"]])
diag_offset <- ifelse(is.null(params[["diag_offset"]]), 1, params[["diag_offset"]])
nb_clust <- params[["nb_clust"]]
chrs <- if(is.null(params[["chrs"]])) unique(Go_dt$seqnames) else params[["chrs"]]
# Load all cm
cm_l <- load_h5_matrix(h5path = h5path)
gg <- list()
dendro_dt <- NULL
for (chr in chrs) {
  dendro_chr_dt <- NULL
  for (condition in names(cm_l)) { # iterate through all contact matrices
    # Turn h5 to matrix
    cm_c <- as.matrix(as.matrix(cm_l[[condition]][[chr]]))
    max_bin <- dim(cm_c)[1] - na_offset # Remove NA's extend use in APA
    cm_c <- cm_c[(na_offset+1):max_bin, (na_offset+1):max_bin] # Remove NA's extend use in APA
    getKDiag(cm_c, -diag_offset:diag_offset) <- NA # Set offdiagonals to NAs
    colRow_sel <- Go_dt[seqnames == chr,chr_bin]
    dendro_chr_dt <- cbind(dendro_chr_dt, data.table::setnames(data.table(as.vector(cm_c[colRow_sel,colRow_sel])), condition))
  }
  dendro_dt <- rbind(dendro_dt,dendro_chr_dt)
}
hc <- t(dendro_dt) %>% scale %>% dist(method="euclidean") %>% hclust
hcdata <- dendro_data_k(hc, nb_clust)
cols    <- RColorBrewer::brewer.pal(nb_clust+1, "Dark2") # c("#a9a9a9", "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
gg[[paste0("Dendrogram.All_chrom.Dendogram")]] <- plot_ggdendro(hcdata,
                                                                direction   = "lr",
                                                                scale.color = cols,
                                                                expand.y    = 0.2) + theme_void()
gg
