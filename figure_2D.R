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

params = list(
  gp_xp = c('N2-old', 'CEC4', 'I158A-old', 'hpl2-lin61-old', 'hpl2-old', 'lin61-old', 'met2-set25-set32-old'),
  norm = 'obs_exp',
  legend = TRUE,
  height = 800,
  res = '25kb',
  width = 800
)



gp_chr_l <- list(chrom_ALL = c("I", "II", "III", "IV", "V", "X"))

res <- params$res
gp_xp <- params$gp_xp
norm <- params$norm

gg <- NULL

chrom_idx <- "chrom_ALL"

chrom.arms <- as.data.frame(import.bed("data/ChIPseq/chr_arms_ce11.bed"))
setDT(chrom.arms)
chrom.arms <- chrom.arms[, .(seqnames, start, end)]

chrom.arms[, ID := 1:.N, by = seqnames]

chrom.size <- chrom.arms[ID==3, end]
names(chrom.size) <- chrom.arms[ID==3, seqnames]

chrom.arms[ID == 1, start := 1]
chrom.arms[, tmp := min(end), by = seqnames]
chrom.arms[ID == 2, start := tmp, by = seqnames]

chrom.arms[ID == 3, end := 0, by = seqnames]
chrom.arms[, tmp := max(end), by = seqnames]
chrom.arms[ID == 3, end := start, by = seqnames]
chrom.arms[ID == 3, start := tmp, by = seqnames]

chrom.arms[ID == 1, class := "arm-left"]
chrom.arms[ID == 2, class := "center"]
chrom.arms[ID == 3, class := "arm-right"]

chrom.arms[, ID := NULL]
chrom.arms[, tmp := NULL]

# make tile genome to keep bin type 
tile.gr <- tileGenome(chrom.size, tilewidth = as.integer(sub("kb", "", res)) * 1000, cut.last.tile.in.chrom = T)
tile.gr <- adf(tile.gr)
setDT(tile.gr)

# add chr_idx and all_idx
tile.gr[, all_idx := 1:.N]
tile.gr[, chr_idx := 1:.N, by = seqnames]

setkey(tile.gr, seqnames, start, end)
setkey(chrom.arms, seqnames, start, end)

# overlap tile and arm anno
tile.gr <- foverlaps(tile.gr, chrom.arms)


# solve duplicate/ambiguous with size of overlap
tile.gr[, ov := min(end, i.end) - max(start, i.start), by = seq_len(nrow(tile.gr))]
tile.gr[order(-ov), n := 1:.N, by = all_idx]
tile.gr <- tile.gr[n == 1]

dat <- NULL
for(xp in gp_xp){
  tadbin_xp.chr.l <- readRDS(paste0("data/HiC/", xp, "_", res,".",norm,".cm.rds"))
  for(i in gp_chr_l[[chrom_idx]]){
    cm <- as.matrix(as.matrix(tadbin_xp.chr.l[[i]]))
    
    # remove outliers column
    Q <- quantile(rowSums(cm, na.rm = T), probs=c(.05, .95), na.rm = T)
    iqr <- IQR(rowSums(cm, na.rm = T), na.rm = T)
    to_remove <- which(!rowSums(cm) > (Q[1] - 1.5*iqr) & rowSums(cm) < (Q[2]+1.5*iqr))
    tile.gr[seqnames == i & chr_idx %in% to_remove, class := "removed"]
    cm[to_remove,] <- NA
    cm[,to_remove] <- NA
    
    getKDiag(cm, -2:2) <- NA
    
    
    # remove outliers point
    Q <- quantile(cm, probs=c(.05, .95), na.rm = T)
    iqr <- IQR(cm, na.rm = T)
    to_remove <- which(cm > (Q[2]+1.5*iqr))
    to_remove <- c(to_remove, which(cm < (Q[1]-1.5*iqr)))
    
    cm[to_remove] <- NA
    
    cm[cm == 0] <- NA
    cm <- log2(cm)
    # pheatmap::pheatmap(cm, cluster_rows = F, cluster_cols = F)
    
    idx.arm.l <- tile.gr[seqnames == i][class == "arm-left", chr_idx]
    idx.arm.r <- tile.gr[seqnames == i][class == "arm-right", chr_idx]
    idx.center <- tile.gr[seqnames == i][class == "center", chr_idx]
    
    inter.arm <- as.matrix(cm[idx.arm.l, idx.arm.r])
    
    intra.arm.l <- as.matrix(as.matrix(cm[idx.arm.l, idx.arm.l]))
    getKDiag(intra.arm.l, -2:2) <- NA
    intra.arm.l[lower.tri(intra.arm.l)] <- NA
    
    intra.arm.r <- as.matrix(as.matrix(cm[idx.arm.r, idx.arm.r]))
    getKDiag(intra.arm.r, -2:2) <- NA
    intra.arm.r[lower.tri(intra.arm.r)] <- NA
    
    arm.center.l <- as.matrix(as.matrix(cm[idx.arm.l, idx.center]))
    arm.center.r <- as.matrix(as.matrix(cm[idx.arm.r, idx.center]))
    
    
    
    dat <- rbind(dat, 
                 data.table(value = as.vector(inter.arm), class = "inter arm", xp = xp),
                 data.table(value = as.vector(intra.arm.l), class = "intra arm", xp = xp),
                 data.table(value = as.vector(intra.arm.r), class = "intra arm", xp = xp),
                 data.table(value = as.vector(arm.center.l), class = "arm to center", xp = xp),
                 data.table(value = as.vector(arm.center.r), class = "intra arm", xp = xp)
    )
  }
}


dat[, value := scale(value), by = .(xp)]
dat <- dat[is.finite(value)]
dat[, ecd := ecdf(value)(value), by = .(class, xp)]

setnames(dat, "value", "zscore")

for (clss in unique(dat$class)) {
  
  sub_dt <- subset(dat, class == clss)
  gg[[paste0("Lineplot.",clss,".",chrom_idx)]] <- ggplot(sub_dt, aes(x = zscore, y = ecd, group = xp, colour = xp)) +
    geom_line()+
    coord_cartesian(xlim = c(-4, 4))+
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray")+
    scale_color_brewer(palette = "Set1")
}


