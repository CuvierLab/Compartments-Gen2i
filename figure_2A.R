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
params = list(
  gp_all_xp = c('N2-old', 'I158A-old', 'hpl2-old', 'hpl2-lin61-old', 'lin61-old', 'met2-set25-set32-old', 'CEC4'),
  r = '25kb',
  norml = 'norm.KR',
  eigen_v = 'eigen_pca1_N2.old',
  legend = TRUE,
  height = 800,
  width = 800
)

r_i <- as.numeric(gsub("(.+)kb","\\1",params$r))*1000
gp_all_xp <- params$gp_all_xp
eigen_v <- params$eigen_v 
norml <- params$norml
res <- as.integer(sub("kb", "", params$r))

all_lst <- NULL
all_df <- NULL
gp_chr_l <- list(I="I",II="II",III="III",IV="IV",V="V",X="X")
adf <- as.data.frame
# Load pca for all xp 
Go$chrom <- seqnames(Go)
mol_df <-  data.table(data.frame(mcols(Go)))
mol_df$comp <- ifelse(mol_df[,get(eigen_v)>=0],"A","B")
# Load ce11 grtile so we have correspondence between pca and CM
for(xp in gp_all_xp) {
  # Load xp contact matrix in oe
  w_xp_r_cml <- readRDS(paste0("data/HiC/",xp,"_",params$r,".",norml,".cm.rds"))
  for(chr in names(gp_chr_l)){
    chr_v <- gp_chr_l[chr]
    # Subset mol
    mol_chr_row_df <- adf(subset(mol_df,chrom == chr)[,c("chr_bin","comp")])
    colnames(mol_chr_row_df) <- c("chr_bin","comp_row")
    mol_chr_col_df <- adf(subset(mol_df,chrom == chr)[,c("chr_bin","comp")])
    colnames(mol_chr_col_df) <- c("chr_bin","comp_col")
    
    # Load CM
    w_xp_r_chr_cm <- log10(dropCMbyChr(cml = w_xp_r_cml,chr = chr,offset = NULL,na2zero = T))
    
    # Create DF ROW/BIN from CM
    hic_df <-  reshape2::melt(w_xp_r_chr_cm)
    colnames(hic_df) <- c("rowbin","colbin","value")
    # Keep only upper triangle matrix
    hic_df <- subset(hic_df,colbin-rowbin >=0)
    hic_df$chr <- chr_v
    # Set compartment value for each bin
    hic_df <- dplyr::left_join(hic_df, mol_chr_row_df, 
                               by = c("rowbin" = "chr_bin")) 
    hic_df <- dplyr::left_join(hic_df, mol_chr_col_df, 
                               by = c("colbin" = "chr_bin")) 
    
    all_lst[[xp]] <- rbind(all_lst[[xp]],hic_df)
  }
}
all_df <- plyr::ldply(all_lst,function(X)data.frame(X))

all_df$genomic_dist <- (all_df$colbin - all_df$rowbin)*r_i
save_df <- all_df[,c("genomic_dist","value",".id","chr","comp_row","comp_col")]
save_df$chr <- unlist(save_df$chr)

if (!is.null(params$comp)) {
  save_itr_df <- subset(save_df,comp_row==params$comp & comp_col==params$comp)
} else {
  save_itr_df <- save_df
}

save_itr_df <- save_itr_df %>% group_by(genomic_dist, .id,chr) %>% summarize(m_val=mean(value,na.rm=T))
save_itr_df <- adf(save_itr_df)
save_itr_df$genomic_dist <- save_itr_df$genomic_dist + 1

rm(all_df,save_df)

col <- RColorBrewer::brewer.pal(length(gp_all_xp),"Set3")
col[2] <- "#BF8475"

gg <- NULL

setDT(save_itr_df)

save_itr_df[!is.finite(m_val), m_val:=0]
save_itr_df[, genomic_dist := genomic_dist/1000]
cut.d <- res 
save_itr_df <- save_itr_df[genomic_dist>cut.d,]

save_itr_df[, m_val := 1-(m_val/min(m_val, na.rm = T)), by = .(chr, .id)]
save_itr_df[genomic_dist < res+1, m_val := 1]

save_itr_df <- save_itr_df[genomic_dist < 6e3]

ggplot(save_itr_df, aes(x=genomic_dist, y=m_val,color=.id)) +
  geom_smooth(method = "loess",se = F,size=1.2,fill="grey", span = .01) +
  scale_color_manual(values=col) +
  scale_x_log10()+
  scale_y_log10()+
  annotation_logticks()  +
  # coord_trans(x="log10",y="log10")+
  theme_classic2() + ggtitle("All chromosomes - Limits 1e6")

# ============================================================= =

