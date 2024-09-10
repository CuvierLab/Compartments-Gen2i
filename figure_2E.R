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

Go <- bw_tile(bw = params$bw_p,
              bed = Go,
              genome = config.genome,
              name = params.aname)


for (condition in paste0(c("hpl2","lin61","hpl2-lin61","met2-set25-set32","I158A"),"-old")) {
  
  dot_cond = sub("-","\\.",condition)
  params.aname = paste0("eigen_pca1_",dot_cond)
  params = list(
    bw_p = paste0('data/HiC/',condition,'_merged.bwa_mem.25kb.pca1.bw'),
    genome = 'ce11'
  )
  
  Go <- bw_tile(bw = params$bw_p,
                bed = Go,
                genome = config.genome,
                name = params.aname)
  
  #### DELTA EIGEN ----
  params.aname = paste0("delta_eigen_pca1_",dot_cond,"_N2.old")
  params.process = 
    params = list(
      operation = 'subtract'
    )
  params.dependencies = c('eigen_pca1_' %+% dot_cond, 'eigen_pca1_N2.old')
  mcols(Go)[params.aname] <- NA
  method <- ifelse(!is.null(params$method),params$method, params$operation)
  mcols(Go)[params.aname] <- Gooperation(unlist(mcols(Go)[[params.dependencies[1]]]),
                                         unlist(mcols(Go)[[params.dependencies[2]]]),
                                         method = method)
}



### GET TAD CLASS ON TAD N2 OLD ----
Go3 <- loadranges("data/HiC/N2-old_merged.bwa_mem.25kb_domains.bed",genome = "ce11")
params.aname = "chr_bin"
params = list(bychr= T)
mcol <- NULL
if(params$bychr) {
  mcol <- unlist(lapply(split(Go3,seqnames(Go3)),function(chr){1:length(chr)}))
} else mcol <- 1:length(Go3)
mcols(Go3)[params.aname] <- NA
mcols(Go3)[params.aname] <- mcol

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
  mcols(Go3)[[params.aname]] <- getCoverageDensity(Go3, Go2, bychr=params$bychr)
}

params.aname = c('H3K4me3_peaks_cov', 'H3K4me1_peaks_cov', 'H3K9me3_peaks_cov', 'hpl2_peaks_cov', 'H3K9me2_peaks_cov', 'H3K27ac_peaks_cov', 'H3K27me3_peaks_cov', 'LEM2_peaks_cov', 'lin61_peaks_cov')
params = list(
  resolution = '25kb',
  height = 800,
  width = 800
)

output <- "results"
gp_res <- params$resolution
# chromosomes
chr_v <- c("I","II","III","IV","V","X")
# Graphical parameters
col_fun_sad <-  circlize::colorRamp2(seq(-1,1,length.out = 10), rev(brewer.pal(n = 10, name = "RdBu")))
col_fun_cov <-  circlize::colorRamp2(seq(0,1,length.out = 9), rev(brewer.pal(n = 9, name = "PuOr")))

go_chr <- subset(Go3,seqnames %in% chr_v)
go_df <- adf(go_chr)
toclust_bin_df <- subset(t(go_df),grepl("cov",colnames(go_df)))
toclust_bin_df <- t(apply(toclust_bin_df,1,as.numeric))
colnames(toclust_bin_df) <- 1:ncol(toclust_bin_df)
# Normalize min max for plotting scale
resC <- caCluster(toclust_bin_df, which="cols", dim=2, opt.part=TRUE)
resCfix <- caCluster(toclust_bin_df,part = 4, which="cols", dim=2, opt.part=TRUE)
names(resC) <- paste0("class_",1:length(resC))
names(resCfix) <- paste0("class_",1:length(resCfix))
col_ord <- ldply(resC,.id = "tad_class",function(X)data.frame(tad_idx_chr=names(X)))
col_ord_fix <- ldply(resCfix,.id = "tad_class_fix",function(X)data.frame(tad_idx_chr=names(X)))


mcols(Go3)["tad_classes_chipdensity_4classes_GEO"] <- NA
GenomicRanges::mcols(Go3[as.numeric(col_ord$tad_idx_chr)])["tad_classes_chipdensity_4classes_GEO"] <- as.character(col_ord_fix$tad_class)



### INHERIT DELTA EIGEN ON TAD SCALE ----

for (condition in paste0(c("hpl2","lin61","hpl2-lin61","met2-set25-set32","I158A"),"-old")) {
  
  dot_cond = sub("-","\\.",condition)
  params.aname = paste0("delta_eigen_pca1_",dot_cond,"_N2.old")
  params = list(
    aname = paste0("delta_eigen_pca1_",dot_cond,"_N2.old"),
    fun = 'mean',
    fun_opts = 'na.rm=T'
  )
  
  Go3 <- aggregate_ranges(Go3, Go, out_mc=params.aname,  subject_mc = params$aname, fun_agr = params$fun,fol_opts = params$fol_opts, fun_opts=params$fun_opts)
}



### PLOT ----

Go <- subset(Go3, tad_classes_chipdensity_4classes_GEO %in% c("class_3","class_4"))


params.aname = c('delta_eigen_pca1_CEC4_N2.old', 'delta_eigen_pca1_hpl2.lin61.old_N2.old', 'delta_eigen_pca1_hpl2.old_N2.old', 'delta_eigen_pca1_I158A.old_N2.old', 'delta_eigen_pca1_lin61.old_N2.old', 'delta_eigen_pca1_met2.set25.set32.old_N2.old')
params = list(
  res = 25000,
  apa_size = 41.0,
  d_bin_min = 2.0,
  d_bin_max = 1000000.0,
  h5path = list(`hpl2-lin61-old` = 'data/HiC/hpl2-lin61-old_merged.bwa_mem.25kb.obs_exp.g2i.h5', 
                `hpl2-old` = 'data/HiC/hpl2-old_merged.bwa_mem.25kb.obs_exp.g2i.h5', 
                `I158A-old` = 'data/HiC/I158A-old_merged.bwa_mem.25kb.obs_exp.g2i.h5', 
                `lin61-old` = 'data/HiC/lin61-old_merged.bwa_mem.25kb.obs_exp.g2i.h5', 
                `met2-set25-set32-old` = 'data/HiC/met2-set25-set32-old_merged.bwa_mem.25kb.obs_exp.g2i.h5', 
                `N2-old` = 'data/HiC/N2-old_merged.bwa_mem.25kb.obs_exp.g2i.h5'),
  cons_nm = 'intra',
  save_dt = FALSE,
  orientation = '5p3p',
  forceSym = FALSE,
  offset = 1.0,
  title = '',
  legend = TRUE,
  threads = 1,
  height = 800,
  width = 800
)

Go_dt <- go2dt(Go)
setkey(Go_dt, name)


#### 3D PARAMETERS
params$orientation <- ifelse(is.null(params$orientation), "both", params$orientation)
params$forceSym <-  ifelse(is.null(params$forceSym), F, params$forceSym)
params$replace <-  ifelse(is.null(params$replace), F, params$replace)
params$nb_iter <-  ifelse(is.null(params$nb_iter), 1, params$nb_iter)
params$metrics_sel <-  ifelse(is.null(params$metrics_sel), c("standard_metrics","diffusion_metrics"), c("standard_metrics"))
params$d_bin_min <- max(params$offset,params$d_bin_min)

if(is.null(params$offset)) params$offset <- 0

#### APA PARAMETERS
mask <- matrix(1:(params$apa_size*params$apa_size), nrow=params$apa_size, byrow=F)
####GRAPHICAL PARAMETERS
params$plot_stat <- ifelse(is.null(params$plot_stat), T, params$plot_stat)
width <- ifelse(is.null(params$width), 800,  as.numeric(params$width))
height <- ifelse(is.null(params$height), 800,  as.numeric(params$height))

# Set anchor and bait as TAD border left and right
l_bait_go_tile <- list()
l_anchor_go_tile <- list()

l_bait_go_tile[["TADborders"]] <- IRanges::resize(Go,fix = "start",1)
l_anchor_go_tile[["TADborders"]] <- IRanges::resize(Go,fix = "end",1)
# get chr_bin equivalent to retrieve hic bins positions for start and end of each TADs
l_bait_go_tile[["TADborders"]]$chr_bin <- ceiling(start(Go)/params$res)
l_anchor_go_tile[["TADborders"]]$chr_bin <- ceiling(end(Go)/params$res)


l_cm <- load_h5_matrix(params$h5path)


#### EXTRACT CONTACT FREQUENCIES FOR EACH COUPLES/EACH MATRICES ----
l_flat_dt <- NULL
nb_cpl_dt <- NULL
dst_cpl_dt <- NULL
out <- list()

# Extract pairwise candidates
bait.nm <- names(l_bait_go_tile)[1]
anchor.nm <- names(l_anchor_go_tile)[1]
pw_dt <- pairwiseInteractionDT(anc_go = l_anchor_go_tile[[anchor.nm]],
                               bait_go =  l_bait_go_tile[[bait.nm]],
                               hic_res = params$res,
                               apa_size =  params$apa_size,
                               d_min = params$d_bin_min,
                               d_max = params$d_bin_max,
                               forceSym = params$forceSym,
                               cons_nm = params$cons_nm,
                               sym = params$orientation)
pw_dt[,keep:=paste0(bait_chrom_bin,"_",bait_chrom,"_",anc_chrom_bin,"_",anc_chrom)]
tokeep <- paste0(l_bait_go_tile[[1]]$chr_bin,"_",seqnames(l_bait_go_tile[[1]]),"_",l_anchor_go_tile[[1]]$chr_bin,"_",seqnames(l_anchor_go_tile[[1]]))
pw_dt <- pw_dt[keep %in% tokeep]
pw_dt$group <- 1

for(cm.nm in names(l_cm)){
  l_flat_dt[[paste0(bait.nm,"\nVS\n",anchor.nm)]][[cm.nm]] <- dumpFlatInteractionsMatrices(pairwise_dt = pw_dt, cm=l_cm[[cm.nm]], intra = T, n = params$apa_size)
}      

#### APA TO FLAT & USABLE MATRICES ----
flat_dt <- aggreg_lflat(l_flat_dt)
remove_diag(flat_dt, params$offset, params$apa_size)
coord_l <- get_apa_metrics_coordinate(params$apa_size, mask,ctloffset = 6)
setApaMetrics(flat_dt, coord_l)

#### APA MATRICES NORMALIZATIONS ----

flat_nopix_dt <- flat_dt[, .SD, .SDcols = !patterns("pix")]
setkey(flat_nopix_dt, bait_name)
scat_dt <- Go_dt[flat_nopix_dt]
scat_dt[,cm:=gsub("\\-",".",cm)]  
tosel <- colnames(scat_dt)[colnames(scat_dt) %like% "CP|BLS"]
scat_dt[,(tosel):=lapply(.SD, function(x) ifelse(x==0,NA,x)),.SDcols = tosel]

scat_n2_dt <- subset(scat_dt, cm %in% c("N2","N2.old"))
l_gg <- list()
for (cond in unique(scat_dt$cm)) {
  if (cond %like% "N2") next()
  scat_s_dt <- subset(scat_dt, cm ==cond)
  scat_s_dt$deltaTADstrength <-  scat_s_dt$BLS/scat_n2_dt$BLS
  eigen_col <- grep(paste0("delta_eigen_pca1_",cond,"_N2"), colnames(scat_s_dt),value = T)
  l_gg[[paste0("Scatterplot.raw.deltaEigen_vs_deltaTADstrength.",gsub("\\.","-",cond))]] <- ggplot(scat_s_dt, aes(x = !!sym(eigen_col), y = deltaTADstrength)) +  geom_point()+ stat_regline_equation() + stat_smooth(method = "lm") + stat_cor(label.x.npc = "center") + theme_classic()
}

