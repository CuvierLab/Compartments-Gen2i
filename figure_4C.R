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

### LOAD BEDFILES ----
bed_l <- list(hpl2 = "hpl2_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
              H3K9me3 = "H3K9me3_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
              H3K9me2 = "H3K9me2_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
              lin61 = "lin61_N2.bwa_aln.rmdup._peaks.narrowPeak_optimal.bed",
              hpl2_lin61 = "hpl2_lin61_N2_sharp_intersect_peakset.bed")

Go <- NULL
for (bed in names(bed_l)) {
  gr <- loadranges(paste0("data/ChIPseq/",bed_l[[bed]]), "ce11")
  gr$rname <- bed
  Go <- c(Go,gr)
}
Go <- do.call("c",Go)


params.rname = names(bed_l)
params = list(
  ranges_ref = 'hpl2_lin61',
  ranges_prof = c('hpl2', 'lin61', 'H3K9me2', 'H3K9me3'),
  res_fix = 'center',
  win_lim = 500000.0,
  binsize = 1000.0,
  width = 1500,
  legend = TRUE
)


### PLOT ----

## FUNCTION PARAMETERS ----
c_ref <- params$ranges_ref
ranges_prof <- params$ranges_prof
win_lim <- params$win_lim
cond_v <- params$cond_v
binsize <- params$binsize
fix <- params$res_fix 
binnumber <- 2*win_lim/binsize + 1
x_seq <- -((binnumber-1)/2):((binnumber-1)/2)

Go$condition <- "N2"
plot_dt <- NULL

## Set GOs ----
Go_ref <- GenomicRanges::trim(GenomicRanges::resize(subset(Go, rname == c_ref),width = binsize, fix = fix))
Go_ref <- Go_ref[width(Go_ref) == binsize]
Go_prof <- subset(Go, rname %in% ranges_prof & condition == unique(Go_ref$condition))

## Set GOs ----
Go_ref_xt <- trim(Go_ref + win_lim)
Go_ref_xt <- Go_ref_xt[width(Go_ref_xt)== binsize + 2*win_lim] # Select complete ranges
Go_ref_xt_tiled <- unlist(GenomicRanges::tile(Go_ref_xt, width = binsize))

all_dt <- NULL
norm_dt <- NULL # Normalize by total peak base spanning
for (ip in unique(Go_prof$rname)){
  Go_prof_s <- subset(Go_prof,rname == ip)
  norm_dt <- rbind(norm_dt,data.table(name = ip, norm = sum(width(Go_prof_s))/1e6))
  col <- countOverlaps(Go_ref_xt_tiled,Go_prof_s)
  all_dt <- rbind(all_dt,
                  data.table(bin = rep(x_seq, length(Go_ref_xt)), value = col,name = ip)
  )
}

val_dt <- all_dt[,.(svalue=sum(value)),by = .(bin,name)]
data.table::setkey(val_dt, name)
data.table::setkey(norm_dt, name)
sub_p_dt <- norm_dt[val_dt]
sub_p_dt[,normalized_density:=svalue/norm]
sub_p_dt[,facetBy:=sub(".+cond_(.+)","\\1",name)]

plot_dt <- rbind(plot_dt,sub_p_dt)
gg <- list()
plot_dt[,ip:=gsub("(.+)_WF.+","\\1",name)]
gg[["Density_profile.All_features"]] <- ggplot(plot_dt, aes(x = bin, y = normalized_density, group = ip)) + geom_step(aes(col = ip)) + xlim(-50,50)
