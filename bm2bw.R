# WORKING DIRECTORY ----
setwd(dir = "config/src/R/github/")

# SOURCE ----
source("manuscript_lib.R")


# RUN ----
for (file in list.files("data/HiC/",pattern = "tad_score.bm")) {
    
    o_n <- sub("bm$","bw",file)
    bm_gr <- bm2GR(bm_p = "data/HiC/"%+%file, genome = "ce11", dropChr =  "MtDNA")
    bm_dt <- data.table::data.table(data.frame(bm_gr))
    
    cols <- colnames(GenomicRanges::mcols(bm_gr))
    # Center reduce before scaling 
    bm_dt[, (cols) := lapply(.SD,function(x) (x - mean(x, na.rm = T))/sd(x,na.rm = T)), by = seqnames, .SDcols = cols]
    bm_dt[, (cols) := lapply(.SD,function(x) rangeMinMax(x, -1, 1)), by = seqnames, .SDcols = cols]
    
    bm_bw <- GenomicRanges::GRanges(bm_dt[,score:= apply(.SD, 1, function(x) mean(x,na.rm = T)), .SDcols = cols][,(cols):=NULL])
    
    bm_bw <- setSeqinfo(bm_bw, "ce11")
    
    rtracklayer::export.bw(object = bm_bw, con = "data/HiC/"%+%o_n)
}

