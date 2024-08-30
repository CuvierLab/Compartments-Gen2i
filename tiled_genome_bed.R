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

# RUN ===========================================================
# create ce11 gr
seqinfo_ce11 <- fai_lst$NCBI$ce11
seqinfo_ce11 <- GenomeInfoDb::dropSeqlevels(seqinfo_ce11,value = "MtDNA")
ce11_gr <- GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(seqinfo_ce11),
                       ranges = IRanges::IRanges(start=1,width=GenomeInfoDb::seqlengths(seqinfo_ce11)),
                       seqlengths = GenomeInfoDb::seqlengths(seqinfo_ce11)) 
for(res in c("1kb","5kb","10kb","25kb")){
  r_i <- as.numeric(gsub("(.*)kb","\\1",res))*1000
  ce11_t_gr <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(ce11_gr), tilewidth = r_i, cut.last.tile.in.chrom = TRUE)
  ce11_t_gr$bin <- 1:length(ce11_t_gr)
  ce11_t_gr$chr_bin <- unlist(lapply(split(ce11_t_gr,GenomeInfoDb::seqnames(ce11_t_gr)),
                                     function(chr){
                                       1:length(chr)
                                     }))
  
  rtracklayer::export.bed(ce11_t_gr,paste0("data/ChIPseq/ce11_tiled_",res,".bed"))
}


# create ce11 gr
ce11_gr <- GenomicRanges::GRanges(seqinfo_ce11)
ce11_noMT_gr <- GenomeInfoDb::dropSeqlevels(ce11_gr,value = c("MtDNA"),pruning.mode="coarse")
ce11_noMT_gr$name <- paste(as.character(GenomeInfoDb::seqnames(ce11_noMT_gr)), GenomicRanges::start(ce11_noMT_gr), GenomicRanges::end(ce11_noMT_gr), GenomicRanges::strand(ce11_noMT_gr),sep = "_")
rtracklayer::export.bed(ce11_noMT_gr,"data/ChIPseq/ce11_noMT.bed")

# ===============================================================