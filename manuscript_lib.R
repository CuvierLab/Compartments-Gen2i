# HEADER ====================================================== =
# GEN2I
# Tue Jun 18 18:22:27 2024
# R version : R version 4.4.0 (2024-04-24)
# System : x86_64, linux-gnu
# ============================================================= =


# LIBRARIES ----
library(crayon)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggpubr)
library(karyoploteR)
library(rtracklayer)
library(corrr)
library(GGally)
library(ComplexHeatmap)
library(plyr)
library(CAinterprTools)
library(RColorBrewer)
library(ggVennDiagram)
library(ggdendro)

# PLOT ----
dendro_data_k <- function(hc, k) {
  
  hcdata    <-  ggdendro::dendro_data(hc, type = "rectangle")
  seg       <-  hcdata$segments
  labclust  <-  cutree(hc, k)[hc$order]
  segclust  <-  rep(0L, nrow(seg))
  heights   <-  sort(hc$height, decreasing = TRUE)
  height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)
  
  for (i in 1:k) {
    xi      <-  hcdata$labels$x[labclust == i]
    idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3    <-  seg$yend < height
    idx     <-  idx1 & idx2 & idx3
    segclust[idx] <- i
  }
  
  idx                    <-  which(segclust == 0L)
  segclust[idx]          <-  segclust[idx + 1L]
  hcdata$segments$clust  <-  segclust
  hcdata$segments$line   <-  as.integer(segclust < 1L)
  hcdata$labels$clust    <-  labclust
  
  hcdata
}

set_labels_params <- function(nbLabels,
                              direction = c("tb", "bt", "lr", "rl"),
                              fan       = FALSE) {
  if (fan) {
    angle       <-  360 / nbLabels * 1:nbLabels + 90
    idx         <-  angle >= 90 & angle <= 270
    angle[idx]  <-  angle[idx] + 180
    hjust       <-  rep(0, nbLabels)
    hjust[idx]  <-  1
  } else {
    angle       <-  rep(0, nbLabels)
    hjust       <-  0
    if (direction %in% c("tb", "bt")) { angle <- angle + 45 }
    if (direction %in% c("tb", "rl")) { hjust <- 1 }
  }
  list(angle = angle, hjust = hjust, vjust = 0.5)
}
plot_ggdendro <- function(hcdata,
                          direction   = c("lr", "rl", "tb", "bt"),
                          fan         = FALSE,
                          scale.color = NULL,
                          branch.size = 1,
                          label.size  = 3,
                          nudge.label = 0.01,
                          expand.y    = 0.1) {
  
  direction <- match.arg(direction) # if fan = FALSE
  ybreaks   <- pretty(ggdendro::segment(hcdata)$y, n = 5)
  ymax      <- max(ggdendro::segment(hcdata)$y)
  
  ## branches
  p <- ggplot() +
    geom_segment(data         =  ggdendro::segment(hcdata),
                 aes(x        =  x,
                     y        =  y,
                     xend     =  xend,
                     yend     =  yend,
                     linetype =  factor(line),
                     colour   =  factor(clust)),
                 lineend      =  "round",
                 show.legend  =  FALSE,
                 size         =  branch.size)
  
  ## orientation
  if (fan) {
    p <- p +
      coord_polar(direction = -1) +
      scale_x_continuous(breaks = NULL,
                         limits = c(0, nrow(ggdendro::label(hcdata)))) +
      scale_y_reverse(breaks = ybreaks)
  } else {
    p <- p + scale_x_continuous(breaks = NULL)
    if (direction %in% c("rl", "lr")) {
      p <- p + coord_flip()
    }
    if (direction %in% c("bt", "lr")) {
      p <- p + scale_y_reverse(breaks = ybreaks)
    } else {
      p <- p + scale_y_continuous(breaks = ybreaks)
      nudge.label <- -(nudge.label)
    }
  }
  
  # labels
  labelParams <- set_labels_params(nrow(hcdata$labels), direction, fan)
  hcdata$labels$angle <- labelParams$angle
  
  p <- p +
    geom_text(data        =  ggdendro::label(hcdata),
              aes(x       =  x,
                  y       =  y,
                  label   =  label,
                  colour  =  factor(clust),
                  angle   =  angle),
              vjust       =  labelParams$vjust,
              hjust       =  labelParams$hjust,
              nudge_y     =  ymax * nudge.label,
              size        =  label.size,
              show.legend =  FALSE)
  
  # colors and limits
  if (!is.null(scale.color)) {
    p <- p + scale_color_manual(values = scale.color)
  }
  
  ylim <- -round(ymax * expand.y, 1)
  p    <- p + expand_limits(y = ylim)
  
  p
}


# UTILS ----
# MACS2 extra columns
extraCols_broadPeak <- c(signalValue="numeric", pValue="numeric", qValue="numeric")
extraCols_narrowPeak <- c(signalValue="numeric", pValue="numeric", qValue="numeric",peak="integer")
`%ni%` <- Negate("%in%")
`%+%` <- function(x, y) {paste0(x, y)}
adf <- as.data.frame
cr_blue <- crayon::make_style("dodgerblue")
cr_green <- crayon::make_style("forestgreen")
cr_fbrick <- crayon::make_style("firebrick")
cr_white <- crayon::make_style("ghostwhite")
cr_orange <- crayon::make_style("darkorange2")

fai_lst <- readRDS("fai_lst.rds")
.xtensions <- c("gff","gff1","gff2","gff3","bed", "narrowPeak", "broadPeak")

fireWarnings <- function(what="",type_warn="..."){
  cat( "\n" %+%
         cr_orange$bold("Warnings") %+%
         cr_white(" : ") %+%
         cr_white$bold(what) %+%
         " " %+%
         cr_blue$bold(type_warn) %+%
         "\n")
}

fireError <- function(what="",type_err="..."){
  stop(cr_white$bold(what) %+%
         " " %+%
         cr_orange$bold(type_err) %+%
         "\n",call. = F)
}

# DATA TRASNFORMATION ----
rangeMinMax <- function(x,min=-1,max=1)
{
  min+((x- min(x,na.rm=T))*(max-min)) /(max(x,na.rm=T)-min(x,na.rm=T))
}

matrix2tibble <- function(m, n_r=41, n_c = 41) {
  colnames(m) <- paste0(seq(1,n_c))
  rownames(m) <- paste0(seq(1,n_r))
  m %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Bait") %>%
    tidyr::pivot_longer(-c(Bait), names_to = "Anchor", values_to = "counts") %>%
    dplyr::mutate(Anchor= forcats::fct_relevel(Anchor,colnames(m))) %>%
    dplyr::mutate(Bait= forcats::fct_relevel(Bait,rev(rownames(m))))
}

redim_matrix <- function(
    mat,
    target_height = 100,
    target_width = 100,
    summary_func = function(x) mean(x, na.rm = TRUE),
    output_type = 0.0, #vapply style
    n_core = 1 # parallel processing
) {
  
  if(target_height > nrow(mat) | target_width > ncol(mat)) {
    stop("Input matrix must be bigger than target width and height.")
  }
  
  seq_height <- round(seq(1, nrow(mat), length.out = target_height + 1))
  seq_width  <- round(seq(1, ncol(mat), length.out = target_width  + 1))
  
  # complicated way to write a double for loop
  do.call(rbind, parallel::mclapply(seq_len(target_height), function(i) { # i is row
    vapply(seq_len(target_width), function(j) { # j is column
      summary_func(
        mat[
          seq(seq_height[i], seq_height[i + 1]),
          seq(seq_width[j] , seq_width[j + 1] )
        ]
      )
    }, output_type)
  }, mc.cores = n_core))
}


# SIGNAL ----
bw_tile <- function(bw, bed, genome="mm10", name = "sig.bw") {
  if(is(bed, "character")){
    anc.gr <- rtracklayer::import.bed(bed)
  } else {
    anc.gr <- bed
  }
  # GenomicRanges::mcols(anc.gr)[name] <- seqsetvis::ssvFetchBigwig(bw, anc.gr,win_size = median(width(anc.gr)))$y
  GenomicRanges::mcols(anc.gr)[name] <- unlist(sevenC::addCovToGR(gr = anc.gr,bwFile= bw, binSize = median(width(anc.gr)))$chip)
  anc.gr
}

bw_to_mcols <- function(bw, bed, genome="mm10", anchor = c("body", "start", "center", "end"), name = "sig.bw", upstream=0, downstream=0){
  
  if(is(bed, "character")){
    anc.gr <- loadranges(bed, genome = genome, seqlevelstyle = "Ensembl",standardcol = T)
  } else {
    anc.gr <- bed
  }
  anc.gr.save <- anc.gr
  
  # keep seqlevelstyle of GR
  # anc.gr.levelStyle <- GenomeInfoDb::seqlevelsStyle(anc.gr)
  anc.gr.levelStyle <- ifelse(all(grepl("chr", (seqnames(anc.gr)))), "UCSC", "Ensembl")
  
  # find seqlevelstyle of bigWig
  bwf <- rtracklayer::BigWigFile(bw)
  bw.levelStyle <- GenomeInfoDb::seqlevelsStyle(bwf)
  
  # Adjust levelstyle of Granges
  if("Ensembl" %in% bw.levelStyle){
    newStyle <- mapSeqlevels(seqlevels(anc.gr), "Ensembl")
    anc.gr <- renameSeqlevels(anc.gr, newStyle[!is.na(newStyle)])
  }
  else{
    newStyle <- mapSeqlevels(seqlevels(anc.gr), "UCSC")
    anc.gr <- renameSeqlevels(anc.gr, newStyle)
  }
  anc.gr <- reframeGR(anc.gr, anchor, upstream, downstream)
  # browser()
  # Compute coverage
  sig.cov <- trackViewer::importData(bw,"BigWig",anc.gr)[[1]]
  # Sometimes, ranges will get out of bound if
  # ranges does not exist in bigwig.
  # Set seqinfo from subset ranges to avoid this behaviour
  common_names <- names(sig.cov)[names(sig.cov) %in% unique(as.vector(seqnames(anc.gr)))]
  gr <- keepSeqlevels(as(sig.cov, "GRanges"), common_names,pruning.mode = "coarse")
  anc.gr <- keepSeqlevels(anc.gr, common_names,pruning.mode = "coarse")
  # seqinfo(gr) <- seqinfo(anc.gr)
  sig.cov <- coverage(gr, weight="score")
  seqinfo(sig.cov) <- seqinfo(anc.gr)
  
  # add mcols to bed Grange
  GenomicRanges::mcols(anc.gr.save)[name] <- rowSums(as.matrix(sig.cov[anc.gr]),na.rm = T)/BiocGenerics::width(anc.gr)
  
  return(anc.gr.save)
}



getCoverageDensity <- function(query, subject, bychr=T) {
  
  missing_chromosomes <- setdiff(unique(seqnames(query)), unique(seqnames(subject)))
  cov_subj <- coverage(subject)
  if (length(missing_chromosomes)) {
    for (chr in missing_chromosomes) {
      cov_subj[[chr]] <- Rle(0)
    } 
  }
  vs_l <- lapply(Views(cov_subj,query), viewSums)
  if(bychr){ # iterate over chr if specified
    s_l <- lapply(vs_l,
                  function(vs) vs/sum(vs))
  }else{
    sumvsl <- sum(unlist(vs_l))
    s_l <- lapply(vs_l,
                  function(vs) vs/sumvsl)
  }
  unlist(s_l)
}

Gooperation <- function(a, b, method){
  if(method == "and")
    res <- a & b
  if(method == "or")
    res <- a | b
  if(method == "add")
    res <- a + b
  if(method == "ratio")
    if(any(b==0))
      res <- (a+1) / (b+1)
  else
    res <- a / b
  if(method == "subtract")
    res <- a - b
  if(method == "mean")
    res <- (a + b)/2
  if(method == "zscore")
    res <- (a - b)/sqrt((a + b)/2)
  if(method == "zscore_norm0"){
    a <- a - min(a, na.rm = T)
    b <- b - min(b, na.rm = T)
    res <- (a - b)/sqrt((a + b)/2)
  }
  if(method == "zscore_norm1"){
    a <- (a - min(a, na.rm = T)) + 1
    b <- (b - min(b, na.rm = T)) + 1
    res <- (a - b)/sqrt((a + b)/2)
  }
  return(res)
}










# PLOT ----
fisher.bool.bool <- function(df, grad_cols, grad_rows, verbose=F, identity = F) {
  mat <- copy(df)
  
  r_df <- NULL
  for(grad_row in grad_rows){
    for(grad_col in grad_cols){
      if (!identity & grad_row==grad_col) next()
      if (sub("(.*?)_(WF_narrowPeak|peaks)","\\1",grad_row) != sub("(.*?)_(WF_narrowPeak|peaks)","\\1",grad_col)) next()
      
      cols.sbset <- c(grad_col, grad_row)
      df <- mat[, ..cols.sbset]
      
      sub_m <- df[, .N, by = c(grad_col, grad_row)]
      sub_m[, pst := apply(.SD, 1, function(x)paste0(x, collapse = "_")), .SDcols = c(names(sub_m)[1], names(sub_m)[2])]
      
      
      if("TRUE_TRUE" %in% sub_m$pst)
        cr <- sub_m[pst=="TRUE_TRUE", N]
      else
        cr <- 0
      
      if("TRUE_FALSE" %in% sub_m$pst)
        cnr <- sub_m[pst=="TRUE_FALSE", N]
      else
        cnr <- 0
      
      if("FALSE_TRUE" %in% sub_m$pst)
        ncr <- sub_m[pst=="FALSE_TRUE", N]
      else
        ncr <- 0
      
      if("FALSE_FALSE" %in% sub_m$pst)
        ncnr <- sub_m[pst=="FALSE_FALSE", N]
      else
        ncnr <- 0
      
      
      f_mat <- matrix(c(L_C=cr,L_nC=ncr,C_nL=cnr,nC_nL=ncnr),nrow=2,ncol=2, byrow = T)
      if(verbose) print(f_mat)
      ft <- fisher.test(f_mat, alternative = "greater")
      mypv <- formatC(as.numeric(ft$p.value, format = "e", digits = 2))
      lfc <- round(log2(ft$estimate),2)
      # deal with inf cases
      lfc <- ifelse(is.finite(lfc),lfc,sign(lfc)*10)
      sign <- ifelse(ft$p.value < 0.05,ifelse(ft$p.value < 0.01,ifelse(ft$p.value < 0.001,"***","**"),"*"),"NS")
      r_df <- rbind.data.frame(r_df,
                               cbind.data.frame(COL=grad_col,ROW=grad_row, intersect = cr,
                                                pv=mypv,lfc=lfc,sign=sign))
    }
  }
  
  
  r_df$ROW <- factor(r_df$ROW, levels = sort(unique(r_df$ROW)))
  r_df$COL <- factor(r_df$COL, levels = sort(unique(r_df$COL)))
  r_df
}

fisher_plot_asym_pdf <- function(DF,main="",norm=NULL,range=c(-2,2)){
  if(!is.null(norm)) DF <- DF %>% mutate(lfc=rangeMinMaxAsym(x = lfc,nmin = -2,nmax=0,pmin = 0,pmax = 2))
  colfunc <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))
  p <- ggplot(DF, aes(x = COL, y = ROW,fill=lfc)) +      geom_tile(col=NA) +
    theme_minimal()+
    coord_equal(ratio=1) +
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(position = "top",expand=c(-1,0)) +
    theme(axis.text.x=element_text(angle = 90,size=15, vjust=0, hjust=0.6,face = "bold",
                                   margin=margin(0,0,0,0)),
          axis.text.y=element_text(size=15, margin=margin(0,0,0,0),face = "bold"),
          axis.title.x = element_blank(),
          legend.title = element_text(size=12, vjust=-2, hjust=0.1,face = "bold"),
          legend.background = element_rect(color="grey",size=0.8, linetype="dashed"),
          axis.title.y = element_blank(),
          legend.title.align=0.1,
          plot.margin = unit(c(0,1,0,1), "cm"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    
    geom_text(aes(label = paste0(sign,"\n",intersect)),size=4) + 
    ggtitle(main) +
    scale_fill_gradientn(colours = colfunc(10) ,
                         guide = guide_colorbar(barwidth = 0.8,
                                                title = 'Log2 \nOddsRatio',
                                                label.theme = element_text(size=10, vjust=1, hjust=.5,face = "bold",angle = 0),
                                                barheight = 10,
                                                nbin = 20,
                                                draw.ulim = FALSE, 
                                                draw.llim = FALSE,
                                                ticks = FALSE),
                         breaks=c(min(range),0,max(range)),
                         labels=c(min(range),0,max(range)),
                         limits=c(min(range),max(range)))
  
  return(p)
}
# RANGES ----
liftover <- function(my.gr,chain,sqLvlstyle="NCBI"){
  seqlevelsStyle(my.gr) = "UCSC"
  my.new.gr <- unlist(liftOver(my.gr,chain))
  seqlevelsStyle(my.new.gr) = sqLvlstyle
  return(my.new.gr)
}

loadranges <- function(con, genome=NULL,seqlevelstyle="NCBI", standardcol = F) {
  
  if(length(con) != 1) fireError("filename is missing or is in wrong format","Please provide valid filename")
  if(!file.exists(con)) fireError("filename does not exist","Please provide valid filename")
  if(tolower(tools::file_ext(con)) %ni% tolower(.xtensions))
    fireError("con file extension not one of",paste(.xtensions,collapse=" , "))
  cur_xt <- tools::file_ext(con)
  if (cur_xt == "bed") {
    if (grepl("Peak_all.bed", con) | grepl("Peak_optimal.bed", con))
      gr <- rtracklayer::import(con, format = "BED", extraCols = extraCols_broadPeak) # chipr like output
    else {
      gr <- data.table::fread(con)
      if(ncol(gr)==3)
        gr <- gr[,.(seqnames = V1, start = V2+1, end = V3)]
      else if(ncol(gr)>=6)
        gr <- gr[,.(seqnames = V1, start = V2+1, end = V3, name = V4, score = V5, strand = V6)]
      gr <- GenomicRanges::GRanges(gr)
    }
    
  } else if (cur_xt == "narrowPeak") {
    gr <-rtracklayer::import(con,
                             format = "BED",
                             extraCols = extraCols_narrowPeak)
  } else if (cur_xt == "broadPeak") {
    gr <- rtracklayer::import(con,
                              format = "BED",
                              extraCols = extraCols_broadPeak)
  } else if (cur_xt == "gff") {
    gr <-rtracklayer::import.gff(con)
  } else if (cur_xt == "gtf1") {
    gr <- rtracklayer::import.gff1(con)
  } else if (cur_xt == "gff2") {
    gr <-rtracklayer::import.gff2(con)
  } else if (cur_xt == "gff3") {
    gr <-rtracklayer::import.gff3(con)
  }
  
  # set names
  if (is.numeric(gr$name) | all(is.na(gr$name)) | is.null(gr$name) | all(gr$name ==".")) {
    names(gr) <- paste(as.vector(GenomicRanges::seqnames(gr)),GenomicRanges::start(gr),GenomicRanges::end(gr), as.vector(GenomicRanges::strand(gr)), seq_along(gr),sep="_")
    GenomicRanges::mcols(gr)[["name"]] <- names(gr)
  } else if (!is.null(gr$name)){
    names(gr) <- gr$name
  }
  
  # Test if genome is loaded
  if(!is.null(genome)){
    gr <- setSeqinfo(gr,genome,seqlevelstyle)
  }else{
    fireWarnings(paste0("genome parameter is NULL"," - Providing a genome parameter can be very useful"))
  }
  
  if (standardcol) {
    cn <- colnames(GenomicRanges::mcols(gr))
    cn <-  cn[!cn %in% c("name")]
    GenomicRanges::mcols(gr)[cn] <- NULL
  }
  
  
  gr
}


reframeGR <- function(gr, anchor, upstream=0, downstream=0){
  is.plus.strand <- BiocGenerics::strand(gr) == "+" | BiocGenerics::strand(gr) == "*"
  
  # resize GRanges (if anchor in TSS/TES/center)
  if(anchor == "start"){
    BiocGenerics::end(gr[is.plus.strand,]) <- BiocGenerics::start(gr[is.plus.strand,])
    BiocGenerics::start(gr[!is.plus.strand,]) <- BiocGenerics::end(gr[!is.plus.strand,])
  }
  if(anchor == "end"){
    BiocGenerics::start(gr[is.plus.strand,]) <- BiocGenerics::end(gr[is.plus.strand,])
    BiocGenerics::end(gr[!is.plus.strand,]) <- BiocGenerics::start(gr[!is.plus.strand,])
  }
  if(anchor == "center"){
    gr <- GenomicRanges::resize(gr, 1, fix = "center")
  }
  
  # extend ranges
  if(upstream != 0 | downstream != 0){
    BiocGenerics::start(gr[is.plus.strand,]) <- BiocGenerics::start(gr[is.plus.strand,])-upstream
    BiocGenerics::end(gr[!is.plus.strand,]) <- BiocGenerics::end(gr[!is.plus.strand,])+upstream
    
    BiocGenerics::end(gr[is.plus.strand,]) <- BiocGenerics::end(gr[is.plus.strand,])+downstream
    BiocGenerics::start(gr[!is.plus.strand,]) <- BiocGenerics::start(gr[!is.plus.strand,])-downstream
    
    GenomicRanges::trim(gr)
  }
  gr
}




go2dt <- function(go) {data.table::data.table(data.frame(seqnames = seqnames(go), strand = strand(go), IRanges::ranges(go,use.mcols = T,use.names = F), check.names = F))}


# SETSEQINFO ----
setSeqinfo <- function(gr,genome,seqlevelstyle="NCBI"){
  seqlevelstyle <- ifelse(seqlevelstyle == "Ensembl","NCBI",seqlevelstyle)
  gnm_sqfo <- fai_lst[[seqlevelstyle]][[genome]]
  
  if (!any(grepl("NCBI",GenomeInfoDb::seqlevelsStyle(gr)))) GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI"
  GenomeInfoDb::seqlevels(gr, pruning.mode="coarse") <- GenomeInfoDb::seqlevels(gnm_sqfo)
  GenomeInfoDb::seqinfo(gr,pruning.mode=c("coarse")) <- gnm_sqfo
  GenomeInfoDb::seqlevels(gr) <-  GenomeInfoDb::seqlevelsInUse(gr)
  gr
}


aggregate_ranges <- function(query=NULL, subject=NULL, fun_agr="mean", fun_opts = "na.rm = T", out_mc=NULL,subject_mc=NULL, fol_opts = NULL){
  if (is.null(query) | !is(query, "GenomicRanges")) fireError("query can't be null and must be a valid GenomicRanges","Please provide a query as GenomicRanges")
  if (is.null(subject)| !is(subject, "GenomicRanges")) fireError("subject can't be null and must be a valid GenomicRanges","Please provide a query as GenomicRanges")
  if (is.null(fun_agr)) fireError("Aggregation funciton can't be NULL","Please provide a functiojn name (mean for instance)")
  if (is.null(out_mc) | out_mc %in% colnames(mcols(query))) fireError("Output mcols for query can't be NULL and can't be already an mcol of query","Please provide a valid outuput mcol name")
  if (is.null(subject_mc) | !(subject_mc %in% colnames(mcols(subject)))) fireError("Subject mcols can't be NULL and must be an existing mcol of subejct","Please provide a valid subject mcol name")
  # browser()
  
  arg_l <- list()
  for (arg in fol_opts) {
    marg <- gsub("[[:space:]]","",strsplit(arg, "=")[[1]][1])
    mval <- gsub("[[:space:]]","",strsplit(arg, "=")[[1]][2])
    if (mval %in% c("T","TRUE","F","FALSE"))mval <-  as.logical(mval)
    else if(!is.na(as.numeric(mval))) mval <- as.numeric(mval)
    arg_l[[marg]] <- mval
  }
  hits <- do.call(
    what = "findOverlaps",
    args = append(
      list(query = query, subject=subject),
      arg_l
    )
  )
  
  # hits <- eval(parse("findOverlaps(query, subject, ...)"))
  hitsByQuery <- as(hits, "List")
  
  lst <- lapply(extractList(mcols(subject)[[subject_mc]], hitsByQuery), unlist)
  mf <- function(type, x, params){
    if (!is.null(params))
      eval(parse(text = paste0(type,"(",quote(x),",",params,")")))
    else
      eval(parse(text = paste0(type,"(",quote(x),")")))
    
  }
  mcols(query)[[out_mc]] <- unlist(lapply(lst, function(x) mf(type = fun_agr, x=x, params=fun_opts)))
  
  query
}



# HIC ----
## APA ----
pairwiseInteractionDT <- function(anc_go=NULL,
                                  bait_go=NULL,
                                  apa_size=41,
                                  d_min=4,
                                  d_max=100,
                                  forceSym = F,
                                  # tile_go=NULL,
                                  cons_nm=NULL,
                                  sym = F,
                                  hic_res=1000){
  chrs <- unique(c(as.vector(seqnames(bait_go)),as.vector(seqnames(anc_go))))
  Go_cons <- NULL
  for (chr in chrs) {
    min_anc <- min(start(subset(anc_go, seqnames == chr)))
    max_anc <- max(end(subset(anc_go, seqnames == chr)))
    min_bait <- min(start(subset(bait_go, seqnames == chr)))
    max_bait <- max(end(subset(bait_go, seqnames == chr)))
    Go_cons <- rbind(Go_cons, data.frame(start = min(min_anc,min_bait), end = max(max_anc,max_bait) ,seqnames = chr))
  }
  
  constraint_go <- GenomicRanges::GRanges(Go_cons)
  names(constraint_go) <- NULL
  #Bait feature
  bait_r_go <- IRanges::resize(bait_go, 1, fix="center") #resize the feature to do the overlap on domains without multi hits
  bait_r_go$cb_b <- ceiling((GenomicRanges::start(bait_r_go)) / hic_res )
  bait_r_go$l_idx_b <- 1:length(bait_r_go)
  mol_dt <- data.table::data.table( data.frame(IRanges::mergeByOverlaps(constraint_go, bait_r_go) ))
  mol_dt$filt <- paste0(mol_dt$constraint_go.seqnames, "_", mol_dt$constraint_go.start)
  mol_b_dt <- mol_dt[ , c( "bait_r_go.name", "chr_bin", "filt", "bait_r_go.seqnames", "bait_r_go.strand" ) ]
  colnames(mol_b_dt) <- c("bait_name", "bait_chrom_bin", "constraint_id", "bait_chrom", "bait_strand")
  data.table::setkey(mol_b_dt, constraint_id)
  #Anchoring feature
  anc_r_go <- IRanges::resize(anc_go, 1, fix="center") #resize the feature to do the overlap on domains without multi hits
  anc_r_go$cb_b <- ceiling((GenomicRanges::start(anc_r_go)) / hic_res )
  anc_r_go$l_idx_b <- 1:length(anc_r_go)
  mol_dt <- data.table::data.table( data.frame(IRanges::mergeByOverlaps(constraint_go, anc_r_go) ))
  mol_dt$filt <- paste0(mol_dt$constraint_go.seqnames, "_", mol_dt$constraint_go.start)
  mol_a_dt <- mol_dt[ , c( "anc_r_go.name", "chr_bin", "filt", "anc_r_go.seqnames", "anc_r_go.strand" ) ]
  colnames(mol_a_dt) <- c("anc_name", "anc_chrom_bin", "constraint_id", "anc_chrom", "anc_strand")
  data.table::setkey(mol_a_dt, constraint_id)
  # Merge
  inter_dt <- mol_b_dt[mol_a_dt, allow.cartesian=TRUE,nomatch=0]
  
  if(nrow(inter_dt) == 0) return(NULL)
  # Mark the one to revert
  inter_dt$rev <- ifelse(inter_dt$anc_chrom_bin < inter_dt$bait_chrom_bin, T, F)
  # Rmove reverse if bait == anchor (duplicated, we don't need it)
  if(sym == "5p3p"){
    inter_dt <- inter_dt[inter_dt$rev==F,]
  } else if (sym == "3p5p") {
    inter_dt <- inter_dt[inter_dt$rev==T,]
  } # else we take both
  # Add distance between int and anc in bin
  inter_dt$dst_bin <- abs(inter_dt$bait_chrom_bin - inter_dt$anc_chrom_bin)
  # Add a key
  inter_dt$key <- 1:nrow(inter_dt)
  # Filter distnace min and max
  inter_dt <- subset(inter_dt, dst_bin >= d_min & dst_bin <= d_max)
  # RETURN DATA.TABLE
  inter_dt
}


dumpFlatInteractionsMatrices <- function (pairwise_dt = NULL,cm = NULL, intra = T, n = 41) {
  if (intra) {
    stack_dt <- NULL
    # We can directly check either bait or anc seqnames as they are the same in intra
    for (chrom in unique(pairwise_dt$bait_chrom)){
      pw_chr_dt <- subset(pairwise_dt, bait_chrom == chrom)
      rv_l <- list()
      for (reverse in unique(pw_chr_dt$rev)){ # For couples that has to be reversed or not
        pw_chr_rev_dt <- subset(pw_chr_dt, rev == reverse)
        rv_l[[paste0("to rev ",reverse)]] <- nrow(pw_chr_rev_dt)
        
        # Subset by chrom
        chr_cm <- cm[[chrom]]
        # Launch fast algorithm
        bb <- pw_chr_rev_dt$bait_chrom_bin
        ba <- pw_chr_rev_dt$anc_chrom_bin
        # Make all bin by repeat
        ba_r <- rep(ba, each=n*n)
        bb_r <- rep(bb,each=n*n)
        # Make mask for extending bins
        mask <- -((n-1)/2):((n-1)/2)
        # apply it to Anch and Bait
        a_mask <- rep(rep(mask, each=n),nrow(pw_chr_rev_dt))
        b_mask <- rep(rep(mask, times=n),nrow(pw_chr_rev_dt))
        # Add to ba & bb
        ba_mask <- ba_r + a_mask
        bb_mask <- bb_r + b_mask
        
        idx <- ceiling(n*n/2)
        j <-  (idx-1)%%n + 1
        i <-  (idx-1)%/%n + 1
        
        bin_i <- i - (n+1)/2 + 20 # The +20 comes from na extend
        bin_j <- j - (n+1)/2 + 20 # The +20 comes from na extend
        itr_vec <- chr_cm[((ba_mask + bin_i)-1) * dim(chr_cm)[1] + bb_mask + bin_j]
        if (reverse) {
          itr_vec <- rev(itr_vec)
          nkey <- rev(pw_chr_rev_dt$key)
        } else {
          nkey <- pw_chr_rev_dt$key
        }
        itr_mat <- matrix(itr_vec, nrow = nrow(pw_chr_rev_dt),ncol=n*n,byrow = T,dimnames = list(paste0("mat",1:nrow(pw_chr_rev_dt)),paste0("pix",1:(n*n))))
        itr_matrices_dt <- data.table::data.table(dtkey=nkey, itr_mat)
        stack_dt <- rbind(stack_dt, itr_matrices_dt)
      }
      gc()
    }
  }
  # # Join tables back to have access to bin info
  data.table::setkey(stack_dt, dtkey)
  data.table::setkey(pairwise_dt, key)
  flat_all_chaos <- pairwise_dt[stack_dt]
  
  flat_all_chaos
}

aggreg_lflat <- function(l_flat_dt){
  
  for(i in names(l_flat_dt)){
    l_flat_dt[[i]] <- rbindlist(l_flat_dt[[i]], idcol = "cm")
  }
  
  flat_dt <- rbindlist(l_flat_dt, idcol = "couple_name")
  flat_dt <- tidyr::separate(flat_dt, "couple_name", into = c("bait","anchor"),sep = "\nVS\n")
  setDT(flat_dt)
}

remove_diag <- function(flat_dt, offset,apa_size){
  dst_rm_offdiag <- unique(subset(flat_dt, dst_bin < (apa_size+offset-1))$dst_bin)
  
  for(i in dst_rm_offdiag){
    # i_rm <-  i - offset
    i_rm <-  abs(i - apa_size - (offset - 1))
    msk <- mask[(apa_size-i_rm):apa_size,1:(i_rm+1)]
    idx <- msk[lower.tri(msk, T)]
    if(length(idx>0))
      flat_dt[dst_bin == i, (paste0("pix", idx)) := NA]
  }
  flat_dt
}


get_apa_metrics_coordinate <- function(apa_size, mask, ctloffset=6, square_size) {
  
  # Check offset consistency
  if(floor(apa_size/2) - (ctloffset + 1) < 0){
    ctloffset <- floor(apa_size/2) - 1
  }
  
  # Coordinate of central pixel
  CP <- (apa_size+1)/2
  CP <- mask[CP, CP]
  
  # Coordinate of 3x3 square centered on central pixel
  m_size <- (apa_size-1)/2 # Middle pixel value - 1
  i_CS <- m_size:(m_size+2)
  j_CS <- m_size:(m_size+2)
  CS <- as.vector(mask[i_CS, j_CS])
  
  # Coordinate of Upper left 3x3 square from 3x3 square at center
  i_ULS <- i_CS-ctloffset
  j_ULS <- j_CS-ctloffset
  ULS <- as.vector(mask[i_ULS, j_ULS])
  
  # Coordinate of Upper Right square starting from central pixel (Compartment)
  i_URS <- 1:floor(apa_size/2)
  j_URS <- (ceiling(apa_size/2)+1):apa_size
  URS <- as.vector(mask[i_URS, j_URS])
  
  # Coordinate of Bottom Right 3x3 square from 3x3 square at center (Compartment)
  i_BRS <- i_CS+ctloffset
  j_BRS <- j_CS+ctloffset
  BRS <- as.vector(mask[i_BRS, j_BRS])
  
  # Coordinate of  Bottom Left square starting from central pixel (TAD)
  i_BLS <- (ceiling(apa_size/2)+1):apa_size
  j_BLS <- 1:floor(apa_size/2)
  BLS <- as.vector(mask[i_BLS, j_BLS])
  
  
  # Outer perimeters squares
  lv2_coord <- floor(apa_size/2) : (floor(apa_size/2) + 2)
  lv3_coord <- (floor(apa_size/2)-1) : (floor(apa_size/2) + 3)
  lv4_coord <- (floor(apa_size/2)-2) : (floor(apa_size/2) + 4)
  dli_9x9 <- (floor(apa_size/2) - 3) : (floor(apa_size/2) + 5)
  
  LV2 <- mask[lv2_coord, lv2_coord][-5]
  LV3 <- mask[lv3_coord, lv3_coord][-c(7:9, 12:14, 17:19)]
  LV4 <- mask[lv4_coord, lv4_coord][-c(9:13, 16:20, 23:27, 30:34,37:41)]
  DLI <- mask[dli_9x9,dli_9x9][-41]
  
  return(list(CP = CP, CS = CS, BLS = BLS, BRS = BRS, ULS = ULS, URS = URS, LV4 = LV4, LV3 = LV3, LV2 = LV2, DLI = DLI))
}


setApaMetrics <- function(apa_dt, coord_l) {
  for (i in names(coord_l)){
    tmp <- apa_dt[, rowMeans(.SD, na.rm = T), .SDcols = paste0("pix",coord_l[[i]])]
    apa_dt[,(i):=tmp]
  }
}

## SADDLE ----
getSaddle2D <- function(bin_go=NULL,#=Go
                        na_offset = 20, # Offset from H5 NA extends to avoid APA batch effect
                        diag_offset = 1, # Offset for removing offdiagonals
                        bait=NULL,#= params.aname, mcols of interest. Name containing score for ranking genomic tiles
                        anchor=NULL,#= params.aname, mcols of interest. Name containing score for ranking genomic tiles
                        bait_quant=30, # 30 digit (size of the saddle) by default
                        anchor_quant=20, # 30 digit (size of the saddle) by default
                        cm_l=NULL, # list of H5 matrix (by chromosomes)
                        f_score="median", # How should be summarized mcols ranker statistics ?
                        f_sad="median", # How should be coomputed each saddle pixel ?
                        f_cm="median", # How should we merge all sub-saddle matrix between chromosomes ?
                        log2=F) {
  
  
  bin_dt <- data.table::data.table(data.frame(bin_go))
  # Compute quantiles by chromosomes
  bin_dt[,quant_bait := dplyr::ntile(.SD,bait_quant), by = "seqnames",.SDcols = bait][,quant_anchor := dplyr::ntile(.SD,anchor_quant), by = "seqnames",.SDcols = anchor]
  # Create summarized value for each quantile class and for each chromosomes
  ranker_df_mergedChr_bait <- bin_dt[,.(tot_bait=get(f_score)(get(bait))), .(quant_bait)]
  ranker_df_mergedChr_anchor <- bin_dt[,.(tot_anchor=get(f_score)(get(anchor))), .(quant_anchor)]
  ranker_df_byChr_bait <- bin_dt[,.(tot_bait=get(f_score)(get(bait))), .(quant_bait,seqnames)]
  ranker_df_byChr_anchor <- bin_dt[,.(tot_anchor=get(f_score)(get(anchor))), .(quant_anchor,seqnames)]
  # Compute saddle pixels
  sadmat_all_cml <- list()
  # get seqnames of the current tile genome
  chr_v <- unique(bin_dt$seqnames)
  
  for(chr in chr_v){
    # Drop the CM for the given chrom and remove off-diagonals -1,0,1
    cm <- cm_l[[chr]]
    # Turn into log2(oe) values and add matrix quantilization
    if(log2) cm <- log2(cm) # Turn obs/exp as linear values
    # cmq <- matrix(ntile(cm,n_tiles),dim(cm)[1],dim(cm)[2])
    # Set given chr saddle matrix
    saddle_m <- matrix(NA,nrow = bait_quant,ncol=anchor_quant)
    # saddle_mq <- matrix(NA,nrow = nb_quant,ncol=nb_quant)
    # Create symetric combinations
    digit_df <- expand.grid(1:bait_quant,1:anchor_quant)
    c_mat <- as.matrix(cm_l[[chr]])
    max_bin <- dim(c_mat)[1] - 20 # Remove NA's extend use in APA
    c_mat <- c_mat[21:max_bin, 21:max_bin] # Remove NA's extend use in APA
    getKDiag(c_mat, (-diag_offset):diag_offset) <- NA 
    c_dt <- data.table(reshape2::melt(c_mat))
    bin_chr_dt <- subset(bin_dt, seqnames == chr)[,c("chr_bin","quant_bait","quant_anchor")]
    setkey(bin_chr_dt, chr_bin)
    setkey(c_dt, Var1)
    c_dt <- c_dt[bin_chr_dt[,.(chr_bin,quant_bait)], nomatch=0]
    setkey(c_dt, Var2)
    c_dt <- c_dt[bin_chr_dt[,.(chr_bin,quant_anchor)], nomatch=0]
    c_dt <- setNames(c_dt, c("bait","anchor","value","bait_quant","anc_quant"))
    sad_mat_chr <- c_dt[,.(tot=get(f_sad)(value,na.rm = T)), .(bait_quant,anc_quant)]
    sadmat_all_cml[[chr]] <- with(sad_mat_chr, tapply(tot, list(bait_quant, anc_quant), FUN = identity))
  }
  
  sadmat_all_cm <- apply(simplify2array(sadmat_all_cml), 1:2, f_cm,na.rm=T)
  sadmat_all_cm[is.na(sadmat_all_cm)] <- 0
  
  return(list(sadmat_cm_aChr=sadmat_all_cml,
              sadmat_cm_mChr=sadmat_all_cm,
              sum_stat_mChr_b=ranker_df_mergedChr_bait,
              sum_stat_mChr_a=ranker_df_mergedChr_anchor,
              sum_stat_aChr_b=ranker_df_byChr_bait,
              sum_stat_aChr_a=ranker_df_byChr_anchor,
              chr_v=chr_v))
}


## MATRICES TRANSFORMATIONS ----
#TODO : GET OR SET K-th diagonal of a matrix
# Input : mat = matrix
#         k = the k-th diagonal (negative for lower diagonals)
# Output : k-th diagonal

getKDiag <- function(mat=NULL,k=0){ # k=0 the first diagonal
  for(i in k){
    print(as.vector(mat[col(mat) - row(mat) == i]))
  }
}

'getKDiag<-' <- function(mat,k,value)
{
  UseMethod('getKDiag<-',mat)
}

'getKDiag<-.matrix' <- function(mat,k,value)
{
  for(i in k){
    mat[col(mat) - row(mat) == i] <- value
  }
  return(mat)
}

## RDS LOADING ----
dropCMbyChr <- function(cml=NULL,chr="",na2zero=T,offset=NULL) {
  if(length(chr)> 1){
    message("Error : Need to specifiy a unique chromosome. E.g chr='chr1'")
    return(F)
  }
  cm <- as.matrix(as.matrix(cml[[chr]]))
  cm[!is.finite(cm)] <- NA # Set 0 to NA
  if(na2zero) cm[cm==0] <- NA # Set 0 to NA
  if(!is.null(offset))
    getKDiag(cm,offset) <- NA # Remove outlier diagonals
  cm
}


## H5 LOADING ----
load_h5_matrix <- function(h5path){
  l_cm <- list()
  # Load H5 matrix
  if (!is.null(h5path)) {
    for(h5.nm in names(h5path)){
      # get h5 name info by h5path
      cm_info <- data.table::data.table(rhdf5::h5ls(h5path[[h5.nm]]))
      if (! nrow(cm_info[otype=="H5I_GROUP","name"])) l_cm[[h5.nm]] <- sapply(cm_info$name, function(nm) { HDF5Array::HDF5Array(h5path[[h5.nm]], nm)})
      else l_cm[[h5.nm]] <- sapply(as.vector(unlist(cm_info[otype=="H5I_GROUP","name"])), function(nm) {print(nm); HDF5Array::TENxMatrix(h5path[[h5.nm]], nm)})
    }
  } else {
    # TODO : put loop here to load cm in l_cm
    mat_nm <- paste0(cond,"_",rep,"_",norm,"_",as.integer(res))
    filt_cm <- filterCmDB(cmname=list(in_=list(value=mat_nm,strict=T)))
    l_cm[["to_dev"]] <- dumpCM(filt_cm)
  }
  l_cm
}


## TADS ----
bm2GR <- function(bm_p, genome=NULL, seqLevelStyle = "Ensembl", dropChr="MtDNA"){
  `%ni%` <- Negate("%in%")
  ts_df <- if(is.null(dropChr))data.table::fread(bm_p,header=F) else data.table::fread(bm_p,header=F) %>% filter(V1%ni%dropChr)
  n_col <- ncol(ts_df)
  x_col <- n_col - 3 # 3 for seqn start end
  ts_v <- paste0("ts",1:x_col)
  colnames(ts_df) <- c("seqnames","start","end",ts_v)
  ts_df$start <- ts_df$start + 1 # avoid overlapping ranges, 0 based here
  ts_gr <- GenomicRanges::GRanges(ts_df,seqinfo = g2i:::fai_lst[[seqLevelStyle]][[genome]])
  if (!is.null(dropChr)) ts_gr <- ts_gr %>% GenomeInfoDb::dropSeqlevels(dropChr,pruning.mode = "coarse")
  ts_gr
}
