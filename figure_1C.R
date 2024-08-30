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
## CONTACT MATRICES HEATMAP ----
gg <- NULL

for (condition in c("CEC4","hpl2","hpl2-lin61","I158A","lin61","met2-set25-set32" )) {
  
  lmat = list()
  lmat[[condition]] = paste0('data/HiC/',condition,'_merged.bwa_mem.25kb.obs_exp.g2i.h5')
  params = list(
    upperMat = list(N2 = 'data/HiC/N2_merged.bwa_mem.25kb.obs_exp.g2i.h5'),
    lowerMat = lmat,
    na_offset = 20,
    diag_offset = 1,
    palette = 'turbo',
    quant = 50,
    legend = TRUE
  )
  
  upperMat <- params[["upperMat"]]
  lowerMat <- params[["lowerMat"]]
  vMin <- ifelse(is.null(params[["vMin"]]), -Inf, params[["vMin"]])
  vMax <- ifelse(is.null(params[["vMax"]]), Inf, params[["vMax"]])
  na_offset <- ifelse(is.null(params[["na_offset"]]), 0, params[["na_offset"]])
  diag_offset <- ifelse(is.null(params[["diag_offset"]]), 1, params[["diag_offset"]])
  ## GRAPHICAL PARAMETERS
  palette <- params[["palette"]]
  quant <- params[["quant"]]
  transform <- params[["transform"]]
  redim_mat <- params[["redim_mat"]]
  
  # Create dt for better handling of each comparisons
  a <- tibble::enframe(params[["upperMat"]]) %>% tidyr::unnest(cols = value) %>% data.table::data.table()
  b <- tibble::enframe(params[["lowerMat"]]) %>% tidyr::unnest(cols = value) %>% data.table::data.table()
  dt <- cbind(a,b)
  data.table::setnames(dt, c("nameUp", "pathUp", "nameLow","pathLow"))
  
  for (i in 1:nrow(dt)) {
    nm_u <- as.character(dt[i, "nameUp"])
    nm_l <- as.character(dt[i, "nameLow"])
    
    h5u_l <- list(as.vector(unlist(dt[i, "pathUp"])))
    h5l_l <- list(as.vector(unlist(dt[i, "pathLow"])))
    
    names(h5u_l) <- nm_u
    names(h5l_l) <- nm_l
    
    cm_u_lbc <- load_h5_matrix(h5path = h5u_l) # lbc = list by chromosomes
    cm_l_lbc <- load_h5_matrix(h5path = h5l_l)
    
    chr_v <- intersect(names(cm_u_lbc[[nm_u]]), names(cm_l_lbc[[nm_l]]))
    
    chr_v <- "I"
    for(chr in chr_v){
      
      cm_u <- as.matrix(as.matrix(cm_u_lbc[[nm_u]][[chr]]))
      cm_l <- as.matrix(as.matrix(cm_l_lbc[[nm_l]][[chr]]))
      
      max_bin <- dim(cm_u)[1] - na_offset # Remove NA's extend use in APA
      cm_u <- cm_u[(na_offset+1):max_bin, (na_offset+1):max_bin] # Remove NA's extend use in APA
      cm_l <- cm_l[(na_offset+1):max_bin, (na_offset+1):max_bin] # Remove NA's extend use in APA
      
      
      
      cm_u[cm_u == 0] <- NA
      cm_l[cm_l == 0] <- NA
      
      getKDiag(cm_u, -diag_offset:diag_offset) <- NA
      getKDiag(cm_l, -diag_offset:diag_offset) <- NA
      
      cm_u[lower.tri(cm_u)] <- 0
      cm_l[upper.tri(cm_l)] <- 0
      getKDiag(cm_l, 0) <- 0
      
      if (!is.null(quant)) {
        cmq_u <- matrix(dplyr::ntile(as.vector(cm_u), quant), nrow = nrow(cm_u), ncol = ncol(cm_u))
        cmq_l <- matrix(dplyr::ntile(as.vector(cm_l), quant), nrow = nrow(cm_l), ncol = ncol(cm_l))
        # remove lower and upper accordingly
        cmq_u[lower.tri(cmq_u)] <- 0
        cmq_l[upper.tri(cmq_l)] <- 0
        cmq <- if(!is.null(redim_mat)) redim_matrix(cmq_u+cmq_l, target_height=redim_mat, target_width=redim_mat, n_core=12) else cmq_u+cmq_l
        
        cmq <- melt(cmq)
        colnames(cmq) <- c("x","y","value")
        
        cmqp <- ggplot(cmq, aes(x = x, y = y, fill = value)) +
          geom_raster() +
          scale_fill_gradientn(colors = viridis::viridis_pal(option =palette)(quant)) +
          theme_void() + ggtitle(paste0("Upper : ",nm_u,"\n", "Lower : ",nm_l,"\n","chr: ", chr))
        gg[[paste0("ContactMatrices.Quantilized.",nm_u,"_vs_",nm_l,".Chromosome_", chr)]] <- cmqp
      }
    }
  }
}
