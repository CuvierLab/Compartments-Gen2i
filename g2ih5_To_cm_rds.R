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

#  READ I/O ----
na_offset <- 20
o_d <- "data/HiC"

# RUN ----
res  <- "25kb"
gp_all_xp <- c('N2-old', 'I158A-old', 'hpl2-old', 'hpl2-lin61-old', 'lin61-old', 'met2-set25-set32-old', 'CEC4')
norml = c('norm.KR','obs_exp')

for(xp in gp_all_xp) {
  
  for (norm in norml) {
    
    h5_path <- list(mat=paste0("data/HiC/",xp,"_",res,".",norm,".g2i.h5"))
    
    h5_mh5_l <- load_h5_matrix(h5_path)[["mat"]]
    
    cm_l <- NULL
    for (chr in names(h5_mh5_l)) {
      cm <- as.matrix(as.matrix(h5_mh5_l[[chr]]))
      # REMOVE NA_OFFSET IF ANY
      max_bin <- dim(cm)[1] - na_offset # Remove NA's extend use in APA
      cm_l[[chr]] <- cm[(na_offset+1):max_bin, (na_offset+1):max_bin]
    }
    saveRDS(cm_l,file = paste0(o_d,'/',xp,"_",res,".",norm,".cm.rds"))
  }
}
# ============================================================= =