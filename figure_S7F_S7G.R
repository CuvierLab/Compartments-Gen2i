# HEADER ====================================================== =
# GEN2I
# Supplementary Figure S7F and S7G
# Repeat-element expression in the compartment-impaired mutants:
#   S7G, family-level log2 fold change versus N2 (heatmap)
#   S7F, the same fold changes for Pao, MULE-MuDR and CMC-Mirage, plus the
#        change of the whole repeat compartment ("repeat total")
# Both panels come from the same analysis, hence a single script.
# ============================================================= =
#
# The mapping and the counting of the repeat copies are NOT redone here: this
# script starts from the count tables produced upstream (featureCounts
# -M --fraction -p -F SAF -s 2 over the RepeatMasker copies of ce11, and the
# DESeq2 size factors estimated on the gene-level counts of the same bam files).
#
#   data/RNAseq/raw_repeat_counts.txt   raw counts, one row per repeat copy
#   data/RNAseq/library_totals.txt      featureCounts library size per sample
#   data/RNAseq/size_factors.txt        gene-level DESeq2 size factors
#   data/RNAseq/rmsk_ce11_annotations.txt  repClass / repFamily of each copy
#   data/RNAseq/rmsk_ce11_genic.txt     whether the copy overlaps an annotated gene
#
# Two quantities are shown, because they answer two different questions:
#   CELLS (S7G) log2FC per family versus N2  -> "which families stand out?"
#   S7F         the same for three families, plus the summed change of the
#               whole compartment -> "does the compartment move as a block?"
# ============================================================= =


# WORKING DIRECTORY ----
# Run the scripts from the root of this repository (where manuscript_lib.R sits):
# setwd("/path/to/Compartments-Gen2i")

# SOURCE ----
source("manuscript_lib.R")
suppressPackageStartupMessages(library(DESeq2))

# PARAMETERS ----
out_d <- "output/figure_S7F_S7G"
dir.create(out_d, showWarnings = FALSE, recursive = TRUE)

raw_p  <- "data/RNAseq/raw_repeat_counts.txt"
lib_p  <- "data/RNAseq/library_totals.txt"
sf_p   <- "data/RNAseq/size_factors.txt"
anno_p <- "data/RNAseq/rmsk_ce11_annotations.txt"
genic_p<- "data/RNAseq/rmsk_ce11_genic.txt"

ref_cond  <- "N2"
grp_col   <- "repFamily"
class_col <- "repClass"

# The four genotypes carried by the manuscript. The allelic series (hpl2a,
# hpl2bc) and I158A are not shown in any figure.
mut_v <- c("hpl2","lin61","hpl2.lin61","met2.set25.set32")
cond_lab <- c("hpl2"             = "hpl-2",
              "lin61"            = "lin-61",
              "hpl2.lin61"       = "hpl-2; lin-61",
              "met2.set25.set32" = "met-2 set-25 set-32")

# Purely technical classes. rRNA is removed BEFORE any normalization: 4 copies
# out of 39,244 carry 70 % of the repeat-assigned reads and their share ranges
# from 17 % to 85 % BETWEEN REPLICATES of the same genotype.
drop_class <- c("rRNA","ARTEFACT")
# Tentative RepeatMasker assignments (the "?" suffix): 0.35 % of the reads.
drop_tentative <- TRUE
# A family is drawn if it carries at least this many INTERGENIC copies. A
# structural threshold, never an expression one, which would condition on the
# denominator of the ratio.
min_copies <- 20L
# Families on which no claim is made whatever the result: the C. elegans SINE
# (tRNA-RTE) is a Pol III, non-polyadenylated transcript, under-sampled 4- to
# 7-fold by poly(A) selection. They keep their tile but lose their stars and
# carry a dagger.
excl_grp <- "tRNA-RTE"
# Bounds of the colour scale; NULL means the range of the data.
fc_lim <- NULL
# The three families reported in S7F.
f_fam <- c("Pao", "MULE-MuDR", "CMC-Mirage")

# SETSEED ----
set.seed(123)

# RUN ----
## 1. ANNOTATION, RAW COUNTS AND SCOPE ----
ann_dt <- merge(fread(anno_p, select = c("name", class_col, grp_col)),
                fread(genic_p), by = "name")
setnames(ann_dt, c(grp_col, class_col), c("grp", "class"))
ann_dt[, genic_ovl := as.logical(genic_ovl)]
ann_dt <- unique(ann_dt[, .(name, grp, class, genic_ovl)])

raw_dt <- fread(raw_p)
lib_dt <- fread(lib_p)

# Intergenic copies only: 61 % of the copies overlap a gene and carry ~80 % of
# the signal. A copy sitting in an intron collects the poly(A) reads of its host.
keep_dt <- ann_dt[!class %in% drop_class & !genic_ovl]
if (drop_tentative)
  keep_dt <- keep_dt[!grepl("\\?$", grp) & !grepl("\\?$", class)]

cat(sprintf("  %d copies kept out of %d (intergenic, without %s%s)\n",
            nrow(keep_dt), nrow(ann_dt), paste(drop_class, collapse = "/"),
            if (drop_tentative) ", without '?' assignments" else ""))

smp_all <- setdiff(colnames(raw_dt), "name")
cond_all <- unique(sub("\\.rep[0-9]+$", "", smp_all))
if (!ref_cond %in% cond_all)
  stop("Reference condition missing from the counts: ", ref_cond)
if (!all(mut_v %in% cond_all))
  stop("Conditions missing from the counts: ",
       paste(setdiff(mut_v, cond_all), collapse = ", "))

smp_v <- smp_all[sub("\\.rep[0-9]+$", "", smp_all) %in% c(ref_cond, mut_v)]
cnt_dt <- merge(raw_dt[, c("name", smp_v), with = FALSE], keep_dt, by = "name")

## 2. IMPOSED GENE-LEVEL SIZE FACTORS ----
# Estimating them on the repeat compartment itself would normalize away the very
# effect being measured; the featureCounts library size counts everything that
# was processed, unassigned and multi-mapping reads included. The gene-level
# DESeq2 size factors of the same bam files are independent of both.
sf_dt <- fread(sf_p)
sf_v  <- setNames(sf_dt$size_factor, sf_dt$sample)[smp_v]
if (anyNA(sf_v)) stop("Missing size factors for: ",
                      paste(smp_v[is.na(sf_v)], collapse = ", "))

copy_m <- round(as.matrix(cnt_dt[, ..smp_v]))
rownames(copy_m) <- cnt_dt$name

## 3. DESeq2 AT THE FAMILY LEVEL ----
fam_dt <- cnt_dt[, lapply(.SD, function(x) sum(as.numeric(x))), by = grp, .SDcols = smp_v]
n_copy <- cnt_dt[, .(n_copies = .N), by = grp]
fam_dt <- merge(fam_dt, n_copy, by = "grp")[n_copies >= min_copies]
if (nrow(fam_dt) < 2)
  stop("Fewer than two families above ", min_copies, " intergenic copies")
cat(sprintf("  %d families with >= %d intergenic copies\n", nrow(fam_dt), min_copies))

fam_m <- round(as.matrix(fam_dt[, ..smp_v]))
rownames(fam_m) <- fam_dt$grp

col_df <- data.frame(
  condition = factor(sub("\\.rep[0-9]+$", "", smp_v), levels = c(ref_cond, mut_v)),
  row.names = smp_v)

dds <- DESeq2::DESeqDataSetFromMatrix(fam_m, col_df, ~ condition)
DESeq2::sizeFactors(dds) <- sf_v[smp_v]

# On ~20 features the parametric dispersion fit often fails. Degrade cleanly,
# but say which fit was used: it belongs in the Methods.
fit_used <- "parametric"
dds <- tryCatch(DESeq2::DESeq(dds, quiet = TRUE), error = function(e) {
  fit_used <<- "local"
  tryCatch(DESeq2::DESeq(dds, fitType = "local", quiet = TRUE), error = function(e2) {
    fit_used <<- "mean"
    DESeq2::DESeq(dds, fitType = "mean", quiet = TRUE)
  })
})
cat("  dispersion fitType : ", fit_used, "\n", sep = "")

de_dt <- rbindlist(lapply(mut_v, function(mu) {
  r <- as.data.frame(DESeq2::results(dds, contrast = c("condition", mu, ref_cond)))
  data.table(grp = rownames(r), condition = mu,
             log2FC = r$log2FoldChange, lfcSE = r$lfcSE,
             pval = r$pvalue, padj = r$padj, baseMean = r$baseMean)
}))

## 4. REPLICATE CONCORDANCE ----
# At n = 2 this is the deciding criterion: a family whose two replicates
# disagree must not be presented as a result, whatever the FDR.
norm_m <- sweep(fam_m, 2, sf_v[smp_v], "/")
rep_v  <- unique(sub("^.*\\.(rep[0-9]+)$", "\\1", smp_v))

conc_dt <- rbindlist(lapply(mut_v, function(mu) {
  rbindlist(lapply(rep_v, function(rp) {
    a <- paste0(mu, ".", rp); b <- paste0(ref_cond, ".", rp)
    if (!all(c(a, b) %in% colnames(norm_m))) return(NULL)
    data.table(grp = rownames(norm_m), condition = mu, rep = rp,
               lr = log2((norm_m[, a] + 1) / (norm_m[, b] + 1)))
  }))
}))
if (nrow(conc_dt) && length(rep_v) >= 2) {
  conc_dt <- dcast(conc_dt, grp + condition ~ rep, value.var = "lr")
  conc_dt[, same_sign := sign(get(rep_v[1])) == sign(get(rep_v[2]))]
  conc_dt[, spread    := abs(get(rep_v[1]) - get(rep_v[2]))]
  de_dt <- merge(de_dt, conc_dt[, .(grp, condition, same_sign, spread)],
                 by = c("grp", "condition"), all.x = TRUE)
} else {
  de_dt[, `:=`(same_sign = NA, spread = NA_real_)]
}

de_dt <- merge(de_dt, unique(keep_dt[, .(grp, class)]), by = "grp")
de_dt <- merge(de_dt, n_copy, by = "grp")

## 5. WEIGHT OF EACH FAMILY IN THE COMPARTMENT ----
# The cells of step 3 are normalized on the median copy, so they say nothing
# about the mass that moved. The weight is expressed in PERCENT OF THE N2 REPEAT
# COMPARTMENT, on counts normalized by the same gene-level size factors, so the
# two quantities are on the same scale. The sum of a column is the change of the
# whole compartment: that is the "repeat total" of S7F.
fam_raw <- cnt_dt[, lapply(.SD, function(x) sum(as.numeric(x))), by = grp, .SDcols = smp_v]
fam_raw <- fam_raw[grp %in% fam_dt$grp]

nrm_m <- sweep(as.matrix(fam_raw[, ..smp_v]), 2, sf_v, "/")
rownames(nrm_m) <- fam_raw$grp

cond_of <- sub("\\.rep[0-9]+$", "", smp_v)
nrm_cond <- vapply(c(ref_cond, mut_v), function(cd)
  rowMeans(nrm_m[, cond_of == cd, drop = FALSE]), numeric(nrow(nrm_m)))

ref_tot <- sum(nrm_cond[, ref_cond])
ctb_dt <- rbindlist(lapply(mut_v, function(mu)
  data.table(grp = rownames(nrm_cond), condition = mu,
             share_ref  = 100 * nrm_cond[, ref_cond] / ref_tot,
             share_mut  = 100 * nrm_cond[, mu]       / ref_tot,
             contrib_pt = 100 * (nrm_cond[, mu] - nrm_cond[, ref_cond]) / ref_tot)))

tot_dt <- ctb_dt[, .(total_pt = sum(contrib_pt)), by = condition]
tot_dt[, fold := 1 + total_pt / 100]

# One descriptive landmark, and only one, comes from the library size: the share
# of the transcriptome taken by the repeats IN N2. No between-genotype
# comparison depends on it.
lib_v <- setNames(lib_dt$lib_total, lib_dt$sample)[smp_v]
ref_smp <- smp_v[cond_of == ref_cond]
ref_pct_lib <- 100 * mean(colSums(copy_m[, ref_smp, drop = FALSE]) / lib_v[ref_smp])
tot_dt[, ref_pct_of_library := ref_pct_lib]
cat(sprintf("  repeat compartment in %s: %.4f %% of library\n", ref_cond, ref_pct_lib))
print(tot_dt)

de_dt <- merge(de_dt, ctb_dt[, .(grp, condition, share_ref, share_mut, contrib_pt)],
               by = c("grp", "condition"), all.x = TRUE)

## 6. FIGURES ----
theme_g2i <- theme_minimal(base_size = 10) +
  theme(plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid   = element_blank(),
        axis.text.x  = element_text(angle = 30, hjust = 1, face = "bold", size = 9),
        axis.text.y  = element_text(size = 8),
        title        = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 7.5, colour = "#555555", face = "plain"))

star_of <- function(p) fifelse(is.na(p), "",
                       fifelse(p < 0.001, "***",
                       fifelse(p < 0.01,  "**",
                       fifelse(p < 0.05,  "*", ""))))

lab_of <- function(x) if (is.null(cond_lab)) x else
  ifelse(x %in% names(cond_lab), cond_lab[x], x)

pl_dt <- copy(de_dt)
# No star on families that cannot be measured.
pl_dt[, star := fifelse(grp %in% excl_grp, "", star_of(padj))]
pl_dt[, condition := factor(lab_of(condition), levels = lab_of(mut_v))]

# Row order: by class, then by decreasing mean effect.
ord <- pl_dt[, .(m = mean(log2FC, na.rm = TRUE)), by = .(class, grp)]
setorder(ord, class, -m)
lv    <- as.character(ord$grp)
lv_sh <- ifelse(lv %in% excl_grp, paste0(lv, " ‡"), lv)
relab <- function(x) factor(ifelse(x %in% excl_grp, paste0(x, " ‡"), x),
                            levels = rev(lv_sh))
pl_dt[, grp   := relab(grp)]
pl_dt[, class := factor(class, levels = sort(unique(as.character(ord$class))))]

lim <- if (!is.null(fc_lim)) fc_lim else max(abs(pl_dt$log2FC), na.rm = TRUE)

sub_txt <- paste0("Intergenic copies, rRNA excluded | * padj<0.05, ** <0.01, *** <0.001",
                  " | x = replicates disagree in sign | n = 2 per genotype",
                  " | dispersion fit: ", fit_used,
                  if (length(excl_grp)) paste0("\n‡ no claim made: Pol III / non-polyadenylated,",
                                               " under-sampled 4-7x by poly(A) selection") else "")

gg <- list()

### S7G : which family moves relative to its peers ----
gg[["Matrix.Cells"]] <-
  ggplot(pl_dt, aes(x = condition, y = grp, fill = log2FC)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  # A cross marks the families whose two replicates disagree in sign. Nothing is
  # removed: the reader sees the family AND its lack of reproducibility.
  geom_point(data = pl_dt[same_sign %in% FALSE], shape = 4, size = 2.6,
             colour = "grey20", stroke = 0.8, show.legend = FALSE) +
  geom_text(aes(label = star), size = 3.6, vjust = 0.75, colour = "grey10") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-lim, lim), oob = scales::squish,
                       name = "log2 FC\nvs N2") +
  facet_grid(class ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(title = "Change: family expression vs N2",
       x = NULL, y = NULL) +
  theme_g2i +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 8))

### S7F : the three families, and the whole compartment ----
f_dt <- de_dt[grp %in% f_fam]
f_dt[, fold := 2^log2FC]
f_dt[, sig  := fifelse(is.na(padj), "na",
               fifelse(padj < 0.001, "***",
               fifelse(padj < 0.01,  "**",
               fifelse(padj < 0.05,  "*", "ns"))))]
f_dt[, cell := sprintf("%.2fx %s", fold, sig)]
f_wide <- dcast(f_dt, condition ~ grp, value.var = "cell")
f_wide <- merge(tot_dt[, .(condition, `Repeat total` = sprintf("%.2fx", fold))],
                f_wide, by = "condition")
f_wide[, condition := lab_of(condition)]
setcolorder(f_wide, c("condition", "Repeat total", "Pao", "MULE-MuDR", "CMC-Mirage"))
f_wide <- f_wide[match(lab_of(mut_v), condition)]

tab_dt <- melt(f_wide, id.vars = "condition", variable.name = "col",
               value.name = "txt", variable.factor = FALSE)
tab_dt[, col := factor(col, levels = c("Repeat total", "Pao", "MULE-MuDR", "CMC-Mirage"))]
tab_dt[, condition := factor(condition, levels = rev(lab_of(mut_v)))]

gg[["Table.Summary"]] <-
  ggplot(tab_dt, aes(x = col, y = condition)) +
  geom_text(aes(label = txt), size = 3.4) +
  scale_x_discrete(position = "top",
                   labels = c("Repeat total" = "Repeat\ntotal",
                              "Pao" = "Pao\n(LTR)",
                              "MULE-MuDR" = "MULE-MuDR\n(DNA)",
                              "CMC-Mirage" = "CMC-Mirage\n(DNA)")) +
  labs(title = "Fold change of repeat-derived reads, mutant vs N2",
       subtitle = paste0("Repeat total: ratio of the summed normalized counts over the ",
                         nrow(fam_dt), " families of S7G (descriptive, not tested)",
                         "\nFamilies: DESeq2 Wald test, BH-corrected over the ",
                         nrow(fam_dt), " families of each contrast"),
       x = NULL, y = NULL) +
  theme_g2i +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 9),
        axis.text.y = element_text(size = 9, face = "italic"),
        axis.ticks = element_blank())

## 7. TABLES ----
setorder(de_dt, condition, -log2FC)
setnames(de_dt, c("grp", "class"), c(grp_col, class_col))
fwrite(de_dt,   file.path(out_d, paste0("Matrix_", grp_col, "_values.csv")))
fwrite(tot_dt,  file.path(out_d, "Compartment_total.csv"))
fwrite(f_wide,  file.path(out_d, "Table_S7F.csv"))
fwrite(data.table(sample = names(sf_v), size_factor = as.numeric(sf_v)),
       file.path(out_d, "SizeFactors_compartment.csv"))

# SAVE ----
ggsave(file.path(out_d, "Matrix.Cells.png"),   gg[["Matrix.Cells"]],   width = 6, height = 8, dpi = 150)
ggsave(file.path(out_d, "Table.Summary.png"),  gg[["Table.Summary"]],  width = 7, height = 3, dpi = 150)
