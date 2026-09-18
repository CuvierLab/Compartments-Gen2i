# HEADER ====================================================== =
# GEN2I
# NAR revision 2026
# R version : R version 4.4.0 (2024-04-24), tested on 4.5.3
# System : x86_64, linux-gnu
# ============================================================= =
#
# Figures 2B and 2C (PC1), Supplementary Figures S4B and S4C (PC2) - B to A
# compartment metaprofiles and the amplitude of the transition.
#
# 2B / S4B : metaprofile of the eigenvector around B -> A compartment
#            transitions, one panel per mutant with wild type recalled in grey,
#            and a bootstrap confidence band computed over the transitions.
# 2C / S4C : amplitude of the transition (A side minus B side), one value per
#            transition, shown genotype by genotype.
#
# Unit of analysis: the transition, not the bin. The bins of a metaprofile are
# consecutive points of a single averaged curve, strongly autocorrelated, and
# their number has nothing to do with the number of observations; a test run
# over them would measure autocorrelation rather than evidence. Every statistic
# below is therefore computed over the individual B -> A transitions.
#
# Bootstrap band. Transitions are resampled with replacement and the mean
# profile is recomputed at each draw, giving a percentile interval at every
# position. The same draw is applied to every genotype at each iteration, so
# the pairing is preserved and the band reflects locus-to-locus variability
# rather than an artificial reshuffling between genotypes.
#
# Baseline alignment. Each profile is brought to zero on its own B side
# (-zone_size to 0), transition by transition and genotype by genotype, so that
# all curves share a baseline and the step can be read directly. The anchor is
# the B side and NOT the far left of the window: at -300 to -175 kb a good share
# of the transitions have already left the B compartment, so that "baseline" is
# noise and subtracting it widens the bands instead of narrowing them (on
# PC1 / 100 kb, the N2 band goes from 0.017 to 0.031). Anchoring on the B side
# also has the merit that the height of the aligned profile on the right IS the
# amplitude quantified in the next panel. Alignment shifts each profile by a
# constant and therefore leaves the amplitude untouched, the amplitude being a
# difference internal to each profile.
#
# Amplitude. The claim is about the SHAPE of the profile, not its level: single
# mutants keep a wild-type-like B -> A step, double and triple mutants are
# flattened. A profile merely shifted lower keeps its amplitude; a flat profile
# has zero amplitude whatever its height. Two tests accompany it in
# Amplitude_Tests.csv: against 0 ("is there still a transition?", the paper's
# actual claim, stated positively for each genotype) and against wild type
# ("is the remaining transition smaller?", where a mutant that is NOT
# significant resembles the wild type). The upper bound of the confidence
# interval is also reported as a percentage of wild type (CI_high_pct_WT): it is
# what tells the reader what is excluded when no transition comes out.
#
# Caveat on min_size. It filters on compartment length and the distribution is
# heavily skewed (median B run: 2 bins, i.e. 50 kb), so it is a very aggressive
# filter: on chromosomes I-III, 200 kb leaves a single transition, against 10 at
# 100 kb and 27 at 50 kb. The published panels use 100 kb. The script prints the
# number of transitions and refuses to compute below five.
#
# Sign convention inherited from the original metaprofile: the eigenvector of
# met2-set25-set32 comes out with the opposite sign and is flipped so that all
# profiles are oriented the same way. The sign of an eigenvector is arbitrary,
# so this changes no structure.
# ============================================================= =

# WORKING DIRECTORY ----
# Run the scripts from the root of this repository (where manuscript_lib.R sits):
# setwd("/path/to/Compartments-Gen2i")

# SOURCE ----
source("manuscript_lib.R")

# PARAMETERS ----
gp_xp      <- c("CEC4", "hpl2-lin61-old", "hpl2-old", "lin61-old",
                "met2-set25-set32-old", "N2-old")            # initial batch, I158A excluded
ref_xp     <- "N2-old"                                       # wild type, defines the transitions
chroms     <- c("I", "II", "III")
bin_size   <- 25000
min_size   <- 100000                                         # minimum length of each compartment
flank_size <- 300000                                         # half-width of the plotted window
zone_size  <- flank_size / 2                                 # half-window defining the two sides
n_boot     <- 2000L
conf_level <- 0.95
boot_seed  <- 123
jitter_seed <- 123                                           # cosmetic only, see note below

## one output directory and one panel name per eigenvector
panels <- list(
  pca1 = list(out_d = "output/figure_2B_2C",
              profile = "figure_2B", amplitude = "figure_2C"),
  pca2 = list(out_d = "output/figure_S4B_S4C",
              profile = "figure_S4B", amplitude = "figure_S4C"))

min_bins   <- ceiling(min_size / bin_size)
flank_bins <- ceiling(flank_size / bin_size)
alpha      <- (1 - conf_level) / 2
star <- function(p) ifelse(is.na(p), "-",
                    ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                    ifelse(p < 0.05, "*", "ns"))))

# RUN ----
for (pca in names(panels)) {

  out_d <- panels[[pca]][["out_d"]]
  dir.create(out_d, showWarnings = FALSE, recursive = TRUE)

  ## ---- GENOMIC OBJECT: bins of 25 kb, one eigenvector column per genotype ---
  ## The serialized column names of the analysis pipeline are reproduced here
  ## ("eigen_pca1_N2.old"), because the genotype label of every table and panel
  ## is read back from them.
  eigen_col <- setNames(paste0("eigen_", pca, "_", gsub("-", ".", gp_xp)), gp_xp)
  ref_col   <- eigen_col[[ref_xp]]

  Go <- loadranges("data/ChIPseq/ce11_tiled_25kb.bed", genome = "ce11")
  for (xp in gp_xp)
    Go <- bw_tile(bw     = sprintf("data/HiC/%s_merged.bwa_mem.%dkb.%s.bw",
                                   xp, bin_size / 1e3, pca),
                  bed    = Go,
                  genome = "ce11",
                  name   = eigen_col[[xp]])

  go_dt <- go2dt(Go)

  ## ---- 1. ONE PROFILE PER B -> A TRANSITION --------------------------------
  cat(sprintf("Step 1 (%s): extracting one profile per B->A transition...\n", pca))

  tr_profiles <- list()
  tr_chr      <- list()

  for (xp in gp_xp) {

    eigencol  <- eigen_col[[xp]]
    condition <- sub("eigen_pca.+_(.+)", "\\1", eigencol)
    torev     <- grepl("met2", condition)                    # see sign convention above

    for (chr in chroms) {

      chr_data <- go_dt[seqnames == chr]

      ## Transitions are defined ONCE, on the reference: every genotype is read
      ## at the same positions, which is what makes the pairing valid.
      rl <- rle(as.vector(ifelse(chr_data[, get(ref_col)] > 0, "A", "B")))
      L  <- rl$lengths
      V  <- rl$values
      ok <- L >= min_bins
      cs <- cumsum(L * bin_size)

      fwd  <- which(V[-length(V)] == "B" & ok[-length(ok)] & V[-1] == "A" & ok[-1])
      ## A -> B transitions are reversed below to become B -> A
      rev_ <- which(V[-length(V)] == "A" & ok[-length(ok)] & V[-1] == "B" & ok[-1])

      pl <- list()
      for (sh in cs[fwd]) {
        wi <- which(chr_data$start >= sh - flank_size & chr_data$end <= sh + flank_size)
        if (length(wi) == 2 * flank_bins)
          pl[[length(pl) + 1]] <- ifelse(torev, -1, 1) * chr_data[wi, get(eigencol)]
      }
      for (sh in cs[rev_]) {
        wi <- which(chr_data$start >= sh - flank_size & chr_data$end <= sh + flank_size)
        if (length(wi) == 2 * flank_bins)
          pl[[length(pl) + 1]] <- ifelse(torev, -1, 1) * rev(chr_data[wi, get(eigencol)])
      }

      if (length(pl)) {
        m <- do.call(cbind, pl)
        tr_profiles[[condition]] <- cbind(tr_profiles[[condition]], m)
        tr_chr[[condition]]      <- c(tr_chr[[condition]], rep(chr, ncol(m)))
      }
    }
  }

  ## Alphabetical order, and not the order of gp_xp: it is the one ggplot
  ## applies implicitly, hence the only one that guarantees that a genotype
  ## keeps the SAME Dark2 colour from one panel to the next.
  conds    <- sort(names(tr_profiles))
  ref_cond <- conds[conds %like% "N2"][1]
  pos      <- setdiff(seq(-flank_size, flank_size, by = bin_size), 0)
  n_tr     <- unique(vapply(tr_profiles, ncol, integer(1)))

  if (length(n_tr) != 1)
    stop("Number of transitions differs between conditions: pairing impossible.")
  if (n_tr < 5 || nrow(tr_profiles[[ref_cond]]) != length(pos))
    stop(sprintf(paste0("Analysis not possible: %d transition(s) retained for ",
                        "compartments >= %g kb. The metaprofile would be a single ",
                        "locus rather than an average."), n_tr, min_size / 1e3))

  cat(sprintf("  %d transitions retained (compartments >= %g kb)\n",
              n_tr, min_size / 1e3))

  ## ---- 2. BOOTSTRAP CONFIDENCE BAND OVER THE TRANSITIONS -------------------
  cat(sprintf("Step 2 (%s): bootstrap (%d resamples over %d transitions)...\n",
              pca, n_boot, n_tr))

  set.seed(boot_seed)
  boot_idx <- lapply(seq_len(n_boot), function(i) sample.int(n_tr, n_tr, replace = TRUE))

  zB <- which(pos >= -zone_size & pos < 0)                   # B side
  zA <- which(pos > 0 & pos <= zone_size)                    # A side

  ## each profile brought to zero on its own B side, transition by transition
  tr_aligned <- lapply(tr_profiles, function(m)
    sweep(m, 2, colMeans(m[zB, , drop = FALSE], na.rm = TRUE), "-"))

  boot_of <- function(prof) {
    b <- rbindlist(lapply(conds, function(cd) {
      m  <- prof[[cd]]
      bm <- vapply(boot_idx, function(ix) rowMeans(m[, ix, drop = FALSE], na.rm = TRUE),
                   numeric(nrow(m)))
      data.table(Condition = cd,
                 Position  = pos,
                 Mean      = rowMeans(m, na.rm = TRUE),
                 Lo        = apply(bm, 1, quantile, probs = alpha,     names = FALSE),
                 Hi        = apply(bm, 1, quantile, probs = 1 - alpha, names = FALSE))
    }))
    b[, Condition := factor(Condition, levels = conds)][]
  }

  boot_al <- boot_of(tr_aligned)
  fwrite(boot_al, file.path(out_d, "Bootstrap_Profile_aligned.csv"))

  muts <- setdiff(conds, ref_cond)

  ## Small multiples: one panel per mutant, wild type recalled in grey behind.
  ## Six superimposed bands would be unreadable; here each comparison reads on
  ## its own.
  ref_dt   <- boot_al[Condition == ref_cond, .(Position, Mean, Lo, Hi)]
  facet_dt <- rbindlist(lapply(muts, function(mu) cbind(boot_al[Condition == mu], Facet = mu)))
  facet_rf <- rbindlist(lapply(muts, function(mu) cbind(copy(ref_dt), Facet = mu)))
  facet_dt[, Facet := factor(Facet, levels = muts)]
  facet_rf[, Facet := factor(Facet, levels = muts)]

  sub_al <- paste0("Band = ", round(100 * conf_level), " % bootstrap CI over ", n_tr,
                   " transitions (", n_boot, " resamplings). In grey: ", ref_cond,
                   ". Profiles set to 0 over ", pos[min(zB)] / 1e3, " to ",
                   pos[max(zB)] / 1e3, " kb (B side), transition by transition")

  g_prof <- ggplot() +
    geom_hline(yintercept = 0, colour = "#DDDDDD", linewidth = 0.4) +
    geom_ribbon(data = facet_rf, aes(x = Position / 1e3, ymin = Lo, ymax = Hi),
                fill = "grey60", alpha = 0.35) +
    geom_line(data = facet_rf, aes(x = Position / 1e3, y = Mean),
              colour = "grey30", linewidth = 0.5) +
    geom_ribbon(data = facet_dt, aes(x = Position / 1e3, ymin = Lo, ymax = Hi, fill = Condition),
                alpha = 0.35, show.legend = FALSE) +
    geom_line(data = facet_dt, aes(x = Position / 1e3, y = Mean, colour = Condition),
              linewidth = 0.7, show.legend = FALSE) +
    geom_vline(xintercept = 0, colour = "#CC3333", linewidth = 0.4, linetype = "dashed") +
    scale_fill_brewer(palette = "Dark2") +
    scale_colour_brewer(palette = "Dark2") +
    facet_wrap(~ Facet, ncol = 2) +
    labs(title = "B-to-A metaprofile aligned on the left baseline", subtitle = sub_al,
         x = "Distance to the transition (kb)", y = "Mean eigenvalue") +
    theme_bw() +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(face = "bold", size = 9),
          panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 8, colour = "#555555"))

  ## ---- 3. AMPLITUDE OF THE TRANSITION, AND ITS TWO TESTS -------------------
  ## No alignment is involved: the amplitude is a difference internal to each
  ## profile, so any constant offset cancels out.
  cat(sprintf("Step 3 (%s): amplitude tests...\n", pca))

  amp <- sapply(conds, function(cd) {
    m <- tr_profiles[[cd]]
    colMeans(m[zA, , drop = FALSE], na.rm = TRUE) -
      colMeans(m[zB, , drop = FALSE], na.rm = TRUE)
  })

  amp_ref <- mean(amp[, ref_cond])

  amp_stats <- rbindlist(lapply(conds, function(cd) {
    a  <- amp[, cd]
    t0 <- t.test(a)                                  # is there still a transition?
    w0 <- suppressWarnings(wilcox.test(a, mu = 0))
    pN <- if (cd == ref_cond) NA_real_ else t.test(a, amp[, ref_cond], paired = TRUE)$p.value
    data.table(
      Condition        = cd,
      n_transitions    = n_tr,
      Amplitude        = mean(a),
      CI_low           = t0$conf.int[1],
      CI_high          = t0$conf.int[2],
      Pct_of_WT        = 100 * mean(a) / amp_ref,
      CI_high_pct_WT   = 100 * t0$conf.int[2] / amp_ref,
      p_vs_zero        = t0$p.value,
      p_vs_zero_wilcox = w0$p.value,
      Transition       = ifelse(t0$p.value < 0.05, "detected", "not detectable"),
      p_vs_ref         = pN,
      Signif_vs_ref    = star(pN))
  }))
  amp_stats[, Signif_vs_zero := star(p_vs_zero)]
  print(amp_stats)
  fwrite(amp_stats, file.path(out_d, "Amplitude_Tests.csv"))

  amp_long <- data.table(
    Condition  = factor(rep(conds, each = n_tr), levels = conds),
    Chromosome = rep(tr_chr[[ref_cond]], times = length(conds)),
    Amplitude  = as.vector(amp))
  fwrite(amp_long, file.path(out_d, "Amplitude_Per_Transition.csv"))

  ## The jitter of the points is given an explicit seed so that the panel is
  ## byte-identical from one run to the next. It is a display offset on the x
  ## axis only: no value is modified, and the boxes are unaffected.
  g_amp <- ggplot(amp_long, aes(x = Condition, y = Amplitude, fill = Condition)) +
    geom_hline(yintercept = 0, colour = "#CC3333", linetype = "dashed", linewidth = 0.5) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.4, colour = "#222222") +
    geom_point(aes(colour = Condition),
               position = position_jitter(width = 0.12, height = 0, seed = jitter_seed),
               size = 1.6, alpha = 0.7, show.legend = FALSE) +
    scale_fill_brewer(palette = "Dark2") +
    scale_colour_brewer(palette = "Dark2") +
    labs(title = "Transition amplitude, transition by transition",
         subtitle = paste0("n = ", n_tr, " transitions"),
         x = NULL, y = "Amplitude") +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10))

  # SAVE ----
  ggsave(file.path(out_d, paste0(panels[[pca]][["profile"]],
                                 ".Metaprofile.Left_aligned.Small_multiples.png")),
         g_prof, width = 11, height = 7.5, dpi = 300)
  ggsave(file.path(out_d, paste0(panels[[pca]][["amplitude"]], ".Amplitude.Boxplot.png")),
         g_amp, width = 6, height = 5, dpi = 300)
}
