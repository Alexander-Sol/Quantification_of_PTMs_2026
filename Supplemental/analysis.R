# analysis.R — Supplemental
# Analyses that feed SI Figures 1–3.
# Run top-to-bottom before sourcing Supplemental/plotting.R.
# Assumes MainManuscript/analysis.R objects are already in the environment
# (y.peptide, y.protein, targets, table_protein_diverse, etc.).

library(tidyverse)
library(limpa)
library(limma)

source("MainManuscript/CustomScripts.R")

PROP_OBS_MIN <- 0.1


# ── SI Figure 1 | TDP-43 Stratification Clustering ───────────────────────────
# Goal: show no meaningful clustering by TDP-43 level and minimal
# protein-level changes across stratification groups.
# Panel A: PCA of all samples, coloured by TDP-43 level.
# Panel B: Volcano of within-ALS TDP-43 severity trend (should show ~0 DA proteins).

# PCA — use proteins observed in >= 50% of samples to keep it robust
min_obs <- ceiling(ncol(y.protein$E) * 0.5)
prot_keep <- rowSums(!is.na(y.protein$E)) >= min_obs
pca_mat   <- y.protein$E[prot_keep, ]

# mean-impute any remaining NAs so prcomp doesn't drop samples
pca_mat <- apply(pca_mat, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE); x
}) %>% t()

pca_res  <- prcomp(t(pca_mat), scale. = TRUE, center = TRUE)
pca_var  <- summary(pca_res)$importance["Proportion of Variance", 1:2] * 100

pca_df <- as.data.frame(pca_res$x[, 1:2]) %>%
  tibble::rownames_to_column("Sample") %>%
  left_join(targets[, c("Sample", "Group", "Donor", "TDP_Lvl")], by = "Sample")

# Within-ALS TDP-43 severity trend (limma, linear contrast over ordered levels)
tdp_order <- c("NON", "MLD", "MOD", "SEV")

als_idx       <- targets$Group == "ALS"
y.protein.als <- y.protein[, als_idx]
targets.als   <- y.protein.als$targets
targets.als$TDP_Lvl <- factor(targets.als$TDP_Lvl, levels = tdp_order, ordered = TRUE)

design_tdp <- model.matrix(~ TDP_Lvl + PMI + Sex, data = targets.als)
fit_tdp    <- lmFit(y.protein.als$E, design_tdp)
fit_tdp    <- eBayes(fit_tdp)
table_tdp  <- topTable(fit_tdp, coef = "TDP_Lvl.L", number = Inf, sort.by = "P")


# ── SI Figure 2 | Carboxymethylation Analysis ─────────────────────────────────
# Code exists in AlsMotorNeuronAnalysis — locate and migrate here.
# TODO: source or paste CML analysis code once identified.


# ── SI Figure 3 | FlashLFQ Protein-Level DA: Diverse vs Limited PTMs ─────────
# Repeat of Figure 4 using FlashLFQ protein-level quantification
# to demonstrate the diverse vs limited PTM pattern is method-agnostic.

PATH_FLASHLFQ_DIVERSE <- "Data/DiversePtms/QuantifiedProteins.tsv"
PATH_FLASHLFQ_LIMITED <- "Data/LimitedPtms/QuantifiedProteins.tsv"

logFC_cutoff <- 1.5
pval_cutoff  <- 0.05

label_de <- function(tbl) {
  tbl %>% mutate(DEStatus = case_when(
    adj.P.Val < pval_cutoff & logFC >=  logFC_cutoff ~ "UP",
    adj.P.Val < pval_cutoff & logFC <= -logFC_cutoff ~ "DOWN",
    TRUE ~ "NotDE"
  ))
}

# TODO: load FlashLFQ protein TSVs, run limma/dpcDE on the protein-level
# quantification, and populate table_protein_diverse_flashlfq and
# table_protein_limited_flashlfq analogously to MainManuscript/analysis.R.
#
# Skeleton:
# flashlfq_div <- read.csv(PATH_FLASHLFQ_DIVERSE, sep = "\t")
# ... (filter, build EList, run dpcDE) ...
# table_protein_diverse_flashlfq <- topTable(...) %>% label_de()
#
# flashlfq_lim <- read.csv(PATH_FLASHLFQ_LIMITED, sep = "\t")
# ... (filter, build EList, run dpcDE) ...
# table_protein_limited_flashlfq <- topTable(...) %>% label_de()
#
# merged_pro_flashlfq <- merge(table_protein_diverse_flashlfq,
#                              table_protein_limited_flashlfq[, c("Gene","logFC","adj.P.Val")],
#                              by = "Gene", suffixes = c("", "_NoMod")) %>%
#   mutate(ModDiff = ...)
