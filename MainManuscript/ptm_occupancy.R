# ptm_occupancy.R — MainManuscript
# Alternative to the diffSplice approach in analysis.R.
# Normalizes each modified peptide's intensity to its parent protein
# (PTM occupancy proxy), then tests for differential abundance with
# standard limma/dpcDE.
#
# Requires analysis.R to have been sourced first (supplies y.peptide,
# y.peptide.imputed, y.protein, design, targets).
#
# Future work: Satterthwaite approximation to propagate uncertainty from
# the combined SE stored in $other$standard.error.

library(tidyverse)
library(limma)
library(limpa)

source("MainManuscript/CustomScripts.R")

PROP_OBS_MIN <- 0.1


# ── Helper functions (from CustomScripts.R, redefined here for isolation) ─────

NormalizeExpressionToPeptides <- function(y.peptide, mod_seq, base_seq, pg) {
  modified_row  <- which(y.peptide$genes$PeptideSequence == mod_seq)
  modified_expr <- y.peptide$E[modified_row, ]

  if (is.null(pg)) return(rep(0, ncol(y.peptide$E)))

  control_peps <- GetControlPeptides(y.peptide, pg, mod_seq)
  if (is.null(control_peps)) return(rep(0, ncol(y.peptide$E)))

  unmodified_expr <- colMeans(control_peps$E)

  stopifnot(all(y.peptide$targets$Sample == control_peps$targets$Sample))

  # Ratio in log-space: (modified - unmod), anchored to median(unmod)
  normed_expr <- median(unmodified_expr) + (modified_expr - unmodified_expr)
  return(normed_expr)
}

# Returns list(expr, se) so GetNormalizedEList can propagate uncertainty.
# se = sqrt(SE_pep^2 + SE_prot^2): subtraction of independent normals adds variances.
# Requires y.peptide.imputed$other$standard.error and y.protein$other$standard.error.
NormalizeExpressionToProteins <- function(peptide_e, protein_e, mod_seq, pg) {
  zero <- list(expr = rep(0, ncol(peptide_e$E)), se = NULL)

  modified_row  <- which(peptide_e$genes$PeptideSequence == mod_seq)
  modified_expr <- peptide_e$E[modified_row, ]
  modified_se   <- peptide_e$other$standard.error[modified_row, ]

  if (is.null(pg)) return(zero)

  protein_row <- which(protein_e$genes$ProteinGroup == pg)
  if (length(protein_row) == 0) return(zero)

  protein_expr <- protein_e$E[protein_row, ] %>% as.vector()
  protein_se   <- protein_e$other$standard.error[protein_row, ] %>% as.vector()

  stopifnot(all(peptide_e$targets$Sample == y.protein$targets$Sample))

  normed_expr <- median(protein_expr) + (modified_expr - protein_expr)
  combined_se <- sqrt(modified_se^2 + protein_se^2)

  return(list(expr = normed_expr, se = combined_se))
}

GetNormalizedEList <- function(peptides_of_interest, y.peptide,
                               y.protein = NULL, normalize_to_protein = TRUE) {
  mod.peptide <- y.peptide[y.peptide$genes$PeptideSequence %in% peptides_of_interest, ]

  if (!normalize_to_protein) {
    # Returns plain vectors; apply() produces an [nSamples × nPeptides] matrix
    raw <- apply(mod.peptide$genes, MARGIN = 1, FUN = function(x)
      NormalizeExpressionToPeptides(y.peptide,
                                    x[["PeptideSequence"]],
                                    x[["BaseSequence"]],
                                    x[["ProteinGroup"]]))
    normed_elist <- if (is.list(raw)) do.call(rbind, raw) else t(raw)
    normed_se    <- NULL

  } else if (!is.null(y.protein)) {
    filtered_protein <- y.protein[y.protein$genes$NPeptides >= 4, ]
    # lapply over row indices so each call returns list(expr, se)
    results <- lapply(seq_len(nrow(mod.peptide$genes)), function(i) {
      x <- mod.peptide$genes[i, , drop = FALSE]
      NormalizeExpressionToProteins(y.peptide, filtered_protein,
                                    x[["PeptideSequence"]], x[["ProteinGroup"]])
    })
    normed_elist <- do.call(rbind, lapply(results, `[[`, "expr"))
    se_list      <- lapply(results, `[[`, "se")
  }

  keep <- !apply(normed_elist, 1, function(x) all(x == 0))

  normed.peptide   <- mod.peptide[keep, ]
  normed.peptide$E <- normed_elist[keep, , drop = FALSE]

  if (exists("se_list")) {
    se_kept <- se_list[keep]
    if (all(!sapply(se_kept, is.null)))
      normed.peptide$other$standard.error <- do.call(rbind, se_kept)
  }

  return(normed.peptide)
}


# ── PTM occupancy: normalize to parent protein ────────────────────────────────

all.mod.normed.pro <- GetNormalizedEList(
  y.peptide.imputed$genes$PeptideSequence,
  y.peptide.imputed,
  y.protein            = y.protein,
  normalize_to_protein = TRUE
)

PROP_OBS_MIN <- 0.12
all.mod.normed.pro <- all.mod.normed.pro[
  all.mod.normed.pro$genes$PropObs >= PROP_OBS_MIN, ]


# ── Differential abundance on occupancy-normalized intensities ────────────────

fit_occ    <- dpcDE(all.mod.normed.pro, block = targets$Donor, design = design, plot = TRUE)
contrast   <- makeContrasts(ALS_vs_CTRL = "ALS - CTRL", levels = design)
fit_occ2   <- contrasts.fit(fit_occ, contrast)
eb_fit_occ <- eBayes(fit_occ2)

table_occ <- topTable(eb_fit_occ, coef = "ALS_vs_CTRL", number = Inf)

logFC_cutoff <- 1.5
pval_cutoff  <- 0.05

table_occ$DEStatus <- case_when(
  table_occ$adj.P.Val < pval_cutoff & table_occ$logFC >=  logFC_cutoff ~ "UP",
  table_occ$adj.P.Val < pval_cutoff & table_occ$logFC <= -logFC_cutoff ~ "DOWN",
  TRUE ~ "NotDE"
)

table(table_occ$DEStatus)

# What changed

all.mod.normed.classic <- all.mod.normed.classic[
  all.mod.normed.classic$genes$PropObs >= PROP_OBS_MIN, ]

all.mod.normed.new$E[1:5, 1:5] - all.mod.normed.classic$E[1:5, 1:5]

all.mod.normed.new$other$standard.error[1:5, 1:5] - all.mod.normed.classic$other$standard.error[1:5, 1:5]

all.mod.normed.new$genes$ProteinGroup %>% head()
all.mod.normed.classic$genes$PeptideSequence %>% head()

all.mod.normed.new$E %>% dim()
all.mod.normed.classic$E %>% dim()
