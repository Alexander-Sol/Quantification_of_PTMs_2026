# Protein-specific peptide heatmap export (DiffSplice → JavaTreeView)
#
# Assumes y.peptide and y.protein EList objects are already in the environment.
# Set the parameters below, then run the whole script.

library(tidyverse)
library(limpa)

source("C:/Users/Alex/Source/Repos/Quantification_of_PTMs_2026/MainManuscript/CustomScripts.R")
source("C:/Users/Alex/Source/Repos/AlsMotorNeuronAnalysis/Supplemental/topSplice.R")
source("C:/Users/Alex/Source/Repos/AlsMotorNeuronAnalysis/Tdp43_Stratified/HclustExport.R")

# ── Parameters ────────────────────────────────────────────────────────────────

protein_of_interest <- "NEFH"      # gene name to export

out_dir <- r"(C:\Users\Alex\Source\Repos\AlsMotorNeuronAnalysis\Tdp43_Stratified)"

# Samples to drop before clustering (set to character(0) to keep all)
outlier_samples <- c("4A3", "4A6", "4A4", "4A9")

# diffSplice filtering thresholds
treat_lfc   <- 0.4    # treat.lfc passed to topSplice
fdr_broad   <- 0.05   # used to build table_splice background
prop_broad  <- 0.1
fdr_strict  <- 0.01   # used to select peptides for this protein
prop_strict <- 0.09
lfc_strict  <- 1.0

# RDS cache paths (saves re-running expensive steps)
y.peptide.imputed_rds_path <- file.path(out_dir, "y.peptide.imputed.rds")

# ── Helper functions ───────────────────────────────────────────────────────────

NormalizeToControl <- function(eList) {
  ctrl <- eList$targets$Group == "CTRL"
  eList$E <- sweep(eList$E, 1, rowMeans(eList$E[, ctrl], na.rm = TRUE), "-")
  eList
}

RemoveOutliers <- function(eList, outlier_samples) {
  keep <- !sapply(eList$targets$Sample,
                  function(s) strsplit(s, "_")[[1]][[3]] %in% outlier_samples)
  eList$E       <- eList$E[, keep, drop = FALSE]
  eList$targets <- eList$targets[keep, , drop = FALSE]
  eList
}

annotate_categories <- function(tbl) {
  category_levels <- c("Biological Mod", "Carboxymethylation", "Other Mod",
                       "Carbamidomethylation", "Unmodified")
  tbl$AllMod   <- GetModsTextWReplacement(tbl$PeptideSequence, all_mods = TRUE) |>
    sapply(paste)
  tbl$Category <- ClassifyPeptideModsSpecific(tbl$AllMod)
  tbl$Category <- factor(tbl$Category, levels = category_levels)
  tbl
}

# ── DPC imputation (cached) ───────────────────────────────────────────────────

if (!exists("y.peptide.imputed") || is.null(y.peptide.imputed)) {
  if (file.exists(y.peptide.imputed_rds_path)) {
    y.peptide.imputed <- readRDS(y.peptide.imputed_rds_path)
  } else {
    dpcfit <- dpc(y.peptide)
    y.peptide.imputed <- dpcImpute(y.peptide, "PeptideSequence", dpc = dpcfit, verbose = TRUE)
    saveRDS(y.peptide.imputed, y.peptide.imputed_rds_path, compress = FALSE)
  }
}

# ── diffSplice analysis ───────────────────────────────────────────────────────

# targets <- y.peptide$targets
# 
# design <- model.matrix(~0 + Group + PMI + Sex, data = targets)
# colnames(design) <- gsub("Group", "", colnames(design))
# colnames(design)[ncol(design)] <- "Sex"
# 
# fit      <- dpcDE(y.peptide.imputed, design, plot = FALSE)
# fit2     <- contrasts.fit(fit, makeContrasts(ALS_vs_CTRL = "ALS - CTRL", levels = design))
# eb_fit   <- eBayes(fit2)
# 
# du <- diffSplice(eb_fit, geneid = "Gene", exonid = "Sequence")
# 
# table_splice <- topSplice(du, coef = "ALS_vs_CTRL", test = "t",
#                            treat.lfc = treat_lfc, number = Inf, FDR = fdr_broad)
# table_splice <- table_splice[table_splice$PropObs >= prop_broad, ]
# table_splice <- annotate_categories(table_splice)
# 
# # ── Protein-specific peptide selection ────────────────────────────────────────
# 
# protein_peptides <- table_splice |>
#   filter(Gene == protein_of_interest,
#          FDR      <= fdr_strict,
#          PropObs  >= prop_strict,
#          abs(logFC) >= lfc_strict) |>
#   pull(PeptideSequence)
# 
# if (length(protein_peptides) == 0)
#   stop("No peptides passed the filters for '", protein_of_interest,
#        "'. Try relaxing fdr_strict, prop_strict, or lfc_strict.")
# 
# # ── Normalization ─────────────────────────────────────────────────────────────

protein_peptides <- table_occ$PeptideSequence[table_occ$Gene == protein_of_interest &
                                                table_occ$PropObs >= 0.12]


protein_data <- GetNormalizedEList(protein_peptides, y.peptide.imputed, y.protein,
                                   normalize_to_protein = TRUE)

if (length(outlier_samples) > 0)
  protein_data <- RemoveOutliers(protein_data, outlier_samples)

protein_data <- NormalizeToControl(protein_data)

# Order columns by donor and rename to Donor-N
protein_ordered <- protein_data[, order(protein_data$targets$Donor)]
donor_table <- table(protein_ordered$targets$Donor)
colnames(protein_ordered) <- unlist(sapply(names(donor_table), function(d)
  paste0(d, "-", seq_len(donor_table[d]))))

protein_ordered$genes <- AddPositionColumns(protein_ordered$genes)

# Subset to unmodified / limited-search peptides only
pep_category  <- ClassifyPeptides(protein_ordered$genes$PeptideSequence)
keep_unmod    <- pep_category %in% c("Unmodified", "Limited-Search Mod")
protein_unmod <- protein_ordered
protein_unmod$E     <- protein_unmod$E[keep_unmod, , drop = FALSE]
protein_unmod$genes <- protein_unmod$genes[keep_unmod, , drop = FALSE]

# ── Export ────────────────────────────────────────────────────────────────────

export_one <- function(elist, label) {
  prefix_cols   <- file.path(out_dir, paste0(protein_of_interest, "_", label, "_HeatmapCols"))
  prefix_nocols <- file.path(out_dir, paste0(protein_of_interest, "_", label, "_HeatmapNoCols"))

  result_cols <- ExportToJavaTreeView(
    elist        = elist,
    out_prefix   = prefix_cols,
    name_col     = "PeptideSequence",
    yorf_col     = "ProteinGroup",
    cluster_cols = TRUE
  )
  result_nocols <- ExportToJavaTreeView(
    elist        = elist,
    out_prefix   = prefix_nocols,
    name_col     = "PeptideSequence",
    yorf_col     = "ProteinGroup",
    cluster_cols = FALSE
  )
  PlotClusteredHeatmap(result_cols,   main = paste(protein_of_interest, label, "(columns clustered)"))
  PlotClusteredHeatmap(result_nocols, main = paste(protein_of_interest, label, "(columns fixed)"))
}

export_one(protein_ordered, "AllPep")
export_one(protein_unmod,   "UnmodLimited")
