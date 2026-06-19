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


# ── SI Figure 1 | FlashLFQ Protein-Level DA: Diverse vs Limited PTMs ─────────
# Repeat of Figure 2 using FlashLFQ protein-level quantification
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

# Helpers ----

# Load QuantifiedProteins.tsv and rename its generic Intensity_Default_N columns
# to the named sample columns from the corresponding peptide TSV.
LoadProteinData <- function(protein_path, peptide_path) {
  pep_header          <- read.csv(peptide_path, sep = "\t", nrows = 0)
  pep_intensity_cols  <- grep("^Intensity_", colnames(pep_header), value = TRUE)
  df                  <- read.csv(protein_path, sep = "\t", row.names = NULL)
  prot_intensity_cols <- grep("^Intensity_Default_", colnames(df), value = TRUE)
  stopifnot(length(prot_intensity_cols) == length(pep_intensity_cols))
  colnames(df)[match(prot_intensity_cols, colnames(df))] <- pep_intensity_cols
  df[grepl("Homo sapiens", df$Organism), ]
}

# Build log2 expression matrix (rows = proteins, columns = samples)
GetProteinExprMatrix <- function(df, intensity_cols) {
  expr <- apply(df[, intensity_cols], 2, as.numeric)
  rownames(expr) <- df$Protein.Groups
  expr[expr == 0] <- NA
  log2(expr)
}

# Limma fit with donor blocking and array quality weights.
# Mirrors dpcDE(sample.weights = TRUE, block = Donor) in the main analysis.
RunLimmaDA <- function(expr, targets) {
  design <- model.matrix(~0 + Group, data = targets)
  colnames(design) <- gsub("Group", "", colnames(design))

  aw     <- arrayWeights(expr, design)
  corfit <- duplicateCorrelation(expr, design, block = targets$Donor, weights = aw)

  fit  <- lmFit(expr, design,
                block = targets$Donor, correlation = corfit$consensus,
                weights = aw)
  fit2 <- contrasts.fit(fit, makeContrasts(ALS_vs_CTRL = "ALS - CTRL", levels = design))
  eBayes(fit2)
}

# Full pipeline: load → filter samples → metadata → filter proteins → limma
RunFlashLFQPipeline <- function(prot_path, pep_path, min_obs_proteins = 300) {
  df             <- LoadProteinData(prot_path, pep_path)
  intensity_cols <- grep("^Intensity_Biogen_", colnames(df), value = TRUE)

  expr <- GetProteinExprMatrix(df, intensity_cols)

  # Drop low-coverage samples (mirrors the >2500 peptide filter in main analysis)
  expr    <- expr[, colSums(!is.na(expr)) > min_obs_proteins]
  targets <- GetMetaDataTdpStratified(expr)
  expr    <- expr[, colnames(expr) %in% targets$Sample]

  # Drop proteins quantified in fewer than 3 samples
  expr <- expr[rowSums(!is.na(expr)) >= 3, ]
  stopifnot(all(targets$Sample == colnames(expr)))

  eb_fit           <- RunLimmaDA(expr, targets)
  tbl              <- topTable(eb_fit, coef = "ALS_vs_CTRL", number = Inf)
  tbl$ProteinGroup <- rownames(tbl)

  gene_map <- df[!duplicated(df$Protein.Groups), c("Protein.Groups", "Gene.Name")]
  result   <- merge(gene_map, tbl, by.x = "Protein.Groups", by.y = "ProteinGroup")
  result$Gene <- CleanGeneNames(result$Gene.Name)
  result
}

# Run analyses ----
message("Running FlashLFQ limma — diverse PTMs...")
table_protein_diverse_flashlfq <- RunFlashLFQPipeline(PATH_FLASHLFQ_DIVERSE,
                                                       "Data/DiversePtms/QuantifiedPeptides.tsv") |>
  label_de()

message("Running FlashLFQ limma — limited PTMs...")
table_protein_limited_flashlfq <- RunFlashLFQPipeline(PATH_FLASHLFQ_LIMITED,
                                                       "Data/LimitedPtms/QuantifiedPeptides.tsv") |>
  label_de()

message(sprintf("FlashLFQ diverse DE proteins: %d",
                sum(table_protein_diverse_flashlfq$DEStatus != "NotDE")))
message(sprintf("FlashLFQ limited DE proteins: %d",
                sum(table_protein_limited_flashlfq$DEStatus != "NotDE")))

# Merged comparison table ----
flashlfq_diverse_de <- table_protein_diverse_flashlfq$Gene[table_protein_diverse_flashlfq$DEStatus != "NotDE"]
flashlfq_limited_de <- table_protein_limited_flashlfq$Gene[table_protein_limited_flashlfq$DEStatus != "NotDE"]

merged_pro_flashlfq <- merge(
  table_protein_diverse_flashlfq,
  table_protein_limited_flashlfq[, c("Gene", "logFC", "adj.P.Val")],
  by = "Gene", suffixes = c("", "_NoMod")
) %>%
  mutate(
    LogFC_Diff = logFC - logFC_NoMod,
    ModDiff = case_when(
      Gene %in% setdiff(flashlfq_limited_de, flashlfq_diverse_de) ~ "NoModDA",
      Gene %in% setdiff(flashlfq_diverse_de, flashlfq_limited_de) ~ "ModDA",
      TRUE ~ "Shared"
    )
  ) %>%
  arrange(desc(ModDiff))

# Write outputs ----
write.table(table_protein_diverse_flashlfq,
            file.path("Supplemental", "LimmaProteinDA_Diverse.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(table_protein_limited_flashlfq,
            file.path("Supplemental", "LimmaProteinDA_Limited.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(merged_pro_flashlfq,
            file.path("Supplemental", "LimmaProteinDA_Comparison.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote FlashLFQ DA results to Supplemental/")


# ── SI Figure 2 | TDP-43 Stratification Clustering ───────────────────────────
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


# ── SI Figure 3 | Carboxymethylation Analysis ─────────────────────────────────
# Code exists in AlsMotorNeuronAnalysis — locate and migrate here.
# TODO: source or paste CML analysis code once identified.
