# Protein-level differential abundance analysis using limma
# Analyzes pre-quantified QuantifiedProteins.tsv files (no peptide rollup needed).
# Runs with and without modified peptides considered, then compares the two.

library(tidyverse)
library(limma)
library(ggrepel)

source("Supplemental/CustomScripts.R")

# Paths ----
pep_path_mod    <- "Tdp43_Stratified/QuantData/DiversePtms/QuantifiedPeptides.tsv"
pep_path_nomod  <- "Tdp43_Stratified/QuantData/LimitedPtms/QuantifiedPeptides.tsv"
prot_path_mod   <- "Tdp43_Stratified/QuantData/DiversePtms/QuantifiedProteins.tsv"
prot_path_nomod <- "Tdp43_Stratified/QuantData/LimitedPtms/QuantifiedProteins.tsv"

logFC_cutoff <- 1.5
pval_cutoff  <- 0.05
text_size    <- 1.1

# Helpers ----

# Load QuantifiedProteins.tsv and rename its generic Intensity_Default_N columns
# to match the named Intensity_Biogen_<SampleID>_... columns from the peptide file.
# FlashLFQ writes proteins with numbered columns in the same sample order as peptides.
LoadProteinData <- function(protein_path, peptide_path) {
  pep_header        <- read.csv(peptide_path, sep = '\t', nrows = 0)
  pep_intensity_cols <- grep("^Intensity_", colnames(pep_header), value = TRUE)

  df <- read.csv(protein_path, sep = '\t', row.names = NULL)
  prot_intensity_cols <- grep("^Intensity_Default_", colnames(df), value = TRUE)
  stopifnot(length(prot_intensity_cols) == length(pep_intensity_cols))
  colnames(df)[match(prot_intensity_cols, colnames(df))] <- pep_intensity_cols

  df <- df[grepl("Homo sapiens", df$Organism), ]
  return(df)
}

# Build log2 protein expression matrix (rownames = Protein.Groups)
GetProteinExprMatrix <- function(df, intensity_cols) {
  expr <- apply(df[, intensity_cols], 2, as.numeric)
  rownames(expr) <- df$Protein.Groups
  expr[expr == 0] <- NA
  expr <- log2(expr)
  return(expr)
}

# Limma fit with donor blocking and sample quality weights
RunLimmaDA <- function(expr, targets) {
  design <- model.matrix(formula(~0 + Group + PMI + Sex), data = targets)
  colnames(design) <- gsub("Group", "", colnames(design))
  colnames(design)[ncol(design)] <- "Sex"

  # Estimate per-sample quality weights (equivalent to sample.weights=T in dpcDE)
  aw <- arrayWeights(expr, design)

  # Estimate inter-donor correlation for repeated-measures blocking
  corfit <- duplicateCorrelation(expr, design, block = targets$Donor, weights = aw)

  fit <- lmFit(expr, design,
               block = targets$Donor, correlation = corfit$consensus,
               weights = aw)

  contrast.matrix <- makeContrasts(ALS_vs_CTRL = "ALS - CTRL", levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  eBayes(fit2)
}

# Full pipeline: load → filter → metadata → limma → annotated results
RunFullPipeline <- function(prot_path, pep_path, min_obs_proteins = 300) {
  df             <- LoadProteinData(prot_path, pep_path)
  intensity_cols <- grep("^Intensity_", colnames(df), value = TRUE)

  expr <- GetProteinExprMatrix(df, intensity_cols)

  # Drop samples with fewer than min_obs_proteins quantified proteins.
  # (Analogous to the >2500 peptide filter in the peptide-level workflow;
  #  adjust threshold based on your observed per-sample coverage.)
  col_keep <- colSums(!is.na(expr)) > min_obs_proteins
  expr     <- expr[, col_keep]

  targets <- GetMetaDataTdpStratified(expr)
  expr    <- expr[, colnames(expr) %in% targets$Sample]

  # Drop proteins observed in fewer than 3 samples
  keep <- rowSums(!is.na(expr)) >= 3
  expr <- expr[keep, ]

  stopifnot(all(targets$Sample == colnames(expr)))

  eb_fit       <- RunLimmaDA(expr, targets)
  tbl          <- topTable(eb_fit, coef = "ALS_vs_CTRL", number = Inf)
  tbl$ProteinGroup <- rownames(tbl)

  # Merge gene names from the protein annotation columns
  unique_df <- df[!duplicated(df$Protein.Groups), c("Protein.Groups", "Gene.Name")]
  result    <- merge(unique_df, tbl, by.x = "Protein.Groups", by.y = "ProteinGroup")
  result$Gene.Name <- CleanGeneNames(result$Gene.Name)
  result
}

AddDEStatus <- function(tbl) {
  tbl$DEStatus <- ifelse(
    tbl$adj.P.Val < pval_cutoff & tbl$logFC >= logFC_cutoff,  "UP",
    ifelse(tbl$adj.P.Val < pval_cutoff & tbl$logFC <= -logFC_cutoff, "DOWN", "NotDE")
  )
  tbl
}

MakeVolcano <- function(tbl, title, genes_highlight = NULL, x_limit = c(-8, 8)) {
  p <- ggplot(tbl, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = DEStatus), alpha = 0.6, size = 2, shape = 19) +
    scale_color_manual(values = c("NotDE" = "grey", "UP" = "red", "DOWN" = "red")) +
    geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(pval_cutoff),            linetype = "dashed", color = "black") +
    xlab("Log2 Fold Change") + ylab("-Log10 Adjusted P-value") +
    ggtitle(title) +
    theme_minimal() +
    theme(plot.title   = element_text(size = 16 * text_size),
          axis.title.x = element_text(size = 14 * text_size),
          axis.title.y = element_text(size = 14 * text_size),
          axis.text.x  = element_text(size = 14 * text_size),
          axis.text.y  = element_text(size = 12 * text_size),
          legend.position = "none") +
    xlim(x_limit)

  if (!is.null(genes_highlight)) {
    highlight <- subset(tbl, Gene.Name %in% genes_highlight)
    p <- p +
      geom_point(data = highlight, aes(x = logFC, y = -log10(adj.P.Val)),
                 color = "darkred", size = 4, shape = 19) +
      geom_text_repel(data = highlight,
                      aes(x = logFC, y = -log10(adj.P.Val), label = Gene.Name),
                      color = "black", size = 5, fontface = "bold")
  }
  p
}

# Run analyses ----
results_mod   <- RunFullPipeline(prot_path_mod,   pep_path_mod)   |> AddDEStatus()
results_nomod <- RunFullPipeline(prot_path_nomod, pep_path_nomod) |> AddDEStatus()

table(results_mod$DEStatus)
table(results_nomod$DEStatus)

# Write outputs ----
write_tsv(results_mod,   "Tdp43_Stratified/LimmaProteinDA_DiversePTMs.tsv")
write_tsv(results_nomod, "Tdp43_Stratified/LimmaProteinDA_LimitedPTMs.tsv")

# Individual volcano plots ----
genes_of_interest <- c("NEFH", "NEFM", "NEFL", "GFAP",
                       "YWHAB", "YWHAG", "YWHAE", "YWHAH", "YWHAQ", "YWHAZ")

p_vol_mod   <- MakeVolcano(results_mod,   "DA Proteins (Diverse PTMs)",  genes_of_interest)
p_vol_nomod <- MakeVolcano(results_nomod, "DA Proteins (Limited PTMs)",  genes_of_interest)
print(p_vol_mod)
print(p_vol_nomod)
ggsave("Tdp43_Stratified/LimmaProteinVolcano_DiversePTMs.png", p_vol_mod,   width = 6, height = 5, dpi = 300)
ggsave("Tdp43_Stratified/LimmaProteinVolcano_LimitedPTMs.png", p_vol_nomod, width = 6, height = 5, dpi = 300)

# Comparison: mod vs no-mod ----
mod_genes   <- results_mod$Gene.Name[results_mod$DEStatus   != "NotDE"]
nomod_genes <- results_nomod$Gene.Name[results_nomod$DEStatus != "NotDE"]

no_mod_da <- setdiff(nomod_genes, mod_genes)  # DA only in the no-mod analysis
mod_da    <- setdiff(mod_genes, nomod_genes)  # DA only when mods are included

comparison <- merge(results_mod,
                    results_nomod[, c("Gene.Name", "logFC", "adj.P.Val")],
                    by = "Gene.Name",
                    suffixes = c("", "_NoMod"))

comparison$LogFC_Diff   <- comparison$logFC - comparison$logFC_NoMod
comparison$ModDiff      <- ifelse(comparison$Gene.Name %in% no_mod_da, "NoModDA",
                            ifelse(comparison$Gene.Name %in% mod_da,   "ModDA", "Shared"))
comparison$Directionality <- ifelse(
  ((comparison$logFC > 0 & comparison$logFC_NoMod < 0) |
   (comparison$logFC < 0 & comparison$logFC_NoMod > 0)) &
    comparison$adj.P.Val      < pval_cutoff &
    comparison$adj.P.Val_NoMod < pval_cutoff,
  "Opposite", "Same"
)
comparison <- comparison[order(comparison$ModDiff, decreasing = TRUE), ]

genes_compare <- c("NEFM", "DRG1", "CLASP1", "EIF3J", "MARCKSL1", "EEF1A1", "USP9X")

# Comparison volcano (mod vs no-mod coloring) ----
highlight <- subset(comparison, Gene.Name %in% genes_compare)

p_moddiff_volcano <- ggplot(comparison, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = ModDiff), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("Shared" = "grey80", "ModDA" = "forestgreen", "NoModDA" = "darkred")) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pval_cutoff),            linetype = "dashed", color = "black") +
  xlab("Log2 Fold Change") + ylab("-Log10 Adjusted P-value") +
  ggtitle("Changes in Differential Abundance when\nModifications are Considered") +
  theme_minimal() +
  theme(plot.title   = element_text(size = 16 * text_size),
        axis.title.x = element_text(size = 14 * text_size),
        axis.title.y = element_text(size = 14 * text_size),
        axis.text.x  = element_text(size = 14 * text_size),
        axis.text.y  = element_text(size = 12 * text_size),
        legend.position = "none") +
  geom_text_repel(data = highlight,
                  aes(x = logFC, y = -log10(adj.P.Val), label = Gene.Name),
                  color = "black", size = 5, fontface = "bold", box.padding = 0.5)
print(p_moddiff_volcano)
ggsave("Tdp43_Stratified/LimmaProteinVolcano_ModDiff.png", p_moddiff_volcano, width = 6, height = 5, dpi = 300)

# logFC scatter: with mods vs without mods ----
highlight <- subset(comparison, Gene.Name %in% genes_compare)

p_logfc_scatter <- ggplot(comparison, aes(x = logFC, y = logFC_NoMod)) +
  geom_point(aes(color = ModDiff), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("Shared" = "grey80", "ModDA" = "forestgreen", "NoModDA" = "darkred")) +
  geom_hline(yintercept = 0, color = "grey20", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "grey20", linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  xlab("Log2 Fold Change (Diverse PTMs)") + ylab("Log2 Fold Change (Limited PTMs)") +
  ggtitle("Comparison of Protein Log2FC With\nand Without Modified Peptides") +
  theme_minimal() +
  theme(plot.title   = element_text(size = 16 * text_size),
        axis.title.x = element_text(size = 14 * text_size),
        axis.title.y = element_text(size = 14 * text_size),
        axis.text.x  = element_text(size = 14 * text_size),
        axis.text.y  = element_text(size = 12 * text_size),
        legend.position = "none") +
  geom_text_repel(data = highlight,
                  aes(x = logFC, y = logFC_NoMod, label = Gene.Name),
                  color = "black", size = 5, fontface = "bold",
                  box.padding = 1.5, point.padding = 0.5,
                  min.segment.length = 0, segment.color = "black")
print(p_logfc_scatter)
ggsave("Tdp43_Stratified/LimmaProteinLogFC_Scatter.png", p_logfc_scatter, width = 6, height = 5, dpi = 300)


# ── DE Overlap: Diverse vs Limited PTMs ────────────────────────────────────----
n_diverse_only <- length(setdiff(mod_genes, nomod_genes))
n_limited_only <- length(setdiff(nomod_genes, mod_genes))
n_both         <- length(intersect(mod_genes, nomod_genes))

cat("DE protein counts:\n")
cat(sprintf("  Diverse PTMs only : %d\n", n_diverse_only))
cat(sprintf("  Limited PTMs only : %d\n", n_limited_only))
cat(sprintf("  Shared (both)     : %d\n", n_both))
cat(sprintf("  Total diverse     : %d\n", length(mod_genes)))
cat(sprintf("  Total limited     : %d\n", length(nomod_genes)))

# Area-proportional Euler diagram (circles scale with set size)
library(eulerr)
fit_euler <- euler(list("Diverse PTMs" = mod_genes, "Limited PTMs" = nomod_genes))

p_venn <- plot(fit_euler,
               quantities = list(cex = 1.3),
               labels     = list(cex = 1.1),
               fills      = list(fill = c("forestgreen", "darkred"), alpha = 0.35),
               edges      = list(col = "grey20", lwd = 1.5),
               main       = "DE Proteins: Diverse vs. Limited PTMs")
print(p_venn)
png("Tdp43_Stratified/LimmaProteinDA_VennDiagram.png", width = 5, height = 4, units = "in", res = 300)
print(p_venn)
dev.off()
