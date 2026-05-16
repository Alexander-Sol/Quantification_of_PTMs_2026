# plotting.R — Supplemental
# All SI figures. Requires Supplemental/analysis.R (and MainManuscript/analysis.R)
# to have been sourced first.
# Figures are saved to Supplemental/Figures/.

library(tidyverse)
library(cowplot)
library(scales)
library(ggrepel)
library(eulerr)

source("MainManuscript/CustomScripts.R")

FIG_DIR   <- "Supplemental/Figures"
TEXT_SIZE <- 1.1


# ── SI Figure 1 | TDP-43 Stratification Clustering ───────────────────────────

tdp_colors <- c(
  "CTRL" = "#4575b4",
  "NON"  = "#fee090",
  "MLD"  = "#fdae61",
  "MOD"  = "#f46d43",
  "SEV"  = "#d73027"
)

# Panel A: PCA ----
p_sifig1_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = TDP_Lvl, shape = Group)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = tdp_colors, name = "TDP-43 Level") +
  scale_shape_manual(values = c("ALS" = 16, "CTRL" = 17), name = "Condition") +
  labs(
    x = sprintf("PC1 (%.1f%%)", pca_var[1]),
    y = sprintf("PC2 (%.1f%%)", pca_var[2]),
    title = "PCA of Protein Intensities"
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 13 * TEXT_SIZE),
    axis.title      = element_text(size = 11 * TEXT_SIZE),
    axis.text       = element_text(size = 10 * TEXT_SIZE),
    legend.text     = element_text(size = 10 * TEXT_SIZE),
    legend.title    = element_text(size = 10 * TEXT_SIZE),
    legend.position = "right"
  )

# Panel B: within-ALS TDP-43 trend volcano ----
pval_cutoff_tdp <- 0.05

table_tdp <- table_tdp %>%
  mutate(Significant = adj.P.Val < pval_cutoff_tdp)

n_sig <- sum(table_tdp$Significant)

p_sifig1_volcano <- ggplot(table_tdp, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = Significant), alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#d73027"),
                     labels = c("Not significant", sprintf("FDR < %.2f", pval_cutoff_tdp)),
                     name = NULL) +
  geom_hline(yintercept = -log10(max(table_tdp$P.Value[table_tdp$Significant],
                                     na.rm = TRUE)),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("DA proteins: %d", n_sig),
           hjust = 1.1, vjust = 1.5, size = 4 * TEXT_SIZE, color = "grey20") +
  labs(
    x = "Log2 Fold Change (per TDP-43 severity unit)",
    y = "-Log10 P-value",
    title = "Within-ALS TDP-43 Severity Trend"
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 13 * TEXT_SIZE),
    axis.title      = element_text(size = 11 * TEXT_SIZE),
    axis.text       = element_text(size = 10 * TEXT_SIZE),
    legend.text     = element_text(size = 10 * TEXT_SIZE),
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.9)
  )

p_sifig1 <- plot_grid(p_sifig1_pca, p_sifig1_volcano, nrow = 1, labels = c("A", "B"))
ggsave(file.path(FIG_DIR, "SIFig1_TDP43_Clustering.png"),
       p_sifig1, width = 12, height = 5, dpi = 300)


# ── SI Figure 2 | Carboxymethylation Analysis ─────────────────────────────────
# TODO: migrate CML plots from AlsMotorNeuronAnalysis once analysis code is in place.
# Likely panels: per-protein CML peptide counts, logFC distribution for CML peptides.


# ── SI Figure 3 | FlashLFQ Protein-Level DA: Diverse vs Limited PTMs ─────────
# Parallel to Figure 4 but using FlashLFQ quantification.
# Requires merged_pro_flashlfq from Supplemental/analysis.R.

# modiff_colors <- c("Shared" = "grey80", "ModDA" = "forestgreen", "NoModDA" = "darkred")

# Diverse PTM volcano (FlashLFQ) ----
# p_sifig3_volcano <- ggplot(table_protein_diverse_flashlfq, ...) + ...
# ggsave(file.path(FIG_DIR, "SIFig3_Volcano_Diverse_FlashLFQ.png"), ...)

# ModDiff volcano (FlashLFQ) ----
# p_sifig3_moddiff <- ggplot(merged_pro_flashlfq, ...) + ...
# ggsave(file.path(FIG_DIR, "SIFig3_ModDiff_FlashLFQ.png"), ...)

# logFC scatter (FlashLFQ) ----
# p_sifig3_scatter <- ggplot(merged_pro_flashlfq, ...) + ...
# ggsave(file.path(FIG_DIR, "SIFig3_LogFC_Scatter_FlashLFQ.png"), ...)

# Venn diagram (FlashLFQ) ----
# diverse_de_flashlfq <- table_protein_diverse_flashlfq$Gene[table_protein_diverse_flashlfq$DEStatus != "NotDE"]
# limited_de_flashlfq <- table_protein_limited_flashlfq$Gene[table_protein_limited_flashlfq$DEStatus != "NotDE"]
# fit_euler_fl <- euler(list("Diverse PTMs" = diverse_de_flashlfq,
#                            "Limited PTMs" = limited_de_flashlfq))
# png(file.path(FIG_DIR, "SIFig3_Venn_FlashLFQ.png"), width = 5, height = 4, units = "in", res = 300)
# print(plot(fit_euler_fl, ...))
# dev.off()
