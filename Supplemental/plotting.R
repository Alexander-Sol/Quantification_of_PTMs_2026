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
# Requires merged_pro_flashlfq from Supplemental/analysis.R and
# table_occ / n_mod_strict from MainManuscript/analysis.R + plotting.R.

modiff_colors_fl <- c("Shared" = "grey80", "ModDA" = "#ef5350", "NoModDA" = "#5c6bc0")
de_status_colors_fl <- c("NotDE" = "grey", "UP" = "#ef5350", "DOWN" = "#ef5350")

genes_highlight_fl <- c("NEFL", "TUBB", "MAP1B", "UCHL1",
                         "CKB", "CTSD", "PSAP", "ENO1",
                         "HSPD1", "HSPA8")

# n_mod_strict: per-protein count of truly modified peptides from the limpa
# analysis. Recompute here from table_occ if not already in the environment
# (table_occ is produced by MainManuscript/analysis.R).
if (!exists("n_mod_strict")) {
  artifact_only <- c("Oxidation", "Deamidation", "Deamidated", "Carbamidomethyl")
  n_mod_strict  <- table_occ %>%
    filter(
      Category %in% c("Enzymatic Mod", "Carboxymethylation") |
      (Category == "Non-Enzymatic Mod" &
         map_lgl(str_split(AllMod, ", "), ~ length(setdiff(.x, artifact_only)) > 0))
    ) %>%
    group_by(Gene) %>%
    summarise(n_mod = n(), .groups = "drop")
}

# Build per-protein data frame for the discordance plots ----
pc_fl <- merged_pro_flashlfq %>%
  select(Gene, LogFC_Diff, ModDiff, logFC_NoMod, adj.P.Val, adj.P.Val_NoMod) %>%
  inner_join(n_mod_strict, by = "Gene") %>%
  filter(n_mod > 10, n_mod < 200) %>%
  mutate(p_diff = -log10(adj.P.Val) - (-log10(adj.P.Val_NoMod)))

highlight_fl      <- merged_pro_flashlfq %>% filter(Gene %in% genes_highlight_fl)
highlight_fl_disc <- pc_fl              %>% filter(Gene %in% genes_highlight_fl)

# Shared theme helpers ----
sifig3_theme <- theme_minimal() +
  theme(
    plot.title      = element_text(size = 14 * TEXT_SIZE),
    axis.title      = element_text(size = 12 * TEXT_SIZE),
    axis.text       = element_text(size = 11 * TEXT_SIZE),
    legend.position = "none"
  )

# Diverse PTM volcano ----
p_sifig3_volcano <- ggplot(table_protein_diverse_flashlfq,
    aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = DEStatus), alpha = 0.6, size = 2) +
  scale_color_manual(values = de_status_colors_fl) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_hline(yintercept = -log10(pval_cutoff),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_point(data = filter(table_protein_diverse_flashlfq, Gene %in% genes_highlight_fl),
             size = 3.5, shape = 21, fill = "#43a047", color = "black", stroke = 0.7) +
  geom_text_repel(data = filter(table_protein_diverse_flashlfq, Gene %in% genes_highlight_fl),
                  aes(label = Gene), color = "black", size = 4.5, fontface = "bold",
                  box.padding = 0.5, min.segment.length = 0, max.overlaps = Inf) +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value",
       title = "DA Proteins: Diverse PTMs\n(FlashLFQ)") +
  sifig3_theme

ggsave(file.path(FIG_DIR, "SIFig3_Volcano_Diverse.png"), p_sifig3_volcano, width = 6, height = 5, dpi = 300)

# ModDiff volcano ----
p_sifig3_moddiff <- ggplot(merged_pro_flashlfq,
    aes(x = logFC, y = -log10(adj.P.Val), color = ModDiff)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(values = modiff_colors_fl) +
  scale_fill_manual(values = modiff_colors_fl, guide = "none") +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_hline(yintercept = -log10(pval_cutoff),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_point(data = highlight_fl, aes(fill = ModDiff),
             size = 4, shape = 21, color = "black", stroke = 0.8) +
  geom_text_repel(data = highlight_fl, aes(label = Gene),
                  color = "black", size = 4.5, fontface = "bold",
                  box.padding = 0.5, min.segment.length = 0,
                  max.overlaps = Inf, force = 3) +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value",
       title = "DA Changes with Modifications\n(FlashLFQ)") +
  sifig3_theme

ggsave(file.path(FIG_DIR, "SIFig3_ModDiff_Volcano.png"), p_sifig3_moddiff, width = 6, height = 5, dpi = 300)

# logFC scatter: diverse vs limited ----
p_sifig3_scatter <- ggplot(merged_pro_flashlfq,
    aes(x = logFC, y = logFC_NoMod, color = ModDiff)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(values = modiff_colors_fl) +
  scale_fill_manual(values = modiff_colors_fl, guide = "none") +
  geom_hline(yintercept = 0, color = "grey20", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey20", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point(data = highlight_fl, aes(fill = ModDiff),
             size = 4, shape = 21, color = "black", stroke = 0.8) +
  geom_text_repel(data = highlight_fl, aes(label = Gene),
                  color = "black", size = 4.5, fontface = "bold",
                  box.padding = 1.5, point.padding = 0.5,
                  min.segment.length = 0, max.overlaps = Inf, force = 3) +
  labs(x = "Log2FC (Diverse PTMs)", y = "Log2FC (Limited PTMs)",
       title = "Protein Log2FC: Diverse vs Limited\n(FlashLFQ)") +
  sifig3_theme

ggsave(file.path(FIG_DIR, "SIFig3_LogFC_Scatter.png"), p_sifig3_scatter, width = 6, height = 5, dpi = 300)

# PctMod vs Discordance ----
# n_mod from limpa peptide analysis; discordance from FlashLFQ limma comparison.
p_sifig3_pctmod <- ggplot(pc_fl, aes(x = LogFC_Diff, y = p_diff, color = n_mod)) +
  geom_point(alpha = 0.8, size = 2.5) +
  scale_color_viridis_c(trans = "log10", name = "Modified\nPeptides",
                        limits = c(NA, 100), oob = scales::squish,
                        breaks = c(10, 30, 100), labels = c("10", "30", "≥100")) +
  geom_hline(yintercept = 0, color = "grey20", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey20", linewidth = 0.4) +
  geom_point(data = highlight_fl_disc,
             size = 4, shape = 21, fill = NA, color = "black", stroke = 0.8) +
  geom_text_repel(data = highlight_fl_disc, aes(label = Gene),
                  color = "black", size = 4.5, fontface = "bold",
                  box.padding = 0.5, min.segment.length = 0,
                  max.overlaps = Inf, force = 3) +
  labs(x = "Log2FC Difference (Diverse − Limited)",
       y = "Δ (−log10 Adj. P-value)",
       title = "Modified Peptide Count vs. DA Discordance\n(FlashLFQ)") +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 14 * TEXT_SIZE),
    axis.title      = element_text(size = 12 * TEXT_SIZE),
    axis.text       = element_text(size = 11 * TEXT_SIZE),
    legend.position = "right",
    legend.text     = element_text(size = 10 * TEXT_SIZE),
    legend.title    = element_text(size = 11 * TEXT_SIZE)
  )

ggsave(file.path(FIG_DIR, "SIFig3_PctMod_vs_Discordance.png"), p_sifig3_pctmod, width = 6.5, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "SIFig3_PctMod_vs_Discordance.svg"), p_sifig3_pctmod, width = 6.5, height = 5)

# Limited PTM volcano colored by ModDiff ----
highlight_fl_lim <- merged_pro_flashlfq %>% filter(Gene %in% genes_highlight_fl)

p_sifig3_lim_volcano <- ggplot(merged_pro_flashlfq,
    aes(x = logFC_NoMod, y = -log10(adj.P.Val_NoMod), color = ModDiff)) +
  geom_point(alpha = 0.8, size = 2.5) +
  scale_color_manual(values = modiff_colors_fl) +
  scale_fill_manual(values = modiff_colors_fl, guide = "none") +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_hline(yintercept = -log10(pval_cutoff),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_point(data = highlight_fl_lim,
             aes(x = logFC_NoMod, y = -log10(adj.P.Val_NoMod), fill = ModDiff),
             size = 4, shape = 21, color = "black", stroke = 0.8) +
  geom_text_repel(data = highlight_fl_lim,
                  aes(x = logFC_NoMod, y = -log10(adj.P.Val_NoMod), label = Gene),
                  color = "black", size = 4.5, fontface = "bold",
                  box.padding = 0.5, min.segment.length = 0,
                  max.overlaps = Inf, force = 3) +
  labs(x = "Log2 Fold Change (Limited PTMs)",
       y = "-Log10 Adjusted P-value (Limited PTMs)",
       title = "DA Proteins: Limited PTMs\n(FlashLFQ)") +
  sifig3_theme

ggsave(file.path(FIG_DIR, "SIFig3_LimitedVolcano.png"), p_sifig3_lim_volcano, width = 6, height = 5, dpi = 300)

# Venn diagram ----
fit_euler_fl <- euler(list(
  "Diverse PTMs" = flashlfq_diverse_de,
  "Limited PTMs" = flashlfq_limited_de
))

p_sifig3_venn <- plot(fit_euler_fl,
                      quantities = list(cex = 1.3),
                      labels     = list(cex = 1.1, box = list(col = NA, fill = alpha("white", 0.6))),
                      fills      = list(fill = c("#ef5350", "#5c6bc0"), alpha = 0.45),
                      edges      = list(col = "grey20", lwd = 1.5),
                      main       = "DE Proteins: Diverse vs. Limited PTMs\n(FlashLFQ)")

png(file.path(FIG_DIR, "SIFig3_Venn.png"), width = 5, height = 4, units = "in", res = 300)
print(p_sifig3_venn)
dev.off()

# Combined panel ----
p_sifig3_combined <- plot_grid(
  p_sifig3_volcano, p_sifig3_moddiff, p_sifig3_scatter,
  p_sifig3_lim_volcano, p_sifig3_pctmod,
  nrow = 2, labels = c("A", "B", "C", "D", "E")
)
ggsave(file.path(FIG_DIR, "SIFig3_Combined.png"), p_sifig3_combined,
       width = 18, height = 10, dpi = 300)
