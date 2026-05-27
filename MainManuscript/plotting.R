# plotting.R — Main Manuscript
# All publication figures. Requires analysis.R to have been sourced first.
# Figures are saved to MainManuscript/Figures/.

library(tidyverse)
library(cowplot)
library(scales)
library(ggrepel)
library(eulerr)

source("MainManuscript/CustomScripts.R")

FIG_DIR  <- "MainManuscript/Figures"
TEXT_SIZE <- 1.1   # global scale factor for all text elements


# ── QC | Intensity Distribution by Condition ──────────────────────────────────
# Density plots of log2 peptide and protein intensities, one curve per sample,
# coloured by ALS vs CTRL. Run before any filtering to confirm no systematic
# batch shift between conditions.

condition_colors <- c("ALS" = "#d73027", "CTRL" = "#4575b4")

# Peptide-level intensities ----
pep_long <- as.data.frame(y.peptide$E) %>%
  rownames_to_column("PeptideSequence") %>%
  pivot_longer(-PeptideSequence, names_to = "Sample", values_to = "Log2Int") %>%
  filter(!is.na(Log2Int)) %>%
  left_join(y.peptide$targets[, c("Sample", "Group", "Donor")], by = "Sample")

pep_mean <- pep_long %>%
  group_by(PeptideSequence, Group) %>%
  summarise(Log2Int = mean(Log2Int), .groups = "drop")

p_qc_peptide <- ggplot(pep_mean, aes(x = Log2Int, color = Group)) +
  geom_density(linewidth = 0.4, alpha = 0.7) +
  scale_color_manual(values = condition_colors) +
  labs(x = "Mean Log2 Peptide Intensity", y = "Density",
       title = "QC: Peptide Intensity Distributions",
       color = "Condition") +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 13 * TEXT_SIZE),
    axis.title      = element_text(size = 11 * TEXT_SIZE),
    axis.text       = element_text(size = 10 * TEXT_SIZE),
    legend.text     = element_text(size = 10 * TEXT_SIZE),
    legend.position = "right"
  )

# Protein-level intensities ----
pro_long <- as.data.frame(y.protein$E) %>%
  rownames_to_column("ProteinGroup") %>%
  pivot_longer(-ProteinGroup, names_to = "Sample", values_to = "Log2Int") %>%
  filter(!is.na(Log2Int)) %>%
  left_join(y.protein$targets[, c("Sample", "Group", "Donor")], by = "Sample")

p_qc_protein <- ggplot(pro_long, aes(x = Log2Int, color = Group)) +
  geom_density(linewidth = 0.4, alpha = 0.7) +
  scale_color_manual(values = condition_colors) +
  labs(x = "Log2 Protein Intensity", y = "Density",
       title = "QC: Protein Intensity Distributions",
       color = "Condition") +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 13 * TEXT_SIZE),
    axis.title      = element_text(size = 11 * TEXT_SIZE),
    axis.text       = element_text(size = 10 * TEXT_SIZE),
    legend.text     = element_text(size = 10 * TEXT_SIZE),
    legend.position = "right"
  )

p_qc_combined <- plot_grid(p_qc_peptide, p_qc_protein, nrow = 1, labels = c("A", "B"))
ggsave(file.path(FIG_DIR, "QC_IntensityDistributions.png"),
       p_qc_combined, width = 12, height = 5, dpi = 300)


# Unfiltered peptide intensity distribution (raw QuantifiedPeptides.tsv) ----
df_raw          <- read.csv(PATH_DIVERSE, sep = "\t", row.names = NULL)
intensity_cols  <- grep("^Intensity_", colnames(df_raw), value = TRUE)
expr_raw        <- GetExpressionMatrix(df_raw, intensity_cols)
targets_raw     <- GetMetaDataTdpStratified(expr_raw)

pep_long_raw <- as.data.frame(expr_raw) %>%
  rownames_to_column("PeptideSequence") %>%
  pivot_longer(-PeptideSequence, names_to = "Sample", values_to = "Log2Int") %>%
  filter(!is.na(Log2Int)) %>%
  left_join(targets_raw[, c("Sample", "Group")], by = "Sample")

pep_mean_raw <- pep_long_raw %>%
  group_by(PeptideSequence, Group) %>%
  summarise(Log2Int = mean(Log2Int), .groups = "drop")

p_qc_unfiltered <- ggplot(pep_mean_raw, aes(x = Log2Int, color = Group)) +
  geom_density(linewidth = 0.4, alpha = 0.7) +
  scale_color_manual(values = condition_colors) +
  labs(x = "Mean Log2 Peptide Intensity", y = "Density",
       title = "QC: Unfiltered Peptide Intensity Distributions",
       color = "Condition") +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 13 * TEXT_SIZE),
    axis.title      = element_text(size = 11 * TEXT_SIZE),
    axis.text       = element_text(size = 10 * TEXT_SIZE),
    legend.text     = element_text(size = 10 * TEXT_SIZE),
    legend.position = "right"
  )

p_qc_unfiltered


# ── Figure 1 | PSM & Peptide Discovery: Diverse vs Limited PTMs ──────────────

search_colors <- c("Search 2: Diverse PTMs" = "#78CFDE",
                   "Search 1: Limited PTMs"  = "#44648E")

# Per-cell PSM counts ----
p_fig1_psm_per_cell <- ggplot(psm_counts_combined, aes(x = File, y = count, fill = Type)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.7) +
  scale_fill_manual(values = search_colors, guide = guide_legend(reverse = TRUE)) +
  labs(title = "PSMs per Cell: Diverse vs. Limited PTMs",
       x = "Sample", y = "PSM Count") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y  = element_text(size = 11 * TEXT_SIZE),
    axis.title   = element_text(size = 12 * TEXT_SIZE),
    plot.title   = element_text(size = 13 * TEXT_SIZE),
    legend.position        = c(0.02, 0.98),
    legend.justification   = c(0, 1),
    legend.background      = element_rect(fill = alpha("white", 0.6), color = NA),
    legend.title           = element_blank(),
    legend.text            = element_text(size = 11 * TEXT_SIZE),
    legend.key.size        = unit(1.4, "lines")
  )

# Total PSM inset bar ----
p_fig1_psm_total <- ggplot(psm_totals, aes(x = 1, y = Total, fill = Type)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.7, width = 0.8) +
  geom_text(aes(label = scales::comma(Total)), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = search_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                     labels = scales::comma) +
  coord_cartesian(clip = "off") +
  labs(title = "Total\nPSMs", x = "", y = "Number of PSMs") +
  theme_minimal() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 10 * TEXT_SIZE),
    axis.title   = element_text(size = 11 * TEXT_SIZE),
    plot.title   = element_text(size = 12 * TEXT_SIZE, hjust = 0.5),
    legend.position = "none"
  )

# Per-cell peptide counts ----
p_fig1_pep_per_cell <- ggplot(pep_counts_combined, aes(x = File, y = count, fill = Type)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.7) +
  scale_fill_manual(values = search_colors, guide = guide_legend(reverse = TRUE)) +
  labs(title = "Peptidoforms per Cell: Diverse vs. Limited PTMs",
       x = "Sample", y = "Unique Peptidoform Count") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
    axis.text.y  = element_text(size = 11 * TEXT_SIZE),
    axis.title   = element_text(size = 12 * TEXT_SIZE),
    plot.title   = element_text(size = 13 * TEXT_SIZE),
    legend.position        = c(0.02, 0.98),
    legend.justification   = c(0, 1),
    legend.background      = element_rect(fill = alpha("white", 0.6), color = NA),
    legend.title           = element_blank(),
    legend.text            = element_text(size = 11 * TEXT_SIZE),
    legend.key.size        = unit(1.4, "lines")
  )

# Total peptide inset bar ----
p_fig1_pep_total <- ggplot(pep_totals, aes(x = 1, y = Total, fill = Type)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.7, width = 0.6) +
  geom_text(aes(label = scales::comma(Total)), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = search_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                     labels = scales::comma) +
  coord_cartesian(clip = "off") +
  labs(title = "Total\nUnique\nPeptidoforms", x = "", y = "Unique Peptidoforms") +
  theme_minimal() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 10 * TEXT_SIZE),
    axis.title   = element_text(size = 11 * TEXT_SIZE),
    plot.title   = element_text(size = 12 * TEXT_SIZE, hjust = 0.5,
                               margin = margin(l = -40, unit = "pt")),
    legend.position = "none"
  )
ggsave(file.path(FIG_DIR, "Fig1_Pep_Total.png"),      p_fig1_pep_total,     width = 1.5, height = 5.5,   dpi = 300)


# Biological modifications detected ----
p_fig1_bio_mods <- ggplot(bio_mods_df, aes(x = all_mods_renamed, y = Freq)) +
  geom_bar(stat = "identity", fill = "#d73027", alpha = 0.9,
           color = "gray20", linewidth = 0.8) +
  coord_flip() +
  scale_y_log10(labels = scales::comma) +
  labs(title = "Enzymatic Modifications Detected",
       x = "", y = "Number of PTMs Observed\nAcross PSMs (log scale)") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 11 * TEXT_SIZE),
    axis.text.y  = element_text(size = 12 * TEXT_SIZE,
                                margin = margin(l = 5, unit = "pt")),
    axis.title   = element_text(size = 12 * TEXT_SIZE),
    plot.title   = element_text(size = 13 * TEXT_SIZE, hjust = 0.5)
  )

# Other modifications detected ----
p_fig1_other_mods <- ggplot(other_mods_df, aes(x = all_mods_renamed, y = Freq)) +
  geom_bar(stat = "identity", fill = "#f46d43", alpha = 0.9,
           color = "gray20", linewidth = 0.8) +
  coord_flip() +
  scale_y_log10(labels = scales::comma) +
  labs(title = "Non-Enzymatic Modifications Detected",
       x = "", y = "Number of PTMs Observed\nAcross PSMs (log scale)") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 10 * TEXT_SIZE),
    axis.text.y  = element_text(size = 12 * TEXT_SIZE),
    axis.title   = element_text(size = 12 * TEXT_SIZE),
    plot.title   = element_text(size = 13 * TEXT_SIZE, hjust = 0.5)
  )

# Save individual panels (for flexible layout in figure assembly)
ggsave(file.path(FIG_DIR, "Fig1_PSM_PerCell.png"),    p_fig1_psm_per_cell,  width = 7,   height = 5.5,   dpi = 300)
ggsave(file.path(FIG_DIR, "Fig1_PSM_Total.png"),      p_fig1_psm_total,     width = 1.75, height = 5.5,   dpi = 300)
ggsave(file.path(FIG_DIR, "Fig1_Pep_PerCell.png"),    p_fig1_pep_per_cell,  width = 7,   height = 5.5,   dpi = 300)
ggsave(file.path(FIG_DIR, "Fig1_Pep_Total.png"),      p_fig1_pep_total,     width = 1.5, height = 5.5,   dpi = 300)
ggsave(file.path(FIG_DIR, "Fig1_BioMods.png"),        p_fig1_bio_mods,      width = 5,   height = 5.5,   dpi = 300)
ggsave(file.path(FIG_DIR, "Fig1_OtherMods.png"),      p_fig1_other_mods,    width = 5,   height = 5.5,   dpi = 300)

# Combined layout (top: PSM+total | pep+total; bottom: bio mods | other mods)
top_row <- plot_grid(p_fig1_psm_per_cell,  p_fig1_psm_total,
                     p_fig1_pep_per_cell,  p_fig1_pep_total,
                     nrow = 1, rel_widths = c(7, 2.0, 7, 1.5),
                     labels = c("A", "", "B", ""))
bot_row <- plot_grid(p_fig1_bio_mods, p_fig1_other_mods,
                     nrow = 1, labels = c("C", "D"))
p_fig1  <- plot_grid(top_row, bot_row, nrow = 2)
ggsave(file.path(FIG_DIR, "Fig1_Combined.png"), p_fig1, width = 14, height = 10, dpi = 300)


# ── Figure 2 | Peptide Category Composition by Stringency ────────────────────
category_levels_ordered <- c("Enzymatic Mod", "Carboxymethylation",
                              "Non-Enzymatic Mod", "Carbamidomethylation", "Unmodified")

stringency_defs <- list(
  list(label = "All Peptides",  lfc = 0,   sig = FALSE),
  list(label = "Significant",   lfc = 0,   sig = TRUE),
  list(label = "|logFC| ≥ 1.5", lfc = 1.5, sig = TRUE)
)

apply_stringency <- function(lfc, sig) {
  if (!sig) return(base_filtered)
  base_filtered %>% filter(adj.P.Val <= pval_cutoff, abs(logFC) >= lfc)
}

bar_df <- map_dfr(stringency_defs, function(lvl) {
  filt    <- apply_stringency(lvl$lfc, lvl$sig)
  n_total <- nrow(filt)
  filt %>%
    group_by(Category) %>%
    summarise(Count = n(), .groups = "drop") %>%
    complete(Category = factor(category_levels_ordered, levels = category_levels_ordered),
             fill = list(Count = 0)) %>%
    mutate(percent = 100 * Count / n_total,
           Group   = lvl$label,
           lfc     = lvl$lfc)
})

bar_df$Category <- factor(bar_df$Category, levels = category_levels_ordered)
bar_df$Group    <- factor(bar_df$Group,
                          levels = map_chr(stringency_defs, ~ .x$label))

small_threshold <- 5  # minimum percent for label display to avoid clutter; adjust as needed

bar_totals <- bar_df %>%
  group_by(Group) %>%
  summarise(n_total = sum(Count), .groups = "drop")

p_fig2_bar <- ggplot(bar_df, aes(x = Group, y = percent, fill = Category)) +
  geom_col(position = "stack", alpha = 0.8, color = "black", linewidth = 0.4) +
  geom_text(
    data = filter(bar_df, percent >= small_threshold),
    aes(label = sprintf("%d%%", round(percent))),
    position = position_stack(vjust = 0.5),
    size = 4, color = "black"
  ) +
  geom_text(
    data = bar_totals,
    aes(x = Group, y = 100, label = paste0("n=", scales::comma(n_total))),
    vjust = -0.4, size = 4, color = "gray20", inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Enzymatic Mod"        = "#d73027",
      "Carboxymethylation"   = "#f46d43",
      "Non-Enzymatic Mod"    = "#fdae61",
      "Carbamidomethylation" = "#fee090",
      "Unmodified"           = "#4575b4"
    ),
    breaks = category_levels_ordered
  ) +
  labs(y = "Percent of Peptides") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12 * TEXT_SIZE),
    axis.text.x  = element_text(size = 10 * TEXT_SIZE, angle = 30, hjust = 1),
    axis.text.y  = element_text(size = 11 * TEXT_SIZE),
    legend.position = "bottom",
    legend.text  = element_text(size = 11 * TEXT_SIZE),
    legend.title = element_blank()
  )

ggsave(file.path(FIG_DIR, "Fig2_CategoryBar.png"), p_fig2_bar, width = 3.5, height = 5.5, dpi = 300)


# ── Figure 2 | Volcano + Modification Breakdown ───────────────────────────────

mod_colors <- c(
  "Unmodified"           = "#4575b4",
  "Carbamidomethylation" = "#fee090",
  "Non-Enzymatic Mod"    = "#fdae61",
  "Carboxymethylation"   = "#f46d43",
  "Enzymatic Mod"        = "#d73027"
)

volcano_df <- base_filtered %>%
  mutate(
    InTop1000 = adj.P.Val <= pval_cutoff & abs(logFC) >= logFC_cutoff,
    Highlight = case_when(
      InTop1000 & Gene %in% c("NEFH", "NEFM", "NEFL")                        ~ "Neurofilaments",
      InTop1000 & Gene == "GFAP"                                               ~ "GFAP",
      InTop1000 & Gene %in% c("EEF1G", "EEF1B2", "EEF1D", "EEF1A1", "EEF1A2", "EEF2") ~ "Elongation Factors",
      TRUE ~ NA_character_
    ),
    y_display = -log10(adj.P.Val)
  )

highlight_df <- volcano_df %>% filter(!is.na(Highlight)) %>%
  mutate(Highlight = factor(Highlight, levels = c("GFAP", "Elongation Factors", "Neurofilaments")))

sig_y <- -log10(pval_cutoff)

upper_compress_trans <- scales::trans_new(
  name      = "upper_compress",
  transform = function(x) ifelse(x <= 5, x, 5 + log1p(x - 5) * 1),
  inverse   = function(x) ifelse(x <= 5, x, 5 + expm1((x - 5) / 1))
)

p_fig3_volcano <- ggplot(volcano_df, aes(x = logFC, y = y_display)) +
  geom_point(data = filter(volcano_df, !InTop1000, is.na(Highlight)),
             color = "grey75", alpha = 0.4, size = 1.5, shape = 19) +
  geom_point(data = filter(volcano_df, InTop1000, is.na(Highlight)),
             aes(color = Category), alpha = 0.8, size = 1.8, shape = 19) +
  geom_point(data = filter(highlight_df, Highlight == "GFAP"),
             aes(fill = Category, shape = Highlight),
             size = 2.3, alpha = 0.9, stroke = 0.6, color = "black") +
  geom_point(data = filter(highlight_df, Highlight == "Elongation Factors"),
             aes(fill = Category, shape = Highlight),
             size = 2.9, alpha = 0.9, stroke = 0.6, color = "black") +
  geom_point(data = filter(highlight_df, Highlight == "Neurofilaments"),
             aes(fill = Category, shape = Highlight),
             size = 3.0, alpha = 0.9, stroke = 0.6, color = "black") +
  geom_hline(yintercept = sig_y,
             color = "grey30", linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff),
             color = "grey30", linewidth = 0.5, linetype = "dashed") +
  scale_color_manual(values = mod_colors, breaks = names(mod_colors)) +
  scale_fill_manual(values = mod_colors, guide = "none") +
  scale_shape_manual(values = c("GFAP" = 24, "Elongation Factors" = 22, "Neurofilaments" = 23), name = NULL) +
  scale_y_continuous(trans   = upper_compress_trans,
                     limits  = c(0, 12.5),
                     breaks  = c(0, 1, 2, 3, 4, 5, 7, 10, 12)) +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value",
       title = "Differing Peptidoform Fractions") +
  theme_minimal() +
  theme(
    plot.title  = element_text(size = 14 * TEXT_SIZE),
    axis.title  = element_text(size = 12 * TEXT_SIZE),
    axis.text   = element_text(size = 11 * TEXT_SIZE),
    legend.position      = c(0.46, 0.98),
    legend.justification = c(0.5, 1),
    legend.background    = element_rect(fill = alpha("white", 0.6), color = NA),
    legend.text          = element_text(size = 11 * TEXT_SIZE),
    legend.title         = element_blank(),
    legend.box           = "vertical"
  ) +
  guides(color = guide_legend(order = 1, reverse = TRUE, override.aes = list(size = 3.5)),
         shape = guide_legend(order = 2, override.aes = list(size = 3.5)))
p_fig3_volcano

ggsave(file.path(FIG_DIR, "Fig3_Volcano.png"), p_fig3_volcano, width = 8.5, height = 6, dpi = 300)

# Non-enzymatic mod bar chart ----
other_mod_df <- table(
  GetModsTextWReplacement(de_splice$PeptideSequence[de_splice$Category == "Non-Enzymatic Mod"]) %>%
    unlist() %>% str_split(., ", ") %>% unlist()
) %>% as.data.frame() %>%
  filter(!grepl("Carboxymeth|Carbamidomethy", Var1)) %>%
  arrange(desc(Freq)) %>%
  mutate(mod = factor(Var1, levels = Var1),
         mod = fct_recode(mod, "dialkylation*" = "Ubiquitination"))

other_mod_df <- other_mod_df[1:12, ]

p_fig3_nonenz_bar <- ggplot(other_mod_df, aes(x = reorder(mod, Freq), y = Freq)) +
  geom_bar(stat = "identity", color = "black", fill = "#fdae61", alpha = 0.9, linewidth = 0.5) +
  coord_flip() +
  labs(title = "Non-Enzymatic Modifications in\nDifferentially Modified Peptides",
       x = "", y = "Instances of Modification Observed") +
  theme_minimal() +
  theme(
    plot.title   = element_text(size = 14 * TEXT_SIZE, hjust = 0),
    axis.title.x = element_text(size = 12 * TEXT_SIZE),
    axis.text    = element_text(size = 11 * TEXT_SIZE)
  )

ggsave(file.path(FIG_DIR, "Fig3_NonEnzBar.png"), p_fig3_nonenz_bar, width = 6, height = 4.5, dpi = 300)

# Enzymatic mod bar chart ----
uninteresting_mods <- c("Carbamidomethyl", "Carbamyl", "Carboxymethylation",
                        "Formylation", "Sodium", "Potassium", "Magnesium",
                        "Calcium", "Deamidation", "Deamidated", "Fe[III", "Iron",
                        "Ammonia", "Water", "Water Loss", "Carboxylation", "Oxidation",
                        "Carboxymethyllysine", "Nitrosylation")

bio_mod_df <- table(
  GetModsTextWReplacement(de_splice$PeptideSequence[de_splice$Category == "Enzymatic Mod"]) %>%
    unlist() %>% str_split(., ", ") %>% unlist() %>%
    .[!(. %in% uninteresting_mods)]
) %>% as.data.frame() %>% arrange(desc(Freq)) %>%
  mutate(mod = factor(Var1, levels = Var1))

bio_mod_df <- bio_mod_df[1:12,] # Keep the top twelve most common enzymatic mods for clarity; adjust as needed

p_fig3_enz_bar <- ggplot(bio_mod_df, aes(x = reorder(mod, Freq), y = Freq)) +
  geom_bar(stat = "identity", color = "black", fill = "#d73027", alpha = 0.9, linewidth = 0.5) +
  coord_flip() +
  labs(title = "Enzymatic Modifications in\nDifferentially Modified Peptides",
       x = "", y = "Instances of Modification Observed") +
  theme_minimal() +
  theme(
    plot.title   = element_text(size = 14 * TEXT_SIZE, hjust = 0),
    axis.title.x = element_text(size = 12 * TEXT_SIZE),
    axis.text    = element_text(size = 11 * TEXT_SIZE)
  )

ggsave(file.path(FIG_DIR, "Fig3_EnzBar.png"), p_fig3_enz_bar, width = 6, height = 4.5, dpi = 300)

# FIGFURE ---- Protein-Level Breakdown of Modifications ----
# Enzymatic mod by protein (top 10) ----
bio_mods_of_interest <- c("Citrullination", "Acetylation", "Phosphorylation",
                           "Methylation", "Dimethylation", "Trimethylation",
                           "Hydroxyproline", "Glutarylation", "Butyrylation")

gene_table_bio_mod <- de_splice$Gene[de_splice$Category == "Enzymatic Mod"] %>%
  table() %>% as.data.frame() %>% arrange(desc(Freq))
top10_genes <- gene_table_bio_mod$.[1:10]

bio_expanded <- de_splice %>%
  filter(Category == "Enzymatic Mod", Gene %in% top10_genes) %>%
  select(Gene, AllMod) %>%
  mutate(Mod = strsplit(AllMod, ", ")) %>%
  unnest(Mod) %>%
  filter(Mod %in% bio_mods_of_interest) %>%
  mutate(Mod = factor(Mod, levels = rev(bio_mods_of_interest)))

bio_gene_counts <- bio_expanded %>% group_by(Gene, Mod) %>% summarise(Count = n(), .groups = "drop")
gene_order_bio  <- bio_gene_counts %>% group_by(Gene) %>%
  summarise(Total = sum(Count), .groups = "drop") %>% arrange(desc(Total)) %>% pull(Gene)
bio_gene_counts$Gene <- factor(bio_gene_counts$Gene, levels = rev(gene_order_bio))

bio_protein_colors <- c(
  "Citrullination"  = "darkblue",  "Acetylation"     = "forestgreen",
  "Methylation"     = "darkorchid2", "Dimethylation"  = "darkorchid3",
  "Trimethylation"  = "darkorchid4", "Phosphorylation" = "cyan3",
  "Hydroxyproline"  = "deeppink3",  "Glutarylation"   = "sienna",
  "Butyrylation"    = "sienna2"
)

p_fig3_enz_protein <- ggplot(arrange(bio_gene_counts, Mod), aes(x = Gene, y = Count, fill = Mod)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "gray20", linewidth = 0.3) +
  scale_fill_manual(values = bio_protein_colors, limits = bio_mods_of_interest) +
  coord_flip() +
  labs(x = "", y = "Number of Mods on Differentially Modified Peptides", fill = "",
       title = "Proteins with Differentially Modified\nPeptides: Enzymatic Modifications") +
  scale_y_continuous(breaks = breaks_width(10), minor_breaks = breaks_width(5)) +
  theme_minimal() +
  theme(
    plot.title         = element_text(size = 14 * TEXT_SIZE),
    axis.title.x       = element_text(size = 12 * TEXT_SIZE),
    axis.text          = element_text(size = 11 * TEXT_SIZE),
    legend.position    = "right",
    legend.justification = "left",
    legend.margin      = margin(l = -40),
    legend.text        = element_text(size = 11 * TEXT_SIZE),
    legend.title       = element_blank(),
    panel.grid.minor   = element_line(color = "grey90", linewidth = 0.25)
  )

ggsave(file.path(FIG_DIR, "Fig3_EnzProteinBar.png"), p_fig3_enz_protein, width = 6.5, height = 4.5, dpi = 300)

# Carboxymethylation by protein ----
cml_mods_of_interest <- c("Carboxymethylation", "Carboxymethyllysine")

cml_expanded <- de_splice %>%
  filter(Category == "Carboxymethylation") %>%
  select(Gene, AllMod) %>%
  mutate(Mod = strsplit(AllMod, ", ")) %>%
  unnest(Mod) %>%
  filter(Mod %in% cml_mods_of_interest) %>%
  mutate(Mod = factor(Mod, levels = cml_mods_of_interest))

cml_gene_counts <- cml_expanded %>%
  group_by(Gene, Mod) %>%
  summarise(Count = n(), .groups = "drop")

gene_order_cml <- cml_gene_counts %>%
  group_by(Gene) %>%
  summarise(Total = sum(Count), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Gene)

top12_cml_genes <- head(gene_order_cml, 12)
cml_gene_counts <- cml_gene_counts %>% filter(Gene %in% top12_cml_genes)
cml_gene_counts$Gene <- factor(cml_gene_counts$Gene, levels = rev(top12_cml_genes))

cml_colors <- c(
  "Carboxymethylation"  = "#fdae61",  # Non-Enzymatic Mod color
  "Carboxymethyllysine" = "#f46d43"   # Carboxymethylation color
)

p_fig3_cml_protein <- ggplot(arrange(cml_gene_counts, Mod),
                              aes(x = Gene, y = Count, fill = Mod)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "gray20", linewidth = 0.3) +
  scale_fill_manual(values = cml_colors, limits = rev(cml_mods_of_interest),
                    labels = c("Carboxymethylation"  = "N-Terminal Carboxymethylation",
                               "Carboxymethyllysine" = "Carboxymethyllysine")) +
  coord_flip() +
  labs(x = "", y = "Number of Mods on Differentially Modified Peptides", fill = "",
       title = "Proteins with Differentially Modified\nPeptides: Carboxymethylation") +
  scale_y_continuous(breaks = breaks_width(10), minor_breaks = breaks_width(5)) +
  theme_minimal() +
  theme(
    plot.title         = element_text(size = 14 * TEXT_SIZE),
    axis.title.x       = element_text(size = 12 * TEXT_SIZE),
    axis.text          = element_text(size = 11 * TEXT_SIZE),
    legend.position    = "bottom",
    legend.justification = "left",
    legend.margin      = margin(l = -20),
    legend.text        = element_text(size = 11 * TEXT_SIZE),
    legend.title       = element_blank(),
    panel.grid.minor   = element_line(color = "grey90", linewidth = 0.25)
  )

ggsave(file.path(FIG_DIR, "Fig3_CMLProteinBar.png"), p_fig3_cml_protein, width = 5.5, height = 4.5, dpi = 300)


# ── Figure 4 | Protein-Level DA: Diverse vs Limited PTMs ─────────────────────
# NOTE: Redo with fasta-database search results (task #27).

modiff_colors        <- c("Shared" = "grey80", "ModDA" = "#ef5350", "NoModDA" = "#5c6bc0")
de_status_colors     <- c("NotDE" = "grey", "UP" = "#ef5350", "DOWN" = "#ef5350")
#highlight_fill_pro   <- "#43a047"  # green

# PTM writers/erasers that are significantly DE in the diverse-PTM search
genes_of_interest_pro <- c("PADI2", "SIRT2", "PRKCA", "CSNK2A1", "STUB1")
# Proteins whose DA status differs between diverse and limited PTM analyses,
# selected for biological relevance and high proportions of diverse-only modified peptides
genes_of_interest_comp <- c(
  # NoModDA: stable modified peptidoforms dilute apparent decrease in ALS
  "NEFL", "TUBB", "MAP1B", "UCHL1",
  # ModDA: modified peptidoforms carry the differential signal
  "CKB", "CTSD", "PSAP", "ENO1",
  # Mitochondrial ATP synthase (NoModDA — representative subunits)
  "USP10", "STUB1"
)

# ModDiff volcano ----
highlight_comp <- merged_pro %>% filter(Gene %in% genes_of_interest_comp)

p_fig4_moddiff <- ggplot(merged_pro, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = ModDiff), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = modiff_colors) +
  scale_fill_manual(values = modiff_colors, guide = "none") +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_point(data = highlight_comp, aes(fill = ModDiff),
             size = 4, shape = 21, color = "black", stroke = 0.8) +
  geom_text_repel(data = highlight_comp, aes(label = Gene),
                  color = "black", size = 5, fontface = "bold",
                  box.padding = 0.5, min.segment.length = 0,
                  max.segment.length = Inf, max.overlaps = Inf,
                  force = 3) +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value",
       title = "Changes in Differential Abundance when\nModifications are Considered") +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 14 * TEXT_SIZE),
    axis.title      = element_text(size = 12 * TEXT_SIZE),
    axis.text       = element_text(size = 11 * TEXT_SIZE),
    legend.position = "none"
  )

ggsave(file.path(FIG_DIR, "Fig4_ModDiff_Volcano.png"), p_fig4_moddiff, width = 6, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_ModDiff_Volcano.svg"), p_fig4_moddiff, width = 6, height = 5)

# logFC scatter ----
p_fig4_scatter <- ggplot(merged_pro, aes(x = logFC, y = logFC_NoMod)) +
  geom_point(aes(color = ModDiff), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = modiff_colors) +
  scale_fill_manual(values = modiff_colors, guide = "none") +
  geom_hline(yintercept = 0, color = "grey20", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "grey20", linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point(data = highlight_comp, aes(fill = ModDiff),
             size = 4, shape = 21, color = "black", stroke = 0.8) +
  geom_text_repel(data = highlight_comp, aes(label = Gene),
                  color = "black", size = 5, fontface = "bold",
                  box.padding = 1.5, point.padding = 0.5,
                  min.segment.length = 0, segment.color = "black",
                  max.segment.length = Inf, max.overlaps = Inf,
                  force = 3) +
  labs(x = "Log2 Fold Change (Diverse PTMs)",
       y = "Log2 Fold Change (Limited PTMs)",
       title = "Comparison of Protein Log2FC With\nand Without Modified Peptides") +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 14 * TEXT_SIZE),
    axis.title      = element_text(size = 12 * TEXT_SIZE),
    axis.text       = element_text(size = 11 * TEXT_SIZE),
    legend.position = "none"
  )

ggsave(file.path(FIG_DIR, "Fig4_LogFC_Scatter.png"), p_fig4_scatter, width = 6, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_LogFC_Scatter.svg"), p_fig4_scatter, width = 6, height = 5)

# Delta logFC vs delta -log10(p), colored by modified peptide count (log scale) ----
highlight_pct_mod <- pc_for_plot %>%
  filter(Gene %in% c(genes_of_interest_comp, "HSPD1", "HSPA8"),
         !Gene %in% c("ATP5F1A", "ATP5F1B"))

p_fig4_pct_mod_scatter <- ggplot(pc_for_plot, aes(x = LogFC_Diff, y = p_diff, color = n_mod)) +
  geom_point(alpha = 0.8, size = 2.5) +
  scale_color_viridis_c(trans = "log10", name = "Modified\nPeptides",
                        limits = c(NA, 100), oob = scales::squish,
                        breaks = c(10, 30, 100),
                        labels = c("10", "30", "≥100")) +
  geom_hline(yintercept = 0, color = "grey20", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey20", linewidth = 0.4) +
  geom_point(data = highlight_pct_mod,
             size = 4, shape = 21, fill = NA, color = "black", stroke = 0.8) +
  geom_text_repel(data = highlight_pct_mod, aes(label = Gene),
                  color = "black", size = 5, fontface = "bold",
                  box.padding = 0.5, min.segment.length = 0,
                  max.overlaps = Inf, force = 3) +
  labs(x = "Log2 FC Difference (Diverse − Limited)",
       y = "-Log10 Adj. P-value Difference (Diverse - Limited)",
       title = "Changes in Differential Abundance when\nModifications are Considered") +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 14 * TEXT_SIZE),
    axis.title      = element_text(size = 12 * TEXT_SIZE),
    axis.text       = element_text(size = 11 * TEXT_SIZE),
    legend.position = "right",
    legend.text     = element_text(size = 10 * TEXT_SIZE),
    legend.title    = element_text(size = 11 * TEXT_SIZE)
  )

ggsave(file.path(FIG_DIR, "Fig4_PctMod_vs_Discordance.png"), p_fig4_pct_mod_scatter, width = 6.5, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_PctMod_vs_Discordance.svg"), p_fig4_pct_mod_scatter, width = 6.5, height = 5)

# Limited PTM volcano, colored by ModDiff ----
highlight_limited <- merged_pro %>% filter(Gene %in% highlight_pct_mod$Gene)

p_fig4_limited_volcano <- ggplot(merged_pro,
    aes(x = logFC_NoMod, y = -log10(adj.P.Val_NoMod), color = ModDiff)) +
  geom_point(alpha = 0.8, size = 2.5) +
  scale_color_manual(values = modiff_colors) +
  scale_fill_manual(values = modiff_colors, guide = "none") +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_hline(yintercept = -log10(pval_cutoff),
             linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_point(data = highlight_limited,
             aes(x = logFC_NoMod, y = -log10(adj.P.Val_NoMod), fill = ModDiff),
             size = 4, shape = 21, color = "black", stroke = 0.8) +
  geom_text_repel(data = highlight_limited,
                  aes(x = logFC_NoMod, y = -log10(adj.P.Val_NoMod), label = Gene),
                  color = "black", size = 5, fontface = "bold",
                  box.padding = 0.5, min.segment.length = 0,
                  max.overlaps = Inf, force = 3) +
  labs(x = "Log2 Fold Change (Limited PTMs)",
       y = "-Log10 Adjusted P-value (Limited PTMs)",
       title = "Protein Differential Abundance\n(Limited PTM Analysis)") +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 14 * TEXT_SIZE),
    axis.title      = element_text(size = 12 * TEXT_SIZE),
    axis.text       = element_text(size = 11 * TEXT_SIZE),
    legend.position = "none"
  )

ggsave(file.path(FIG_DIR, "Fig4_LimitedVolcano_NMod.png"), p_fig4_limited_volcano, width = 6.5, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_LimitedVolcano_NMod.svg"), p_fig4_limited_volcano, width = 6.5, height = 5)


# Venn diagram ----
fit_euler <- euler(list("Diverse PTMs" = diverse_de_genes,
                        "Limited PTMs" = limited_de_genes))

p_fig4_venn <- plot(fit_euler,
                    quantities = list(cex = 1.3),
                    labels     = list(cex = 1.1, box = list(col = NA, fill = alpha("white", 0.6))),
                    fills      = list(fill = c("#ef5350", "#5c6bc0"), alpha = 0.45),
                    edges      = list(col = "grey20", lwd = 1.5),
                    main       = list(label = "DE Proteins: Diverse vs. Limited PTMs",
                                     just = "left", x = grid::unit(0.02, "npc")))

png(file.path(FIG_DIR, "Fig4_Venn.png"), width = 5, height = 4, units = "in", res = 300)
print(p_fig4_venn)
dev.off()

# Venn diagram — proteins with >= 2 modified peptides ----
genes_with_2mods <- n_mod_strict %>% filter(n_mod >= 1) %>% pull(Gene)

fit_euler_mod <- euler(list(
  "Diverse PTMs" = intersect(diverse_de_genes, genes_with_2mods),
  "Limited PTMs" = intersect(limited_de_genes, genes_with_2mods)
))

p_fig4_venn_mod <- plot(fit_euler_mod,
                        quantities = list(cex = 1.3),
                        labels     = list(cex = 1.1, box = list(col = NA, fill = alpha("white", 0.6))),
                        fills      = list(fill = c("#ef5350", "#5c6bc0"), alpha = 0.45),
                        edges      = list(col = "grey20", lwd = 1.5),
                        main       = list(label = "DE Proteins: Diverse vs. Limited PTMs\n(≥1 Modified Peptides)",
                                         just = "left", x = grid::unit(0.02, "npc")))

png(file.path(FIG_DIR, "Fig4_Venn_2ModPep.png"), width = 4.5, height = 3.5, units = "in", res = 300)
grid::grid.newpage()
grid::pushViewport(grid::viewport(y = 0, height = 0.88, just = "bottom"))
grid::grid.draw(p_fig4_venn_mod)
grid::popViewport()
dev.off()

