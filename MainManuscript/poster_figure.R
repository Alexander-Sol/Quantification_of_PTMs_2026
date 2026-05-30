# poster_figure.R
# Centerpiece poster figure: peptidoform discovery gain from diverse PTM search.
# Standalone — reads raw QuantifiedPeptides.tsv files directly, no analysis.R needed.
#
# Outputs:
#   Poster_PeptideDiscovery_Log.png   — log x-axis; both panels A + B
#   Poster_PeptideDiscovery_Facet.png — free-scale facets by family; both panels A + B

library(tidyverse)
library(cowplot)

setwd("C:/Users/Alex/Source/Repos/Quantification_of_PTMs_2026")
source("MainManuscript/CustomScripts.R")

FIG_DIR <- "MainManuscript/Figures"

# ── Target proteins and family groupings ──────────────────────────────────────
target_genes <- c("NEFL", "NEFM", "NEFH",
                  "EEF1A1", "EEF1A2", "EEF1B2", "EEF1D", "EEF1G",
                  "GFAP")

protein_family <- c(
  NEFL   = "Neurofilaments", NEFM = "Neurofilaments", NEFH = "Neurofilaments",
  EEF1A1 = "EEF1 Family",   EEF1A2 = "EEF1 Family", EEF1B2 = "EEF1 Family",
  EEF1D  = "EEF1 Family",   EEF1G  = "EEF1 Family",
  GFAP   = "GFAP"
)

# Display order: bottom → top within combined y-axis
gene_order <- c("GFAP", "EEF1G", "EEF1D", "EEF1B2", "EEF1A2", "EEF1A1",
                "NEFH", "NEFM", "NEFL")

# ── Data loading ──────────────────────────────────────────────────────────────
load_pep <- function(path, search_label) {
  df <- read.csv(path, sep = "\t", stringsAsFactors = FALSE)
  df %>%
    filter(grepl("Homo sapiens", Organism)) %>%
    # Use primary gene assignment; shared paralogs (e.g. EEF1A1/EEF1A2) count
    # for whichever protein is listed first in Gene.Names.
    mutate(Gene = CleanGeneNames(Gene.Names)) %>%
    filter(Gene %in% target_genes) %>%
    group_by(Gene) %>%
    summarise(
      Peptidoforms  = n_distinct(Sequence),
      BaseSequences = n_distinct(Base.Sequence),
      .groups = "drop"
    ) %>%
    mutate(Search = search_label)
}

diverse_counts <- load_pep("Data/DiversePtms/QuantifiedPeptides.tsv", "Diverse PTMs")
limited_counts <- load_pep("Data/LimitedPtms/QuantifiedPeptides.tsv", "Limited PTMs")

# ── Combined dumbbell data ────────────────────────────────────────────────────
dumbbell <- diverse_counts %>%
  select(Gene, DiversePep = Peptidoforms, DiverseBase = BaseSequences) %>%
  left_join(
    limited_counts %>% select(Gene, LimitedPep = Peptidoforms, LimitedBase = BaseSequences),
    by = "Gene"
  ) %>%
  replace_na(list(LimitedPep = 0, LimitedBase = 0)) %>%
  mutate(
    FoldChange     = DiversePep  / pmax(LimitedPep,  1),
    BaseFoldChange = DiverseBase / pmax(LimitedBase, 1),
    Family = factor(protein_family[Gene],
                    levels = c("Neurofilaments", "EEF1 Family", "GFAP")),
    Gene   = factor(Gene, levels = gene_order)
  )

# ── Shared colors / theme ─────────────────────────────────────────────────────
# UW-Madison brand colors
# Dark Red (#9B0000): secondary digital
# Med Gray Blue (#6B9999): accent, converted from CMYK 30,0,0,40
# Yellow (#FFBF00): accent, converted from CMYK 0,25,100,0
family_colors <- c(
  "Neurofilaments" = "#9B0000",
  "EEF1 Family"    = "#6B9999",
  "GFAP"           = "#FFBF00"
)

base_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_line(color = "gray93"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "gray88"),
    axis.text.y        = element_text(size = 13, face = "bold"),
    axis.text.x        = element_text(size = 11),
    axis.title.x       = element_text(size = 12, margin = margin(t = 6)),
    plot.title         = element_text(size = 14, face = "bold"),
    plot.subtitle      = element_text(size = 10, color = "gray40",
                                      margin = margin(t = 3, b = 6)),
    legend.position    = "none"
  )

# Subtitle with open/filled circle legend
pep_subtitle <- "○  Limited PTM search     ●  Diverse PTM search"
base_subtitle <- "○  Limited     ●  Diverse  (backbone sequences only)"


# ══════════════════════════════════════════════════════════════════════════════
# VERSION 1: Log x-axis (single combined panel, log10 scale)
# ══════════════════════════════════════════════════════════════════════════════

make_dumbbell <- function(data, x_start, x_end, fold_col,
                          x_label, title_txt, subtitle_txt,
                          log_scale = FALSE,
                          show_y = TRUE,
                          fold_threshold = 1.1,
                          x_expand = expansion(mult = c(0.02, 0.12))) {

  p <- ggplot() +
    geom_segment(
      data = data,
      aes(y = Gene, yend = Gene,
          x = .data[[x_start]], xend = .data[[x_end]],
          color = Family),
      linewidth = 1.4, alpha = 0.55
    ) +
    geom_point(
      data = data,
      aes(y = Gene, x = .data[[x_start]], color = Family),
      shape = 21, size = 4, fill = "white", stroke = 1.8
    ) +
    geom_point(
      data = data,
      aes(y = Gene, x = .data[[x_end]], color = Family),
      shape = 19, size = 5.5
    ) +
    geom_text(
      data = data %>% filter(.data[[fold_col]] > fold_threshold),
      aes(y = Gene, x = .data[[x_end]],
          label = sprintf("%.1f×", .data[[fold_col]]),
          color = Family),
      hjust = -0.25, fontface = "bold", size = 4
    ) +
    scale_color_manual(values = family_colors) +
    scale_x_continuous(name = x_label, expand = x_expand) +
    labs(y = NULL, title = title_txt, subtitle = subtitle_txt) +
    base_theme

  if (log_scale)
    p <- p + scale_x_log10(name = x_label, expand = x_expand,
                           labels = scales::comma)

  if (!show_y)
    p <- p + theme(axis.text.y = element_blank())

  p
}

# Panel A: peptidoforms, log scale
p_A_log <- make_dumbbell(
  dumbbell, "LimitedPep", "DiversePep", "FoldChange",
  x_label    = "Unique Peptidoforms Detected (log scale)",
  title_txt  = "More Peptidoforms with Modification-Aware Search",
  subtitle_txt = pep_subtitle,
  log_scale  = TRUE
)

# Panel B: base sequences, log scale
p_B_log <- make_dumbbell(
  dumbbell, "LimitedBase", "DiverseBase", "BaseFoldChange",
  x_label    = "Unique Base Sequences Detected (log scale)",
  title_txt  = "Backbone Sequence Coverage",
  subtitle_txt = base_subtitle,
  log_scale  = TRUE,
  show_y     = FALSE,
  fold_threshold = 1.05
)

p_log <- plot_grid(p_A_log, p_B_log, nrow = 1, rel_widths = c(1.5, 1),
                   labels = c("A", "B"), label_size = 16)

ggsave(file.path(FIG_DIR, "Poster_PeptideDiscovery_Log.png"),
       p_log, width = 14, height = 6, dpi = 300)
ggsave(file.path(FIG_DIR, "Poster_PeptideDiscovery_Log.svg"),
       p_log, width = 14, height = 6)


# ══════════════════════════════════════════════════════════════════════════════
# VERSION 2: Faceted by protein family (free x-scale per family)
# ══════════════════════════════════════════════════════════════════════════════

make_dumbbell_facet <- function(data, x_start, x_end, fold_col,
                                x_label, title_txt, subtitle_txt,
                                show_y = TRUE,
                                fold_threshold = 1.1,
                                x_expand = expansion(mult = c(0.02, 0.18))) {

  ggplot() +
    geom_segment(
      data = data,
      aes(y = Gene, yend = Gene,
          x = .data[[x_start]], xend = .data[[x_end]],
          color = Family),
      linewidth = 1.4, alpha = 0.55
    ) +
    geom_point(
      data = data,
      aes(y = Gene, x = .data[[x_start]], color = Family),
      shape = 21, size = 4, fill = "white", stroke = 1.8
    ) +
    geom_point(
      data = data,
      aes(y = Gene, x = .data[[x_end]], color = Family),
      shape = 19, size = 5.5
    ) +
    geom_text(
      data = data %>% filter(.data[[fold_col]] > fold_threshold),
      aes(y = Gene, x = .data[[x_end]],
          label = sprintf("%.1f×", .data[[fold_col]]),
          color = Family),
      hjust = -0.3, fontface = "bold", size = 4
    ) +
    facet_wrap(~ Family, ncol = 1, scales = "free",
               strip.position = "left") +
    scale_color_manual(values = family_colors) +
    scale_x_continuous(name = x_label,
                       expand = x_expand,
                       labels = scales::comma) +
    labs(y = NULL, title = title_txt, subtitle = subtitle_txt) +
    base_theme +
    theme(
      strip.text.y.left  = element_text(size = 12, face = "bold", angle = 90),
      strip.placement    = "outside",
      panel.spacing.y    = unit(12, "pt"),
      panel.border       = element_rect(color = "gray80", fill = NA,
                                        linewidth = 0.5)
    ) +
    { if (!show_y) theme(axis.text.y = element_blank()) }
}

p_A_facet <- make_dumbbell_facet(
  dumbbell, "LimitedPep", "DiversePep", "FoldChange",
  x_label    = "Unique Peptidoforms Detected",
  title_txt  = "More Peptidoforms with Modification-Aware Search",
  subtitle_txt = pep_subtitle
)

p_B_facet <- make_dumbbell_facet(
  dumbbell, "LimitedBase", "DiverseBase", "BaseFoldChange",
  x_label    = "Unique Base Sequences Detected",
  title_txt  = "Backbone Sequence Coverage",
  subtitle_txt = base_subtitle,
  show_y     = FALSE,
  fold_threshold = 1.05,
  x_expand   = expansion(mult = c(0.10, 0.18))
)

p_facet <- plot_grid(p_A_facet, p_B_facet, nrow = 1, rel_widths = c(1.5, 1),
                     labels = c("A", "B"), label_size = 16)

ggsave(file.path(FIG_DIR, "Poster_PeptideDiscovery_Facet.png"),
       p_facet, width = 14, height = 7, dpi = 300)
ggsave(file.path(FIG_DIR, "Poster_PeptideDiscovery_Facet.svg"),
       p_facet, width = 14, height = 7)

message("Saved:")
message("  Poster_PeptideDiscovery_Log.png/.svg")
message("  Poster_PeptideDiscovery_Facet.png/.svg")
