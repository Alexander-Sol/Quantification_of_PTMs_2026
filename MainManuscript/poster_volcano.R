# poster_volcano.R
# Poster versions of the peptidoform-fraction figures for the main poster figure.
#
#   Panel 1: Peptidoform volcano ("Differing Peptidoform Fractions") — wider and
#            larger than the manuscript Fig3_Volcano, with bigger points, text,
#            and legend for at-a-distance legibility.
#   Panel 2: Peptide category composition by stringency — the same stacked bar as
#            the manuscript Fig2_CategoryBar, but laid out HORIZONTALLY (groups as
#            rows) so it sits naturally beside the wide volcano on the poster.
#
# Depends on objects built by analysis.R (base_filtered, pval_cutoff, logFC_cutoff),
# exactly like plotting.R. analysis.R is sourced automatically if those objects are
# not already in the session; its slow steps are .rds-cached, so reruns are cheap.
#
# Outputs (MainManuscript/Figures/):
#   Poster_Volcano.png / .svg
#   Poster_CategoryBar_Horizontal.png / .svg

library(tidyverse)
library(scales)

setwd("C:/Users/Alex/Source/Repos/Quantification_of_PTMs_2026")
source("MainManuscript/CustomScripts.R")

FIG_DIR <- "MainManuscript/Figures"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Ensure analysis objects are available ─────────────────────────────────────
# analysis.R caches its slow intermediates as .rds, so this is cheap on reruns.
if (!exists("base_filtered") || !exists("pval_cutoff") || !exists("logFC_cutoff")) {
  message("Analysis objects not found in session; sourcing analysis.R ...")
  source("MainManuscript/analysis.R")
}

# ── Shared category palette (matches plotting.R) ──────────────────────────────
mod_colors <- c(
  "Unmodified"           = "#4575b4",
  "Carbamidomethylation" = "#fee090",
  "Non-Enzymatic Mod"    = "#fdae61",
  "Carboxymethylation"   = "#f46d43",
  "Enzymatic Mod"        = "#d73027"
)

# Legend display labels (Carbamidomethylation shown as "Artifactual")
mod_labels <- c(
  "Unmodified"           = "Unmodified",
  "Carbamidomethylation" = "Artifactual",
  "Non-Enzymatic Mod"    = "Non-Enzymatic Mod",
  "Carboxymethylation"   = "Carboxymethylation",
  "Enzymatic Mod"        = "Enzymatic Mod"
)


# ══════════════════════════════════════════════════════════════════════════════
# PANEL 1 — Peptidoform volcano (wider + larger for the poster)
# ══════════════════════════════════════════════════════════════════════════════

volcano_df <- base_filtered %>%
  mutate(
    InTop1000 = adj.P.Val <= pval_cutoff & abs(logFC) >= logFC_cutoff,
    Highlight = case_when(
      InTop1000 & Gene %in% c("NEFH", "NEFM", "NEFL")                                ~ "Neurofilaments",
      InTop1000 & Gene == "GFAP"                                                       ~ "GFAP",
      InTop1000 & Gene %in% c("EEF1G", "EEF1B2", "EEF1D", "EEF1A1", "EEF1A2", "EEF2")  ~ "Elongation Factors",
      TRUE ~ NA_character_
    ),
    y_display = -log10(adj.P.Val)
  )

highlight_df <- volcano_df %>%
  filter(!is.na(Highlight)) %>%
  mutate(Highlight = factor(Highlight,
                            levels = c("GFAP", "Elongation Factors", "Neurofilaments")))

sig_y <- -log10(pval_cutoff)

# Compress the upper tail so a few very small p-values don't stretch the axis.
upper_compress_trans <- scales::trans_new(
  name      = "upper_compress",
  transform = function(x) ifelse(x <= 5, x, 5 + log1p(x - 5) * 1),
  inverse   = function(x) ifelse(x <= 5, x, 5 + expm1((x - 5) / 1))
)

p_poster_volcano <- ggplot(volcano_df, aes(x = logFC, y = y_display)) +
  geom_point(data = filter(volcano_df, !InTop1000, is.na(Highlight)),
             color = "grey75", alpha = 0.4, size = 2.4, shape = 19) +
  geom_point(data = filter(volcano_df, InTop1000, is.na(Highlight)),
             aes(color = Category), alpha = 0.85, size = 2.9, shape = 19) +
  geom_point(data = filter(highlight_df, Highlight == "GFAP"),
             aes(fill = Category, shape = Highlight),
             size = 3.8, alpha = 0.9, stroke = 0.7, color = "black") +
  geom_point(data = filter(highlight_df, Highlight == "Elongation Factors"),
             aes(fill = Category, shape = Highlight),
             size = 4.6, alpha = 0.9, stroke = 0.7, color = "black") +
  geom_point(data = filter(highlight_df, Highlight == "Neurofilaments"),
             aes(fill = Category, shape = Highlight),
             size = 4.8, alpha = 0.9, stroke = 0.7, color = "black") +
  geom_hline(yintercept = sig_y,
             color = "grey30", linewidth = 0.6, linetype = "dashed") +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff),
             color = "grey30", linewidth = 0.6, linetype = "dashed") +
  scale_color_manual(values = mod_colors, breaks = names(mod_colors),
                     labels = mod_labels) +
  scale_fill_manual(values = mod_colors, guide = "none") +
  scale_shape_manual(values = c("GFAP" = 24, "Elongation Factors" = 22,
                                "Neurofilaments" = 23), name = NULL) +
  scale_y_continuous(trans  = upper_compress_trans,
                     limits = c(0, 12.5),
                     breaks = c(0, 1, 2, 3, 4, 5, 7, 10, 12)) +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value",
       title = "Differing Peptidoform Fractions") +
  theme_minimal(base_size = 18) +
  theme(
    plot.title           = element_text(size = 22, face = "bold"),
    axis.title           = element_text(size = 18),
    axis.text            = element_text(size = 16),
    panel.grid.minor     = element_blank(),
    legend.position      = c(0.46, 0.99),
    legend.justification = c(0.5, 1),
    legend.background    = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.text          = element_text(size = 15),
    legend.title         = element_blank(),
    legend.box           = "vertical",
    legend.spacing.y     = unit(2, "pt")
  ) +
  guides(color = guide_legend(order = 1, reverse = TRUE,
                              override.aes = list(size = 4.5)),
         shape = guide_legend(order = 2, override.aes = list(size = 4.5)))

ggsave(file.path(FIG_DIR, "Poster_Volcano.png"),
       p_poster_volcano, width = 13, height = 8, dpi = 300)
ggsave(file.path(FIG_DIR, "Poster_Volcano.svg"),
       p_poster_volcano, width = 13, height = 8)


# ══════════════════════════════════════════════════════════════════════════════
# PANEL 2 — Peptide category composition by stringency (HORIZONTAL)
# ══════════════════════════════════════════════════════════════════════════════

category_levels_ordered <- c("Enzymatic Mod", "Carboxymethylation",
                             "Non-Enzymatic Mod", "Carbamidomethylation",
                             "Unmodified")

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
# Reverse the group order so, in a horizontal layout, the first stringency level
# ("All Peptides") reads at the TOP and the strictest is at the bottom.
group_levels <- map_chr(stringency_defs, ~ .x$label)
bar_df$Group <- factor(bar_df$Group, levels = rev(group_levels))

small_threshold <- 5  # minimum percent for an in-bar label, to avoid clutter

bar_totals <- bar_df %>%
  group_by(Group) %>%
  summarise(n_total = sum(Count), .groups = "drop")

p_poster_bar <- ggplot(bar_df, aes(x = percent, y = Group, fill = Category)) +
  geom_col(position = "stack", alpha = 0.85, color = "black", linewidth = 0.4) +
  geom_text(
    data = filter(bar_df, percent >= small_threshold),
    aes(label = sprintf("%d%%", round(percent))),
    position = position_stack(vjust = 0.5),
    size = 5.5, color = "black"
  ) +
  geom_text(
    data = bar_totals,
    aes(x = 100, y = Group, label = paste0("n=", scales::comma(n_total))),
    hjust = -0.08, size = 5, color = "gray20", inherit.aes = FALSE
  ) +
  scale_fill_manual(values = mod_colors, breaks = category_levels_ordered,
                    labels = mod_labels) +
  scale_x_continuous(name   = "Percent of Peptides",
                     limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.13))) +
  labs(y = NULL, title = "Peptide Category Composition by Stringency") +
  theme_minimal(base_size = 18) +
  theme(
    plot.title       = element_text(size = 22, face = "bold"),
    axis.title.x     = element_text(size = 18),
    axis.text.y      = element_text(size = 17, face = "bold"),
    axis.text.x      = element_text(size = 15),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 15),
    legend.title     = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

ggsave(file.path(FIG_DIR, "Poster_CategoryBar_Horizontal.png"),
       p_poster_bar, width = 12, height = 4.5, dpi = 300)
ggsave(file.path(FIG_DIR, "Poster_CategoryBar_Horizontal.svg"),
       p_poster_bar, width = 12, height = 4.5)

message("Saved:")
message("  Poster_Volcano.png / .svg                 (13 x 8 in)")
message("  Poster_CategoryBar_Horizontal.png / .svg  (12 x 4.5 in)")
