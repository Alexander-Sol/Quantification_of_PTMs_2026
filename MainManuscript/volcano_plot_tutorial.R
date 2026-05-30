# =============================================================================
# Volcano Plot Tutorial
# Builds from a minimal example to a fully customized publication figure.
# Requires: ggplot2, ggrepel
#   install.packages(c("ggplot2", "ggrepel"))
# =============================================================================

library(ggplot2)
library(ggrepel)


# =============================================================================
# SIMULATED DATA
# Replace `df` with your own data frame. It must have at least two columns:
#   - log2FC : numeric, log2 fold-change (x-axis)
#   - pvalue  : numeric, raw p-value (y-axis, -log10 transformed)
# Optional columns used in later examples:
#   - gene    : character, feature label
#   - padj    : numeric, adjusted p-value (used instead of pvalue if preferred)
# =============================================================================

set.seed(42)
n <- 500

df <- data.frame(
  gene   = paste0("Gene", seq_len(n)),
  log2FC = c(rnorm(400, mean = 0, sd = 1),          # most genes near zero
              rnorm(50,  mean =  2.5, sd = 0.5),     # upregulated
              rnorm(50,  mean = -2.5, sd = 0.5)),    # downregulated
  pvalue = c(runif(400, 0.05, 1),
              runif(50,  1e-10, 0.01),
              runif(50,  1e-10, 0.01))
)

# Add an adjusted p-value column (here just a simple BH-like scaling for demo)
df$padj <- p.adjust(df$pvalue, method = "BH")


# =============================================================================
# EXAMPLE 1 — Minimal: p-value cutoff only
# Points above the dashed line are nominally significant (p < 0.05).
# =============================================================================

p_cutoff <- 0.05

# Color by significance: TRUE = significant, FALSE = not
df$sig <- df$pvalue < p_cutoff

ggplot(df, aes(x = log2FC, y = -log10(pvalue), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +

  # Horizontal dashed line at the p-value threshold
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black") +

  # Simple two-color scale: grey = NS, red = significant
  scale_color_manual(
    values = c("FALSE" = "grey70", "TRUE" = "firebrick"),
    labels = c("NS", paste0("p < ", p_cutoff)),
    name   = NULL
  ) +

  labs(
    title = "Example 1 — p-value cutoff only",
    x     = expression(log[2]~Fold~Change),
    y     = expression(-log[10]~(p~value))
  ) +
  theme_bw()


# =============================================================================
# EXAMPLE 2 — Add a log2 fold-change cutoff
# Points are only called significant if they pass BOTH thresholds.
# The two vertical dashed lines mark ±fc_cutoff.
# =============================================================================

fc_cutoff <- 1.0   # corresponds to 2-fold change on the original scale

# Classify each point into one of three groups
df$direction <- with(df, ifelse(
  pvalue >= p_cutoff | abs(log2FC) < fc_cutoff, "NS",
  ifelse(log2FC > 0, "Up", "Down")
))

ggplot(df, aes(x = log2FC, y = -log10(pvalue), color = direction)) +
  geom_point(alpha = 0.6, size = 1.5) +

  # Significance thresholds
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +

  scale_color_manual(
    values = c("NS" = "grey70", "Up" = "firebrick", "Down" = "steelblue"),
    name   = NULL
  ) +

  labs(
    title = "Example 2 — p-value + fold-change cutoffs",
    x     = expression(log[2]~Fold~Change),
    y     = expression(-log[10]~(p~value))
  ) +
  theme_bw()


# =============================================================================
# EXAMPLE 3 — Custom highlighting + labels
# Specify a vector of gene names to highlight in a distinct color with labels.
# ggrepel automatically moves labels to avoid overlap.
# =============================================================================

# Genes you want to call out regardless of statistical thresholds
genes_of_interest <- c("Gene451", "Gene452", "Gene453", "Gene501", "Gene12")

# Add a highlight column: "highlight" if in the list, otherwise use direction
df$highlight <- ifelse(df$gene %in% genes_of_interest, "Highlight", df$direction)

# Only label the highlighted genes (passing a label column to aes)
df$label <- ifelse(df$gene %in% genes_of_interest, df$gene, NA)

ggplot(df, aes(x = log2FC, y = -log10(pvalue), color = highlight, label = label)) +
  geom_point(alpha = 0.6, size = 1.5) +

  # Draw highlighted points on top and slightly larger
  geom_point(
    data  = subset(df, highlight == "Highlight"),
    size  = 3,
    shape = 21,           # filled circle with border
    color = "black",
    fill  = "gold"
  ) +

  # Non-overlapping labels for highlighted genes
  geom_text_repel(
    na.rm         = TRUE,      # skip rows where label == NA
    size          = 3,
    box.padding   = 0.4,
    point.padding = 0.3,
    segment.color = "grey40",
    max.overlaps  = Inf        # show all labels even in dense regions
  ) +

  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +

  scale_color_manual(
    values = c(
      "NS"        = "grey70",
      "Up"        = "firebrick",
      "Down"      = "steelblue",
      "Highlight" = "gold"    # legend entry; actual color applied by geom_point above
    ),
    name = NULL
  ) +

  labs(
    title = "Example 3 — custom highlighting + labels",
    x     = expression(log[2]~Fold~Change),
    y     = expression(-log[10]~(p~value))
  ) +
  theme_bw()


# =============================================================================
# EXAMPLE 4 — Publication-ready: full customization
# Combines all of the above plus refined theme, axis limits, count annotations,
# and a cleaner legend.
# =============================================================================

# --- Count significant hits per direction for annotation text ----------------
n_up   <- sum(df$direction == "Up",   na.rm = TRUE)
n_down <- sum(df$direction == "Down", na.rm = TRUE)

# --- Custom color palette ----------------------------------------------------
pal <- c(
  "NS"        = "grey80",
  "Up"        = "#D7263D",    # vivid red
  "Down"      = "#1B4F8A",    # dark blue
  "Highlight" = "#F4A417"     # amber
)

# --- Plot --------------------------------------------------------------------
ggplot(df, aes(x = log2FC, y = -log10(pvalue), color = highlight)) +

  # All background points first
  geom_point(
    data  = subset(df, highlight != "Highlight"),
    alpha = 0.55,
    size  = 1.5
  ) +

  # Highlighted points drawn last so they sit on top
  geom_point(
    data  = subset(df, highlight == "Highlight"),
    size  = 3,
    shape = 21,
    stroke = 0.8,
    color  = "black",
    fill   = pal["Highlight"]
  ) +

  # Labels for highlighted genes
  geom_text_repel(
    data          = subset(df, highlight == "Highlight"),
    aes(label = gene),
    size          = 3,
    fontface      = "italic",
    box.padding   = 0.5,
    point.padding = 0.4,
    segment.color = "grey40",
    segment.size  = 0.3,
    max.overlaps  = Inf,
    color         = "black"   # override color aes so labels are always black
  ) +

  # Threshold lines
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "grey30", linewidth = 0.4) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey30", linewidth = 0.4) +

  # Annotate hit counts in the upper corners
  annotate("text",
    x = max(df$log2FC, na.rm = TRUE) * 0.95,
    y = max(-log10(df$pvalue), na.rm = TRUE) * 0.97,
    label = paste0("Up: ", n_up),
    hjust = 1, size = 3.5, color = pal["Up"], fontface = "bold"
  ) +
  annotate("text",
    x = min(df$log2FC, na.rm = TRUE) * 0.95,
    y = max(-log10(df$pvalue), na.rm = TRUE) * 0.97,
    label = paste0("Down: ", n_down),
    hjust = 0, size = 3.5, color = pal["Down"], fontface = "bold"
  ) +

  # Colors and legend labels
  scale_color_manual(
    values = pal,
    breaks = c("Up", "Down", "Highlight", "NS"),
    labels = c(
      paste0("Up (n=", n_up, ")"),
      paste0("Down (n=", n_down, ")"),
      "Highlight",
      "NS"
    ),
    name = NULL
  ) +

  # Axis labels using plotmath for subscripts
  labs(
    title    = "Example 4 — publication-ready",
    subtitle = paste0("Cutoffs: p < ", p_cutoff, "  |  |log₂FC| > ", fc_cutoff),
    x        = expression(log[2]~Fold~Change),
    y        = expression(-log[10]~(p~value))
  ) +

  # Clean theme
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "grey93"),
    legend.position   = "right",
    legend.key.size   = unit(0.4, "cm"),
    plot.title        = element_text(face = "bold"),
    plot.subtitle     = element_text(color = "grey40", size = 10)
  ) +

  # Lock x-axis to be symmetric around zero for visual balance
  scale_x_continuous(limits = c(-max(abs(df$log2FC)), max(abs(df$log2FC))) * 1.05)
