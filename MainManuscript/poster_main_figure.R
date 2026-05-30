# poster_main_figure.R
# Centerpiece poster figure: PTM discovery gain
#
# Panel A: Per-protein peptide sequence tracks (NEFM, NEFH) — peptides as colored
#           blocks along the protein backbone, faceted by Limited vs. Diverse PTM
#           search. Uses Diverse search data; Limited-compatible peptides
#           identified by IsLimitedSearchPeptide() + N-terminal acetylation rule.
#
# Panel B: Full-dataset stacked bar — same mod-category color palette shows
#           both the scale increase and the compositional shift.
#           Limited bar = actual LimitedPtms/QuantifiedPeptides.tsv counts.
#           Diverse bar = actual DiversePtms/QuantifiedPeptides.tsv counts.
#
# Outputs (saved separately so they can be placed independently on the poster):
#   Poster_NEFM_Track.png
#   Poster_NEFH_Track.png
#   Poster_Dataset_Bars.png

library(tidyverse)
library(cowplot)

setwd("C:/Users/Alex/Source/Repos/Quantification_of_PTMs_2026")
source("MainManuscript/CustomScripts.R")

FIG_DIR   <- "MainManuscript/Figures"
CACHE_DIR <- "MainManuscript/cache"
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Input file paths ───────────────────────────────────────────────────────────
diverse_tsv <- "Data/DiversePtms/QuantifiedPeptides.tsv"
limited_tsv <- "Data/LimitedPtms/QuantifiedPeptides.tsv"
seq_pos_tsv <- "Data/SequencePositionTable.tsv"
scripts_r   <- "MainManuscript/CustomScripts.R"

# Returns TRUE if cache_path exists and is newer than all source_paths
cache_valid <- function(cache_path, source_paths) {
  file.exists(cache_path) &&
    file.mtime(cache_path) > max(file.mtime(source_paths))
}

# ── Shared color palette (matches plotting.R + adds Limited-Search Mod) ────────
mod_colors <- c(
  "Unmodified"           = "#4575b4",
  "Limited-Search Mod"   = "#abd9e9",   # oxidation on M, deamidation on N/Q, N-term acetyl
  "Carbamidomethylation" = "#fee090",
  "Non-Enzymatic Mod"    = "#fdae61",
  "Carboxymethylation"   = "#f46d43",
  "Enzymatic Mod"        = "#d73027"
)

# Display order for legend and stacking (left → right / bottom → top)
mod_order <- names(mod_colors)

# Legend display labels (Carbamidomethylation shown as "Artifactual")
mod_labels <- c(
  "Unmodified"           = "Unmodified",
  "Limited-Search Mod"   = "Limited-Search Mod",
  "Carbamidomethylation" = "Artifactual",
  "Non-Enzymatic Mod"    = "Non-Enzymatic Mod",
  "Carboxymethylation"   = "Carboxymethylation",
  "Enzymatic Mod"        = "Enzymatic Mod"
)

# Categories that a limited PTM search could have found
# (Carbamidomethylation excluded: it is an artifactual fixed mod, not biologically selective)
limited_cats <- c("Unmodified", "Limited-Search Mod")

# ── Shared theme ───────────────────────────────────────────────────────────────
poster_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(face = "bold", size = 15),
    plot.subtitle      = element_text(color = "gray40", size = 11,
                                      margin = margin(t = 2, b = 6)),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 11)
  )


# ══════════════════════════════════════════════════════════════════════════════
# PANEL A — per-protein peptide sequence tracks (NEFM, NEFH)
# ══════════════════════════════════════════════════════════════════════════════

# ── Greedy lane assignment (no two overlapping peptides share a lane) ──────────
assign_lanes <- function(df) {
  if (nrow(df) == 0) return(mutate(df, Lane = integer(0)))
  df        <- df[order(df$Start), ]
  lane_ends <- numeric(0)
  df$Lane   <- NA_integer_
  for (i in seq_len(nrow(df))) {
    avail <- which(lane_ends < df$Start[i])
    if (length(avail) == 0) {
      df$Lane[i]  <- length(lane_ends) + 1
      lane_ends   <- c(lane_ends, df$End[i])
    } else {
      df$Lane[i]        <- avail[1]
      lane_ends[avail[1]] <- df$End[i]
    }
  }
  df
}

# Build one horizontal track panel along a protein backbone
make_track_panel <- function(data, strip_label, axis_label, vlines, backbone,
                             prot_len, show_x_axis = FALSE, show_legend = FALSE) {
  p <- ggplot() +
    geom_vline(xintercept = vlines, color = "gray87",
               linewidth = 0.35, linetype = "solid") +
    geom_rect(data    = backbone,
              mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill    = "gray80", color = NA) +
    geom_rect(data    = data,
              mapping = aes(xmin = Start - 0.5, xmax = End + 0.5,
                            ymin = Lane - 0.4,  ymax = Lane + 0.4,
                            fill = Category),
              color   = NA, alpha = 0.9) +
    annotate("text",
             x = prot_len * 0.98, y = max(data$Lane) + 0.5,
             label = sprintf("n = %d", nrow(data)),
             hjust = 1, vjust = 0, size = 3.8, fontface = "bold", color = "gray25") +
    scale_fill_manual(values = mod_colors, breaks = mod_order,
                      labels = mod_labels, drop = FALSE) +
    scale_x_continuous(labels = scales::comma,
                       expand = expansion(mult = 0.01),
                       name   = if (show_x_axis) axis_label else NULL) +
    scale_y_continuous(expand = expansion(add = 0.8)) +
    labs(y = NULL, title = strip_label) +
    poster_theme +
    theme(
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      plot.title      = element_text(face = "bold", size = 12,
                                     margin = margin(b = 4)),
      panel.border    = element_rect(color = "gray80", fill = NA,
                                     linewidth = 0.4),
      legend.position = if (show_legend) "bottom" else "none",
      plot.margin     = margin(t = 4, r = 4,
                               b = if (show_x_axis) 4 else 0, l = 4)
    )

  if (!show_x_axis)
    p <- p + theme(axis.text.x  = element_blank(),
                   axis.ticks.x = element_blank())
  p
}

# Classify one gene's peptidoforms, assign lanes, and save its Limited vs.
# Diverse sequence track. Applies the same per-peptide rules as Panel B.
build_protein_track <- function(gene,
                                axis_label = paste(gene, "Amino Acid Position")) {
  cache <- file.path(CACHE_DIR, sprintf("%s_pep.rds", tolower(gene)))

  if (cache_valid(cache, c(diverse_tsv, seq_pos_tsv, scripts_r))) {
    message(sprintf("Loading %s peptides from cache...", gene))
    pep <- readRDS(cache)
  } else {
    message(sprintf("Classifying %s peptides (cache miss)...", gene))
    seq_pos     <- read.csv(seq_pos_tsv, sep = "\t", stringsAsFactors = FALSE)
    diverse_raw <- read.csv(diverse_tsv,  sep = "\t", stringsAsFactors = FALSE)

    pep <- diverse_raw %>%
      filter(grepl("Homo sapiens", Organism),
             !grepl("CONTAMINANT", Protein.Groups, ignore.case = TRUE)) %>%
      mutate(Gene = CleanGeneNames(Gene.Names)) %>%
      filter(Gene == gene) %>%
      distinct(Sequence, .keep_all = TRUE) %>%
      mutate(
        BaseSeq  = Base.Sequence,
        AllMod   = GetModsTextWReplacement(Sequence, all_mods = TRUE) %>% sapply(paste),
        Category = ClassifyPeptideModsSpecific(AllMod),
        Category = ifelse(
          Category == "Non-Enzymatic Mod" & IsLimitedSearchPeptide(GetMods(Sequence)),
          "Limited-Search Mod",
          Category
        )
      ) %>%
      left_join(seq_pos %>% select(Base.Sequence, Start, End),
                by = c("BaseSeq" = "Base.Sequence")) %>%
      filter(!is.na(Start)) %>%
      mutate(
        # N-terminal acetylation is found in any reasonable search; treat as limited-compatible
        Category = ifelse(
          Category == "Enzymatic Mod" & grepl("Acetylation", AllMod) & Start %in% c(1, 2),
          "Limited-Search Mod",
          Category
        )
      ) %>%
      select(Sequence, BaseSeq, Category, Start, End)

    saveRDS(pep, cache)
  }

  limited_track <- pep %>%
    filter(Category %in% limited_cats) %>%
    assign_lanes() %>%
    mutate(Panel = "Limited PTM Search")

  diverse_track <- pep %>%
    assign_lanes() %>%
    mutate(Panel = "Diverse PTM Search")

  prot_len <- max(pep$End, na.rm = TRUE)
  vlines   <- seq(100, floor(prot_len / 100) * 100, by = 100)   # guides every 100 aa
  backbone <- data.frame(xmin = 1, xmax = prot_len, ymin = 0.25, ymax = 0.75)

  p_div_track <- make_track_panel(diverse_track, "Diverse PTM Search",
                                  axis_label, vlines, backbone, prot_len)
  p_lim_track <- make_track_panel(limited_track, "Limited PTM Search",
                                  axis_label, vlines, backbone, prot_len,
                                  show_x_axis = TRUE, show_legend = TRUE)

  p <- plot_grid(
    p_div_track, p_lim_track,
    ncol = 1, align = "v", axis = "lr",
    rel_heights = c(1, 1)
  )

  ggsave(file.path(FIG_DIR, sprintf("Poster_%s_Track.png", gene)),
         p, width = 13, height = 6, dpi = 300)

  message(sprintf("Panel A (%s): %d limited / %d diverse peptidoforms (%.1f×)",
                  gene, nrow(limited_track), nrow(diverse_track),
                  nrow(diverse_track) / nrow(limited_track)))
  invisible(p)
}

build_protein_track("NEFM")
build_protein_track("NEFH")


# ══════════════════════════════════════════════════════════════════════════════
# PANEL B — Full-dataset stacked bar
# ══════════════════════════════════════════════════════════════════════════════

# Classify all human, non-contaminant peptides; applies same rules as Panel A
classify_peptides <- function(df) {
  df %>%
    filter(grepl("Homo sapiens", Organism),
           !grepl("CONTAMINANT", Protein.Groups, ignore.case = TRUE)) %>%
    distinct(Sequence, .keep_all = TRUE) %>%
    mutate(
      BaseSeq  = Base.Sequence,
      AllMod   = GetModsTextWReplacement(Sequence, all_mods = TRUE) %>% sapply(paste),
      Category = ClassifyPeptideModsSpecific(AllMod),
      Category = ifelse(
        Category == "Non-Enzymatic Mod" & IsLimitedSearchPeptide(GetMods(Sequence)),
        "Limited-Search Mod",
        Category
      )
    ) %>%
    left_join(seq_pos %>% select(Base.Sequence, Start),
              by = c("BaseSeq" = "Base.Sequence")) %>%
    mutate(
      Category = ifelse(
        Category == "Enzymatic Mod" & grepl("Acetylation", AllMod) &
          !is.na(Start) & Start %in% c(1, 2),
        "Limited-Search Mod",
        Category
      )
    ) %>%
    count(Category, name = "n")
}

counts_cache <- file.path(CACHE_DIR, "poster_counts.rds")

if (cache_valid(counts_cache, c(diverse_tsv, limited_tsv, seq_pos_tsv, scripts_r))) {
  message("Loading peptide counts from cache...")
  .counts        <- readRDS(counts_cache)
  limited_counts <- .counts$limited
  diverse_counts <- .counts$diverse
} else {
  message("Classifying all peptides (cache miss)...")
  if (!exists("seq_pos"))     seq_pos     <- read.csv(seq_pos_tsv, sep = "\t", stringsAsFactors = FALSE)
  if (!exists("diverse_raw")) diverse_raw <- read.csv(diverse_tsv,  sep = "\t", stringsAsFactors = FALSE)
  limited_raw    <- read.csv(limited_tsv, sep = "\t", stringsAsFactors = FALSE)

  limited_counts <- classify_peptides(limited_raw) %>% mutate(Panel = "Limited PTM Search")
  diverse_counts <- classify_peptides(diverse_raw) %>% mutate(Panel = "Diverse PTM Search")

  saveRDS(list(limited = limited_counts, diverse = diverse_counts), counts_cache)
}

bar_data <- bind_rows(limited_counts, diverse_counts) %>%
  mutate(
    Panel    = factor(Panel, levels = c("Limited PTM Search", "Diverse PTM Search")),
    Category = factor(Category, levels = rev(mod_order))
  )

totals_B <- bar_data %>%
  group_by(Panel) %>%
  summarise(Total = sum(n), .groups = "drop")

p_B <- ggplot(bar_data, aes(x = n, y = Panel, fill = Category)) +
  geom_col(position = "stack", width = 0.55) +
  geom_text(
    data    = totals_B,
    mapping = aes(x = Total, y = Panel,
                  label = paste0(scales::comma(Total), " peptidoforms"),
                  fill  = NULL),
    hjust   = -0.05, size = 4.5, fontface = "bold", color = "gray20"
  ) +
  scale_fill_manual(values = mod_colors,
                    breaks = mod_order,
                    labels = mod_labels,
                    drop   = FALSE) +
  scale_x_continuous(
    name   = "Unique Peptidoforms Detected",
    expand = expansion(mult = c(0, 0.28)),
    labels = scales::comma
  ) +
  labs(
    y        = NULL,
    title    = "Peptidoform Discovery",
    subtitle = "Diverse PTM search recovers more peptidoforms across every modification class"
  ) +
  poster_theme +
  theme(
    axis.text.y  = element_text(size = 13, face = "bold"),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.4)
  )

ggsave(file.path(FIG_DIR, "Poster_Dataset_Bars.png"),
       p_B, width = 11, height = 3.5, dpi = 300)

fold_pB <- totals_B$Total[totals_B$Panel == "Diverse PTM Search"] /
           totals_B$Total[totals_B$Panel == "Limited PTM Search"]
message(sprintf("Panel B: %d limited / %d diverse total peptidoforms (%.1f×)",
                totals_B$Total[1], totals_B$Total[2], fold_pB))
message("Saved:")
message("  Poster_NEFM_Track.png")
message("  Poster_NEFH_Track.png")
message("  Poster_Dataset_Bars.png")
