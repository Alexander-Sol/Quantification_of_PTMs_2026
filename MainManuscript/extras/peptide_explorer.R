# peptide_explorer.R
# Diagnostic plots for inspecting all peptides from a single protein.
# Requires analysis.R to have been sourced first.

library(tidyverse)
library(cowplot)

mod_colors <- c(
  "Unmodified"           = "#4575b4",
  "Carbamidomethylation" = "#fee090",
  "Non-Enzymatic Mod"    = "#fdae61",
  "Carboxymethylation"   = "#f46d43",
  "Enzymatic Mod"        = "#d73027"
)

bio_mods_of_interest <- c(
  "Citrullination", "Acetylation", "Phosphorylation",
  "Methylation", "Dimethylation", "Trimethylation",
  "Hydroxyproline", "Glutarylation", "Butyrylation",
  "Nitrosylation"
)

bio_protein_colors <- c(
  "Citrullination"      = "darkblue",
  "Acetylation"         = "forestgreen",
  "Methylation"         = "darkorchid1",
  "Dimethylation"       = "darkorchid",
  "Trimethylation"      = "darkorchid4",
  "Phosphorylation"     = "cyan3",
  "Hydroxyproline"      = "deeppink3",
  "Glutarylation"       = "sienna",
  "Butyrylation"        = "sienna2",
  "Nitrosylation"       = "firebrick4",
  "Oxidation"           = "steelblue4",
  "Carboxymethyllysine" = "#f46d43",
  "Other"               = "grey88"
)


# ── Plot parameters ───────────────────────────────────────────────────────────
PEP_X_LAB          <- "Log2 Fold Change"
PEP_Y_LAB_PVAL     <- "-Log10 Adjusted P-value"
PEP_Y_LAB_S2       <- "Posterior variance (s2.post)"
PEP_PT_SIZE        <- 2.75
PEP_TITLE_SIZE     <- 13   # subplot title (pt)
PEP_AXIS_TITLE_SIZE <- 10   # axis label (pt)
PEP_AXIS_TEXT_SIZE  <- 9   # axis tick text (pt)


# ── EnzMod assignment (shared) ────────────────────────────────────────────────

.assign_enzmod <- function(all_mod_str, category) {
  if (category == "Enzymatic Mod") {
    mods <- trimws(unlist(strsplit(all_mod_str, ", ")))
    hit  <- bio_mods_of_interest[bio_mods_of_interest %in% mods]
    if (length(hit) > 0) return(hit[1])
  }
  if (grepl("Carboxymethyllysine", all_mod_str, fixed = TRUE)) return("Carboxymethyllysine")
  if (grepl("Nitrosylation",       all_mod_str, fixed = TRUE)) return("Nitrosylation")
  if (grepl("Oxidation",           all_mod_str, fixed = TRUE)) return("Oxidation")
  "Other"
}

.enzmod_scale <- scale_color_manual(
  values = bio_protein_colors, drop = TRUE,
  breaks = c("Acetylation", "Butyrylation", "Carboxymethyllysine",
             "Citrullination", "Dimethylation", "Glutarylation",
             "Hydroxyproline", "Methylation", "Nitrosylation",
             "Oxidation", "Phosphorylation", "Trimethylation",
             "Other")
)


# ── Shared data prep ───────────────────────────────────────────────────────────

.prep_pep_rows <- function(gene_name, table, eb_fit, min_prop_obs = 0) {
  pep_rows <- table[table$Gene == gene_name, ]
  pep_rows <- pep_rows[!is.na(pep_rows$logFC), ]
  if (min_prop_obs > 0) pep_rows <- pep_rows[pep_rows$PropObs >= min_prop_obs, ]
  if (nrow(pep_rows) == 0) stop("No peptides found for: ", gene_name)
  pep_rows$s2.post        <- eb_fit$s2.post[match(rownames(pep_rows), rownames(eb_fit))]
  pep_rows$neg_log10_pval <- -log10(pep_rows$adj.P.Val)
  pep_rows
}

.base_aes <- list(
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60"),
  labs(x = PEP_X_LAB),
  theme_bw(),
  theme(
    plot.title = element_text(size = PEP_TITLE_SIZE),
    axis.title = element_text(size = PEP_AXIS_TITLE_SIZE),
    axis.text  = element_text(size = PEP_AXIS_TEXT_SIZE)
  )
)


# ── Individual plots ───────────────────────────────────────────────────────────

plot_protein_peptides <- function(gene_name,
                                  table        = table_occ,
                                  eb_fit       = eb_fit_occ,
                                  use_pval     = FALSE,
                                  min_prop_obs = 0) {
  pep_rows <- .prep_pep_rows(gene_name, table, eb_fit, min_prop_obs)
  y_var    <- if (use_pval) "neg_log10_pval" else "s2.post"
  y_label  <- if (use_pval) PEP_Y_LAB_PVAL else PEP_Y_LAB_S2

  ggplot(pep_rows, aes(x = logFC, y = .data[[y_var]],
                       color = PropObs, label = PeptideSequence)) +
    geom_point(size = PEP_PT_SIZE, alpha = 0.85) +
    .base_aes +
    scale_color_viridis_c(option = "plasma", limits = c(0, 1)) +
    labs(title = "Proportion Observed", color = "PropObs", y = y_label)
}


plot_protein_peptides_mod <- function(gene_name,
                                      table        = table_occ,
                                      eb_fit       = eb_fit_occ,
                                      use_pval     = FALSE,
                                      min_prop_obs = 0) {
  pep_rows <- .prep_pep_rows(gene_name, table, eb_fit, min_prop_obs)
  y_var    <- if (use_pval) "neg_log10_pval" else "s2.post"
  y_label  <- if (use_pval) PEP_Y_LAB_PVAL else PEP_Y_LAB_S2

  ggplot(pep_rows, aes(x = logFC, y = .data[[y_var]],
                       color = Category, label = PeptideSequence)) +
    geom_point(size = PEP_PT_SIZE, alpha = 0.85) +
    .base_aes +
    scale_color_manual(values = mod_colors, na.value = "gray70",
                       breaks = c("Enzymatic Mod", "Carboxymethylation",
                                  "Non-Enzymatic Mod", "Carbamidomethylation",
                                  "Unmodified")) +
    labs(title = "Modification Category", color = "Category", y = y_label)
}


plot_protein_peptides_intensity <- function(gene_name,
                                            table        = table_occ,
                                            eb_fit       = eb_fit_occ,
                                            y.pep        = y.peptide,
                                            use_pval     = FALSE,
                                            min_prop_obs = 0) {
  pep_rows <- .prep_pep_rows(gene_name, table, eb_fit, min_prop_obs)
  y_var    <- if (use_pval) "neg_log10_pval" else "s2.post"
  y_label  <- if (use_pval) PEP_Y_LAB_PVAL else PEP_Y_LAB_S2

  idx <- match(pep_rows$PeptideSequence, y.pep$genes$PeptideSequence)
  pep_rows$mean_intensity <- rowMeans(y.pep$E[idx, , drop = FALSE], na.rm = TRUE)

  ggplot(pep_rows, aes(x = logFC, y = .data[[y_var]],
                       color = mean_intensity, label = PeptideSequence)) +
    geom_point(size = PEP_PT_SIZE, alpha = 0.85) +
    .base_aes +
    scale_color_viridis_c(option = "viridis", na.value = "gray70") +
    labs(title = "Mean Intensity", color = "Mean\nIntensity", y = y_label)
}


plot_protein_peptides_groups <- function(gene_name,
                                         table        = table_occ,
                                         eb_fit       = eb_fit_occ,
                                         use_pval     = FALSE,
                                         min_prop_obs = 0) {
  pep_rows <- .prep_pep_rows(gene_name, table, eb_fit, min_prop_obs)
  y_var    <- if (use_pval) "neg_log10_pval" else "s2.post"
  y_label  <- if (use_pval) PEP_Y_LAB_PVAL else PEP_Y_LAB_S2

  seqs <- pep_rows$BaseSequence
  n    <- length(seqs)
  adj  <- matrix(FALSE, n, n)
  for (i in seq_len(n))
    for (j in seq_len(n))
      if (i != j)
        adj[i, j] <- grepl(seqs[i], seqs[j], fixed = TRUE) |
                     grepl(seqs[j], seqs[i], fixed = TRUE)

  group_id <- integer(n)
  gid <- 0L
  for (i in seq_len(n)) {
    if (group_id[i] == 0L) {
      gid <- gid + 1L
      queue <- i
      while (length(queue) > 0) {
        node  <- queue[1]; queue <- queue[-1]
        if (group_id[node] == 0L) {
          group_id[node] <- gid
          queue <- c(queue, which(adj[node, ] & group_id == 0L))
        }
      }
    }
  }

  pep_rows$group_id <- group_id

  group_ranks <- names(sort(table(group_id), decreasing = TRUE))
  top6        <- as.integer(group_ranks[seq_len(min(6, length(group_ranks)))])
  has_other   <- length(group_ranks) > 6

  pep_rows$SeqGroup <- factor(
    ifelse(pep_rows$group_id %in% top6,
           as.character(match(pep_rows$group_id, top6)),
           "Other"),
    levels = c(if (has_other) "Other", as.character(seq_len(length(top6))))
  )

  pep_rows <- pep_rows[order(pep_rows$SeqGroup, decreasing = FALSE), ]

  n_top <- length(top6)
  pal   <- c(if (has_other) c("Other" = "gray70"),
             setNames(scales::hue_pal()(n_top), as.character(seq_len(n_top))))

  ggplot(pep_rows, aes(x = logFC, y = .data[[y_var]],
                       color = SeqGroup, label = PeptideSequence)) +
    geom_point(size = PEP_PT_SIZE, alpha = 0.85) +
    .base_aes +
    scale_color_manual(values = pal) +
    labs(title = "Sequence Groups", color = "Sequence\nGroup", y = y_label)
}


plot_protein_peptides_enzmod <- function(gene_name,
                                          table        = table_occ,
                                          eb_fit       = eb_fit_occ,
                                          use_pval     = FALSE,
                                          min_prop_obs = 0) {
  pep_rows <- .prep_pep_rows(gene_name, table, eb_fit, min_prop_obs)
  y_var    <- if (use_pval) "neg_log10_pval" else "s2.post"
  y_label  <- if (use_pval) PEP_Y_LAB_PVAL else PEP_Y_LAB_S2

  pep_rows$EnzMod <- factor(
    mapply(.assign_enzmod, pep_rows$AllMod, pep_rows$Category),
    levels = c("Other", bio_mods_of_interest, "Carboxymethyllysine", "Oxidation")
  )

  # Draw Other behind colored points
  pep_rows <- pep_rows[order(pep_rows$EnzMod), ]

  ggplot(pep_rows, aes(x = logFC, y = .data[[y_var]],
                       color = EnzMod, label = PeptideSequence)) +
    geom_point(data = \(d) d[d$EnzMod == "Other", ], size = PEP_PT_SIZE, alpha = .7) +
    geom_point(data = \(d) d[d$EnzMod != "Other", ], size = PEP_PT_SIZE, alpha = 0.85) +
    .base_aes +
    .enzmod_scale +
    labs(title = "Modification Type", color = "Mod", y = y_label)
}


# ── Combined grid ──────────────────────────────────────────────────────────────

plot_protein_all <- function(gene_name,
                             table        = table_occ,
                             eb_fit       = eb_fit_occ,
                             y.pep        = y.peptide,
                             use_pval     = FALSE,
                             min_prop_obs = 0) {
  sig_line <- if (use_pval)
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "firebrick", alpha = 0.6)

  p1 <- plot_protein_peptides(gene_name, table, eb_fit, use_pval, min_prop_obs) + sig_line
  p2 <- plot_protein_peptides_mod(gene_name, table, eb_fit, use_pval, min_prop_obs) + sig_line
  p3 <- plot_protein_peptides_intensity(gene_name, table, eb_fit, y.pep, use_pval, min_prop_obs) + sig_line
  p4 <- plot_protein_peptides_enzmod(gene_name, table, eb_fit, use_pval, min_prop_obs) + sig_line

  title <- ggdraw() +
    draw_label(paste("Peptides for", gene_name), fontface = "bold", size = 14)

  grid <- plot_grid(p1, p2, p3, p4, nrow = 2, ncol = 2)
  plot_grid(title, grid, ncol = 1, rel_heights = c(0.04, 1))
}

p_supp_eef1g <- plot_protein_all("EEF1G", use_pval = TRUE, min_prop_obs = 0.12)
ggsave("Supplemental/Figures/SupplFig_PeptideExplorer_EEF1G.png",
       p_supp_eef1g, width = 14, height = 10, dpi = 300)


# ── Generic 3-gene cat/mod panel ──────────────────────────────────────────────

.plot_gene_panel <- function(genes, panel_title,
                             table        = table_occ,
                             eb_fit       = eb_fit_occ,
                             use_pval     = FALSE,
                             min_prop_obs = 0) {
  sig_line <- if (use_pval)
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "firebrick", alpha = 0.6)

  # Shared y limits across all genes and both plot types
  y_var  <- if (use_pval) "neg_log10_pval" else "s2.post"
  all_y  <- unlist(lapply(genes, function(g) {
    rows <- tryCatch(.prep_pep_rows(g, table, eb_fit, min_prop_obs), error = function(e) NULL)
    if (is.null(rows)) return(numeric(0))
    rows[[y_var]][is.finite(rows[[y_var]])]
  }))
  y_pad  <- diff(range(all_y)) * 0.05
  y_lim  <- range(all_y) + c(-y_pad, y_pad)

  cat_plots <- lapply(genes, function(g)
    plot_protein_peptides_mod(g, table, eb_fit, use_pval, min_prop_obs) +
      labs(title = paste0("Peptidoforms from ", g)) + sig_line +
      coord_cartesian(ylim = y_lim))

  mod_plots <- lapply(genes, function(g)
    plot_protein_peptides_enzmod(g, table, eb_fit, use_pval, min_prop_obs) +
      labs(title = NULL) + sig_line +
      coord_cartesian(ylim = y_lim))

  # Build mod legend from combined data so all present mods are represented
  all_rows <- do.call(rbind, lapply(genes, function(g)
    tryCatch(.prep_pep_rows(g, table, eb_fit, min_prop_obs), error = function(e) NULL)
  ))
  all_rows$EnzMod <- factor(
    mapply(.assign_enzmod, all_rows$AllMod, all_rows$Category),
    levels = c("Other", bio_mods_of_interest, "Carboxymethyllysine", "Oxidation")
  )
  mod_legend <- get_legend(
    ggplot(all_rows, aes(x = logFC, y = neg_log10_pval, color = EnzMod)) +
      geom_point() + .enzmod_scale + labs(color = "Mod") + theme_bw()
  )

  cat_legend <- get_legend(cat_plots[[1]])

  cat_plots <- lapply(cat_plots, function(p) p + theme(legend.position = "none"))
  mod_plots <- lapply(mod_plots, function(p) p + theme(legend.position = "none"))

  top_row <- plot_grid(plotlist = cat_plots, nrow = 1)
  bot_row <- plot_grid(plotlist = mod_plots, nrow = 1)

  top_row <- plot_grid(top_row, cat_legend, nrow = 1, rel_widths = c(1, 0.15))
  bot_row <- plot_grid(bot_row, mod_legend, nrow = 1, rel_widths = c(1, 0.15))

  row_label <- function(txt)
    ggdraw() + draw_label(txt, angle = 90, fontface = "bold", size = 11)

  row1 <- plot_grid(row_label("Modification Category"), top_row,
                    nrow = 1, rel_widths = c(0.06, 1))
  row2 <- plot_grid(row_label("Modification Type"),     bot_row,
                    nrow = 1, rel_widths = c(0.06, 1))

  title <- ggdraw() +
    draw_label(panel_title, fontface = "bold", size = 14)

  plot_grid(title, row1, row2, ncol = 1, rel_heights = c(0.04, 1, 1))
}


# ── Neurofilament panel (NEFL / NEFM / NEFH) ──────────────────────────────────

plot_neurofilament_panel <- function(table        = table_occ,
                                     eb_fit       = eb_fit_occ,
                                     use_pval     = FALSE,
                                     min_prop_obs = 0) {
  .plot_gene_panel(c("NEFL", "NEFM", "NEFH"), "Neurofilament Peptides",
                   table, eb_fit, use_pval, min_prop_obs)
}

p_nef <- plot_neurofilament_panel(use_pval = TRUE, min_prop_obs = 0.12)
ggsave("MainManuscript/Figures/NEF_peptide_panel.png", p_nef, width = 12.8, height = 8, dpi = 300)


# ── Translation elongation factor panel (EEF1G / EEF1B2 / EEF1D) ─────────────

plot_eef1_panel <- function(table        = table_occ,
                            eb_fit       = eb_fit_occ,
                            use_pval     = FALSE,
                            min_prop_obs = 0) {
  .plot_gene_panel(c("EEF1A1", "EEF1D", "EEF1G"), "Translation Elongation Factor Peptides",
                   table = table, eb_fit = eb_fit,
                   use_pval = use_pval, min_prop_obs = min_prop_obs)
}

p_eef1 <- plot_eef1_panel(use_pval = TRUE, min_prop_obs = 0.12)
ggsave("MainManuscript/Figures/EEF1_peptide_panel.png", p_eef1, width = 12.8, height = 8, dpi = 300)


# tables for protein specific peptides
