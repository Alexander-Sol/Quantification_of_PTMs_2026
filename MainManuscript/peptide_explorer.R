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


# ── Shared data prep ───────────────────────────────────────────────────────────

.prep_pep_rows <- function(gene_name, table, eb_fit) {
  pep_rows <- table[table$Gene == gene_name, ]
  pep_rows <- pep_rows[!is.na(pep_rows$logFC), ]
  if (nrow(pep_rows) == 0) stop("No peptides found for: ", gene_name)
  pep_rows$s2.post <- eb_fit$s2.post[match(rownames(pep_rows), rownames(eb_fit))]
  pep_rows
}

.base_aes <- list(
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60"),
  labs(x = "log2 FC (ALS vs CTRL)", y = "Posterior variance (s2.post)"),
  theme_bw()
)


# ── Individual plots ───────────────────────────────────────────────────────────

plot_protein_peptides <- function(gene_name,
                                  table  = table_occ,
                                  eb_fit = eb_fit_occ) {
  pep_rows <- .prep_pep_rows(gene_name, table, eb_fit)

  ggplot(pep_rows, aes(x = logFC, y = s2.post,
                       color = PropObs, label = PeptideSequence)) +
    geom_point(size = 3, alpha = 0.85) +
    .base_aes +
    scale_color_viridis_c(option = "plasma", limits = c(0, 1)) +
    labs(title = "Proportion Observed", color = "PropObs")
}


plot_protein_peptides_mod <- function(gene_name,
                                      table  = table_occ,
                                      eb_fit = eb_fit_occ) {
  pep_rows <- .prep_pep_rows(gene_name, table, eb_fit)

  ggplot(pep_rows, aes(x = logFC, y = s2.post,
                       color = Category, label = PeptideSequence)) +
    geom_point(size = 3, alpha = 0.85) +
    .base_aes +
    scale_color_manual(values = mod_colors, na.value = "gray70") +
    labs(title = "Modification Category", color = "Category")
}


plot_protein_peptides_intensity <- function(gene_name,
                                            table  = table_occ,
                                            eb_fit = eb_fit_occ,
                                            y.pep  = y.peptide) {
  pep_rows <- .prep_pep_rows(gene_name, table, eb_fit)

  idx <- match(pep_rows$PeptideSequence, y.pep$genes$PeptideSequence)
  pep_rows$mean_intensity <- rowMeans(y.pep$E[idx, , drop = FALSE], na.rm = TRUE)

  ggplot(pep_rows, aes(x = logFC, y = s2.post,
                       color = mean_intensity, label = PeptideSequence)) +
    geom_point(size = 3, alpha = 0.85) +
    .base_aes +
    scale_color_viridis_c(option = "viridis", na.value = "gray70") +
    labs(title = "Mean Intensity", color = "Mean\nIntensity")
}


plot_protein_peptides_groups <- function(gene_name,
                                         table  = table_occ,
                                         eb_fit = eb_fit_occ) {
  pep_rows <- .prep_pep_rows(gene_name, table, eb_fit)

  # Build adjacency on BaseSequence containment
  seqs <- pep_rows$BaseSequence
  n    <- length(seqs)
  adj  <- matrix(FALSE, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j)
        adj[i, j] <- grepl(seqs[i], seqs[j], fixed = TRUE) |
                     grepl(seqs[j], seqs[i], fixed = TRUE)
    }
  }

  # Connected components (BFS)
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

  # Rank groups by size, cap at 6
  group_ranks <- names(sort(table(group_id), decreasing = TRUE))
  top6        <- as.integer(group_ranks[seq_len(min(6, length(group_ranks)))])
  has_other   <- length(group_ranks) > 6

  pep_rows$SeqGroup <- factor(
    ifelse(pep_rows$group_id %in% top6,
           as.character(match(pep_rows$group_id, top6)),
           "Other"),
    levels = c(if (has_other) "Other", as.character(seq_len(length(top6))))
  )

  # Draw Other first so numbered groups render on top
  pep_rows <- pep_rows[order(pep_rows$SeqGroup, decreasing = FALSE), ]

  n_top <- length(top6)
  pal   <- c(if (has_other) c("Other" = "gray70"),
             setNames(scales::hue_pal()(n_top), as.character(seq_len(n_top))))

  ggplot(pep_rows, aes(x = logFC, y = s2.post,
                       color = SeqGroup, label = PeptideSequence)) +
    geom_point(size = 3, alpha = 0.85) +
    .base_aes +
    scale_color_manual(values = pal) +
    labs(title = "Sequence Groups", color = "Sequence\nGroup")
}


# ── Combined grid ──────────────────────────────────────────────────────────────

plot_protein_all <- function(gene_name,
                             table  = table_occ,
                             eb_fit = eb_fit_occ,
                             y.pep  = y.peptide) {
  p1 <- plot_protein_peptides(gene_name, table, eb_fit)
  p2 <- plot_protein_peptides_mod(gene_name, table, eb_fit)
  p3 <- plot_protein_peptides_intensity(gene_name, table, eb_fit, y.pep)
  p4 <- plot_protein_peptides_groups(gene_name, table, eb_fit)

  title <- ggdraw() +
    draw_label(paste("Peptides for", gene_name), fontface = "bold", size = 14)

  grid <- plot_grid(p1, p2, p3, p4, nrow = 2, ncol = 2)
  plot_grid(title, grid, ncol = 1, rel_heights = c(0.05, 1))
}

plot_protein_all("NEFM")
