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

# Modified peptide count vs logFC difference ----
# Classify all peptides in y.peptide$genes by mod category (not table_occ, which
# is restricted to proteins with ≥4 peptides). A new "Limited-Search Mod" category
# captures peptides whose only modifications are Oxidation on M or Deamidation on
# N/Q (detectable in a limited-PTM search); these are excluded from n_mod_strict.
# IsLimitedSearchPeptide checks raw mod strings from GetMods() before text
# replacement, which is required because .replaceMods collapses residue info
# (e.g. "Hydroxylation on N" and "Oxidation on M" both become "Oxidation").

pep_classified <- y.peptide$genes %>%
  mutate(
    AllMod   = GetModsTextWReplacement(PeptideSequence, all_mods = TRUE) %>% sapply(paste),
    Category = ClassifyPeptideModsSpecific(AllMod),
    Category = ifelse(
      Category == "Non-Enzymatic Mod" &
        IsLimitedSearchPeptide(GetMods(PeptideSequence)),
      "Limited-Search Mod",
      Category
    )
  )

n_mod_strict <- pep_classified %>%
  filter(Category %in% c("Enzymatic Mod", "Carboxymethylation", "Non-Enzymatic Mod")) %>%
  group_by(Gene) %>%
  summarise(n_mod = n(), .groups = "drop")
nrow(n_mod_strict)
sum(n_mod_strict$n_mod >= 2)

pc_for_plot <- protein_cats %>%
  inner_join(n_mod_strict, by = "Gene") %>%
  inner_join(
    merged_pro[, c("Gene", "LogFC_Diff", "ModDiff", "adj.P.Val", "adj.P.Val_NoMod", "logFC_NoMod")],
    by = "Gene"
  ) %>%
  filter(n_mod > 5, n_mod < 2000) %>%
  mutate(p_diff = -log10(adj.P.Val) - (-log10(adj.P.Val_NoMod)))

highlight_nmod <- pc_for_plot %>% filter(Gene %in% genes_of_interest_comp)

shared_nmod_aes <- list(
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linewidth = 0.7, linetype = "dashed"),
  scale_color_manual(values = modiff_colors),
  scale_fill_manual(values = modiff_colors, guide = "none"),
  geom_point(data = highlight_nmod, aes(fill = ModDiff),
             size = 4, shape = 21, color = "black", stroke = 0.8),
  geom_text_repel(data = highlight_nmod, aes(label = Gene),
                  color = "black", size = 5, fontface = "bold",
                  box.padding = 0.5, min.segment.length = 0,
                  max.overlaps = Inf, force = 3),
  labs(x = "Modified Peptides Detected"),
  theme_minimal(),
  theme(
    plot.title      = element_text(size = 14 * TEXT_SIZE),
    axis.title      = element_text(size = 12 * TEXT_SIZE),
    axis.text       = element_text(size = 11 * TEXT_SIZE),
    legend.position = "none"
  )
)

p_fig4_nmod_signed <- ggplot(pc_for_plot, aes(x = n_mod, y = LogFC_Diff, color = ModDiff)) +
  geom_point(alpha = 0.7, size = 2.5) +
  shared_nmod_aes +
  geom_hline(yintercept = 0, color = "grey20", linewidth = 0.4) +
  labs(y = "Log2FC Difference (Diverse − Limited)",
       title = "Modified Peptide Detection\nvs. Log2FC Discordance")

p_fig4_nmod_abs <- ggplot(pc_for_plot, aes(x = n_mod, y = abs(LogFC_Diff), color = ModDiff)) +
  geom_point(alpha = 0.7, size = 2.5) +
  shared_nmod_aes +
  labs(y = "|Log2FC Difference| (Diverse − Limited)",
       title = "Modified Peptide Detection\nvs. Log2FC Discordance (Magnitude)")

p_fig4_nmod_combined <- plot_grid(p_fig4_nmod_signed, p_fig4_nmod_abs,
                                  nrow = 1, labels = c("A", "B"))

ggsave(file.path(FIG_DIR, "Fig4_NMod_vs_LogFCDiff_Signed.png"),   p_fig4_nmod_signed,   width = 6, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_NMod_vs_LogFCDiff_Abs.png"),      p_fig4_nmod_abs,      width = 6, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_NMod_vs_LogFCDiff_Combined.png"), p_fig4_nmod_combined, width = 12, height = 5, dpi = 300)

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
  labs(x = "Log2FC Difference (Diverse − Limited)",
       y = "-Log10 Adj. P-value Difference (Diverse - Limited)",
       title = "Changes in Differential Abundance\nwhen Modifications are Considered") +
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
                    main       = "DE Proteins: Diverse vs. Limited PTMs\n(All Proteins)")

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
                        main       = "DA Proteins: Diverse vs. Limited PTMs\n(≥2 Modified Peptides)")

png(file.path(FIG_DIR, "Fig4_Venn_2ModPep.png"), width = 5, height = 5, units = "in", res = 300)
print(p_fig4_venn_mod)
dev.off()


# ── Figure 4 | Protein Intensity: Diverse vs Limited PTMs ────────────────────
get_mean_int <- function(elist, samples = NULL) {
  E <- elist$E
  if (!is.null(samples)) E <- E[, colnames(E) %in% samples, drop = FALSE]
  data.frame(ProteinGroup = elist$genes$ProteinGroup,
             Gene         = elist$genes$Gene,
             MeanInt      = rowMeans(E, na.rm = TRUE),
             stringsAsFactors = FALSE)
}

build_int_df <- function(samples_div, samples_lim) {
  inner_join(
    get_mean_int(y.protein.filt, samples_div) %>% select(ProteinGroup, Gene, Int_Diverse = MeanInt),
    get_mean_int(y.protein.lim,  samples_lim)  %>% select(ProteinGroup,       Int_Limited = MeanInt),
    by = "ProteinGroup"
  ) %>%
    filter(is.finite(Int_Diverse), is.finite(Int_Limited)) %>%
    left_join(n_mod_strict[, c("Gene", "n_mod")], by = "Gene")
}

als_samples_div  <- y.protein.filt$targets$Sample[y.protein.filt$targets$Group == "ALS"]
als_samples_lim  <- y.protein.lim$targets$Sample[y.protein.lim$targets$Group  == "ALS"]
ctrl_samples_div <- y.protein.filt$targets$Sample[y.protein.filt$targets$Group == "CTRL"]
ctrl_samples_lim <- y.protein.lim$targets$Sample[y.protein.lim$targets$Group  == "CTRL"]

int_all  <- build_int_df(NULL,            NULL)
int_als  <- build_int_df(als_samples_div,  als_samples_lim)
int_ctrl <- build_int_df(ctrl_samples_div, ctrl_samples_lim)

int_nmod_scale <- scale_color_viridis_c(trans = "log10", name = "Modified\nPeptides",
                                        limits = c(NA, 100), oob = scales::squish,
                                        breaks = c(10, 30, 100),
                                        labels = c("10", "30", "≥100"),
                                        na.value = "grey70")

int_scatter_theme <- list(
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black"),
  int_nmod_scale,
  theme_minimal(),
  theme(
    plot.title   = element_text(size = 14 * TEXT_SIZE),
    axis.title   = element_text(size = 12 * TEXT_SIZE),
    axis.text    = element_text(size = 11 * TEXT_SIZE),
    legend.text  = element_text(size = 10 * TEXT_SIZE),
    legend.title = element_text(size = 11 * TEXT_SIZE)
  )
)

make_int_plot <- function(df, title) {
  lims <- range(c(df$Int_Diverse, df$Int_Limited))
  ggplot(df, aes(x = Int_Limited, y = Int_Diverse, color = n_mod)) +
    geom_point(alpha = 0.6, size = 1.8) +
    int_scatter_theme +
    coord_equal(xlim = lims, ylim = lims) +
    labs(x = "Mean Log2 Intensity (Limited PTMs)",
         y = "Mean Log2 Intensity (Diverse PTMs)",
         title = title)
}

p_fig4_int_all  <- make_int_plot(int_all,  "Protein Intensity: All Samples")
p_fig4_int_als  <- make_int_plot(int_als,  "Protein Intensity: ALS Cells")
p_fig4_int_ctrl <- make_int_plot(int_ctrl, "Protein Intensity: CTRL Cells")

p_fig4_int_combined <- plot_grid(p_fig4_int_all, p_fig4_int_als, p_fig4_int_ctrl,
                                 nrow = 1, labels = c("A", "B", "C"))

ggsave(file.path(FIG_DIR, "Fig4_Intensity_All.png"),      p_fig4_int_all,      width = 5.5, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_Intensity_ALS.png"),      p_fig4_int_als,      width = 5.5, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_Intensity_CTRL.png"),     p_fig4_int_ctrl,     width = 5.5, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_Intensity_Combined.png"), p_fig4_int_combined, width = 15,  height = 5, dpi = 300)


# ── Figure 4 | Protein Intensity from QuantifiedProteins.tsv ─────────────────
PATH_DIVERSE_PRO <- "Data/DiversePtms/QuantifiedProteins.tsv"
PATH_LIMITED_PRO <- "Data/LimitedPtms/QuantifiedProteins.tsv"

# FlashLFQ names protein TSV columns Intensity_Default_N; rename them to the
# actual sample names by matching position with the peptide TSV header.
read_protein_tsv <- function(protein_path, peptide_path) {
  pro          <- read.csv(protein_path, sep = "\t", row.names = NULL, check.names = FALSE)
  pep_header   <- colnames(read.csv(peptide_path, sep = "\t", nrows = 0))
  pep_int_cols <- grep("^Intensity_", pep_header, value = TRUE)
  pro_int_idx  <- grep("^Intensity_Default_", colnames(pro))
  colnames(pro)[pro_int_idx] <- pep_int_cols
  pro %>%
    rename(ProteinGroup = `Protein Groups`, Gene = `Gene Name`) %>%
    filter(grepl("Homo sapiens", Organism))
}

mean_int_from_tsv <- function(pro_df, sample_names) {
  int_cols <- intersect(sample_names, colnames(pro_df))
  expr     <- apply(pro_df[, int_cols, drop = FALSE], 2, as.numeric)
  expr[expr == 0] <- NA
  expr <- log2(expr)
  data.frame(ProteinGroup = pro_df$ProteinGroup,
             Gene         = pro_df$Gene,
             MeanInt      = rowMeans(expr, na.rm = TRUE),
             stringsAsFactors = FALSE)
}

pro_tsv_div <- read_protein_tsv(PATH_DIVERSE_PRO, PATH_DIVERSE)
pro_tsv_lim <- read_protein_tsv(PATH_LIMITED_PRO, PATH_LIMITED)

build_int_tsv_df <- function(samples_div, samples_lim) {
  inner_join(
    mean_int_from_tsv(pro_tsv_div, samples_div) %>% select(ProteinGroup, Gene, Int_Diverse = MeanInt),
    mean_int_from_tsv(pro_tsv_lim, samples_lim)  %>% select(ProteinGroup,       Int_Limited = MeanInt),
    by = "ProteinGroup"
  ) %>%
    filter(is.finite(Int_Diverse), is.finite(Int_Limited)) %>%
    left_join(n_mod_strict[, c("Gene", "n_mod")], by = "Gene")
}

targets_div <- y.protein.filt$targets
targets_lim <- y.protein.lim$targets

int_tsv_all  <- build_int_tsv_df(targets_div$Sample,
                                 targets_lim$Sample)
int_tsv_als  <- build_int_tsv_df(targets_div$Sample[targets_div$Group == "ALS"],
                                 targets_lim$Sample[targets_lim$Group == "ALS"])
int_tsv_ctrl <- build_int_tsv_df(targets_div$Sample[targets_div$Group == "CTRL"],
                                 targets_lim$Sample[targets_lim$Group == "CTRL"])

p_fig4_tsv_all  <- make_int_plot(int_tsv_all,  "Protein Intensity (TSV): All Samples")
p_fig4_tsv_als  <- make_int_plot(int_tsv_als,  "Protein Intensity (TSV): ALS Cells")
p_fig4_tsv_ctrl <- make_int_plot(int_tsv_ctrl, "Protein Intensity (TSV): CTRL Cells")

p_fig4_tsv_combined <- plot_grid(p_fig4_tsv_all, p_fig4_tsv_als, p_fig4_tsv_ctrl,
                                 nrow = 1, labels = c("A", "B", "C"))

ggsave(file.path(FIG_DIR, "Fig4_IntensityTSV_All.png"),      p_fig4_tsv_all,      width = 5.5, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_IntensityTSV_ALS.png"),      p_fig4_tsv_als,      width = 5.5, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_IntensityTSV_CTRL.png"),     p_fig4_tsv_ctrl,     width = 5.5, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_IntensityTSV_Combined.png"), p_fig4_tsv_combined, width = 15,  height = 5, dpi = 300)





# -----------
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

# Diverse PTM volcano ----
highlight_pro <- table_protein_diverse %>% filter(Gene %in% genes_of_interest_pro)

p_fig4_pro_volcano <- ggplot(table_protein_diverse, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = DEStatus), alpha = 0.6, size = 2, shape = 19) +
  scale_color_manual(values = de_status_colors) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_point(data = highlight_pro,
             fill = highlight_fill_pro, color = "black",
             size = 4, shape = 21, stroke = 0.8) +
  geom_text_repel(data = highlight_pro, aes(label = Gene),
                  color = "black", size = 5, fontface = "bold",
                  box.padding = 0.5, min.segment.length = 0) +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value",
       title = "Differentially Abundant Proteins") +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 14 * TEXT_SIZE),
    axis.title      = element_text(size = 12 * TEXT_SIZE),
    axis.text       = element_text(size = 11 * TEXT_SIZE),
    legend.position = "none"
  )

ggsave(file.path(FIG_DIR, "Fig4_ProteinVolcano_Diverse.png"), p_fig4_pro_volcano, width = 6, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_ProteinVolcano_Diverse.svg"), p_fig4_pro_volcano, width = 6, height = 5)

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

# Adjusted p-value scatter ----
p_fig4_pval_scatter <- ggplot(merged_pro,
                              aes(x = -log10(adj.P.Val), y = -log10(adj.P.Val_NoMod))) +
  geom_point(aes(color = ModDiff), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = modiff_colors) +
  scale_fill_manual(values = modiff_colors, guide = "none") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = -log10(pval_cutoff), linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point(data = highlight_comp, aes(fill = ModDiff),
             size = 4, shape = 21, color = "black", stroke = 0.8) +
  geom_text_repel(data = highlight_comp, aes(label = Gene),
                  color = "black", size = 5, fontface = "bold",
                  box.padding = 1.5, point.padding = 0.5,
                  min.segment.length = 0, segment.color = "black",
                  max.segment.length = Inf, max.overlaps = Inf,
                  force = 3) +
  labs(x = "-Log10 Adjusted P-value (Diverse PTMs)",
       y = "-Log10 Adjusted P-value (Limited PTMs)",
       title = "Comparison of Protein Significance With\nand Without Modified Peptides") +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 14 * TEXT_SIZE),
    axis.title      = element_text(size = 12 * TEXT_SIZE),
    axis.text       = element_text(size = 11 * TEXT_SIZE),
    legend.position = "none"
  )

ggsave(file.path(FIG_DIR, "Fig4_PVal_Scatter.png"), p_fig4_pval_scatter, width = 6, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_PVal_Scatter.svg"), p_fig4_pval_scatter, width = 6, height = 5)

# Modified peptide count vs logFC difference ----
# Classify all peptides in y.peptide$genes by mod category (not table_occ, which
# is restricted to proteins with ≥4 peptides). A new "Limited-Search Mod" category
# captures peptides whose only modifications are Oxidation on M or Deamidation on
# N/Q (detectable in a limited-PTM search); these are excluded from n_mod_strict.
# IsLimitedSearchPeptide checks raw mod strings from GetMods() before text
# replacement, which is required because .replaceMods collapses residue info
# (e.g. "Hydroxylation on N" and "Oxidation on M" both become "Oxidation").

pep_classified <- y.peptide$genes %>%
  mutate(
    AllMod   = GetModsTextWReplacement(PeptideSequence, all_mods = TRUE) %>% sapply(paste),
    Category = ClassifyPeptideModsSpecific(AllMod),
    Category = ifelse(
      Category == "Non-Enzymatic Mod" &
        IsLimitedSearchPeptide(GetMods(PeptideSequence)),
      "Limited-Search Mod",
      Category
    )
  )

n_mod_strict <- pep_classified %>%
  filter(Category %in% c("Enzymatic Mod", "Carboxymethylation", "Non-Enzymatic Mod")) %>%
  group_by(Gene) %>%
  summarise(n_mod = n(), .groups = "drop")

pc_for_plot <- protein_cats %>%
  inner_join(n_mod_strict, by = "Gene") %>%
  inner_join(
    merged_pro[, c("Gene", "LogFC_Diff", "ModDiff", "adj.P.Val", "adj.P.Val_NoMod", "logFC_NoMod")],
    by = "Gene"
  ) %>%
  filter(n_mod > 5, n_mod < 2000) %>%
  mutate(p_diff = -log10(adj.P.Val) - (-log10(adj.P.Val_NoMod)))

highlight_nmod <- pc_for_plot %>% filter(Gene %in% genes_of_interest_comp)

shared_nmod_aes <- list(
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linewidth = 0.7, linetype = "dashed"),
  scale_color_manual(values = modiff_colors),
  scale_fill_manual(values = modiff_colors, guide = "none"),
  geom_point(data = highlight_nmod, aes(fill = ModDiff),
             size = 4, shape = 21, color = "black", stroke = 0.8),
  geom_text_repel(data = highlight_nmod, aes(label = Gene),
                  color = "black", size = 5, fontface = "bold",
                  box.padding = 0.5, min.segment.length = 0,
                  max.overlaps = Inf, force = 3),
  labs(x = "Modified Peptides Detected"),
  theme_minimal(),
  theme(
    plot.title      = element_text(size = 14 * TEXT_SIZE),
    axis.title      = element_text(size = 12 * TEXT_SIZE),
    axis.text       = element_text(size = 11 * TEXT_SIZE),
    legend.position = "none"
  )
)

p_fig4_nmod_signed <- ggplot(pc_for_plot, aes(x = n_mod, y = LogFC_Diff, color = ModDiff)) +
  geom_point(alpha = 0.7, size = 2.5) +
  shared_nmod_aes +
  geom_hline(yintercept = 0, color = "grey20", linewidth = 0.4) +
  labs(y = "Log2FC Difference (Diverse − Limited)",
       title = "Modified Peptide Detection\nvs. Log2FC Discordance")

p_fig4_nmod_abs <- ggplot(pc_for_plot, aes(x = n_mod, y = abs(LogFC_Diff), color = ModDiff)) +
  geom_point(alpha = 0.7, size = 2.5) +
  shared_nmod_aes +
  labs(y = "|Log2FC Difference| (Diverse − Limited)",
       title = "Modified Peptide Detection\nvs. Log2FC Discordance (Magnitude)")

p_fig4_nmod_combined <- plot_grid(p_fig4_nmod_signed, p_fig4_nmod_abs,
                                  nrow = 1, labels = c("A", "B"))

ggsave(file.path(FIG_DIR, "Fig4_NMod_vs_LogFCDiff_Signed.png"),   p_fig4_nmod_signed,   width = 6, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_NMod_vs_LogFCDiff_Abs.png"),      p_fig4_nmod_abs,      width = 6, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_NMod_vs_LogFCDiff_Combined.png"), p_fig4_nmod_combined, width = 12, height = 5, dpi = 300)

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
  labs(x = "Log2FC Difference (Diverse − Limited)",
       y = "Δ (−log10 Adj. P-value)",
       title = "Modified Peptide Count\nvs. DA Discordance") +
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
                    main       = "DE Proteins: Diverse vs. Limited PTMs")

png(file.path(FIG_DIR, "Fig4_Venn.png"), width = 5, height = 4, units = "in", res = 300)
print(p_fig4_venn)
dev.off()

# Venn diagram — proteins with >= 2 modified peptides ----
genes_with_2mods <- n_mod_strict %>% filter(n_mod >= 2) %>% pull(Gene)

fit_euler_mod <- euler(list(
  "Diverse PTMs" = intersect(diverse_de_genes, genes_with_2mods),
  "Limited PTMs" = intersect(limited_de_genes, genes_with_2mods)
))

p_fig4_venn_mod <- plot(fit_euler_mod,
                        quantities = list(cex = 1.3),
                        labels     = list(cex = 1.1, box = list(col = NA, fill = alpha("white", 0.6))),
                        fills      = list(fill = c("#ef5350", "#5c6bc0"), alpha = 0.45),
                        edges      = list(col = "grey20", lwd = 1.5),
                        main       = "DE Proteins: Diverse vs. Limited PTMs\n(≥2 Modified Peptides)")

png(file.path(FIG_DIR, "Fig4_Venn_2ModPep.png"), width = 5, height = 4, units = "in", res = 300)
print(p_fig4_venn_mod)
dev.off()


# ── Figure 4 | Protein Intensity: Diverse vs Limited PTMs ────────────────────
get_mean_int <- function(elist, samples = NULL) {
  E <- elist$E
  if (!is.null(samples)) E <- E[, colnames(E) %in% samples, drop = FALSE]
  data.frame(ProteinGroup = elist$genes$ProteinGroup,
             Gene         = elist$genes$Gene,
             MeanInt      = rowMeans(E, na.rm = TRUE),
             stringsAsFactors = FALSE)
}

build_int_df <- function(samples_div, samples_lim) {
  inner_join(
    get_mean_int(y.protein.filt, samples_div) %>% select(ProteinGroup, Gene, Int_Diverse = MeanInt),
    get_mean_int(y.protein.lim,  samples_lim)  %>% select(ProteinGroup,       Int_Limited = MeanInt),
    by = "ProteinGroup"
  ) %>%
    filter(is.finite(Int_Diverse), is.finite(Int_Limited)) %>%
    left_join(n_mod_strict[, c("Gene", "n_mod")], by = "Gene")
}

als_samples_div  <- y.protein.filt$targets$Sample[y.protein.filt$targets$Group == "ALS"]
als_samples_lim  <- y.protein.lim$targets$Sample[y.protein.lim$targets$Group  == "ALS"]
ctrl_samples_div <- y.protein.filt$targets$Sample[y.protein.filt$targets$Group == "CTRL"]
ctrl_samples_lim <- y.protein.lim$targets$Sample[y.protein.lim$targets$Group  == "CTRL"]

int_all  <- build_int_df(NULL,            NULL)
int_als  <- build_int_df(als_samples_div,  als_samples_lim)
int_ctrl <- build_int_df(ctrl_samples_div, ctrl_samples_lim)

int_nmod_scale <- scale_color_viridis_c(trans = "log10", name = "Modified\nPeptides",
                                        limits = c(NA, 100), oob = scales::squish,
                                        breaks = c(10, 30, 100),
                                        labels = c("10", "30", "≥100"),
                                        na.value = "grey70")

int_scatter_theme <- list(
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black"),
  int_nmod_scale,
  theme_minimal(),
  theme(
    plot.title   = element_text(size = 14 * TEXT_SIZE),
    axis.title   = element_text(size = 12 * TEXT_SIZE),
    axis.text    = element_text(size = 11 * TEXT_SIZE),
    legend.text  = element_text(size = 10 * TEXT_SIZE),
    legend.title = element_text(size = 11 * TEXT_SIZE)
  )
)

make_int_plot <- function(df, title) {
  lims <- range(c(df$Int_Diverse, df$Int_Limited))
  ggplot(df, aes(x = Int_Limited, y = Int_Diverse, color = n_mod)) +
    geom_point(alpha = 0.6, size = 1.8) +
    int_scatter_theme +
    coord_equal(xlim = lims, ylim = lims) +
    labs(x = "Mean Log2 Intensity (Limited PTMs)",
         y = "Mean Log2 Intensity (Diverse PTMs)",
         title = title)
}

p_fig4_int_all  <- make_int_plot(int_all,  "Protein Intensity: All Samples")
p_fig4_int_als  <- make_int_plot(int_als,  "Protein Intensity: ALS Cells")
p_fig4_int_ctrl <- make_int_plot(int_ctrl, "Protein Intensity: CTRL Cells")

p_fig4_int_combined <- plot_grid(p_fig4_int_all, p_fig4_int_als, p_fig4_int_ctrl,
                                 nrow = 1, labels = c("A", "B", "C"))

ggsave(file.path(FIG_DIR, "Fig4_Intensity_All.png"),      p_fig4_int_all,      width = 5.5, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_Intensity_ALS.png"),      p_fig4_int_als,      width = 5.5, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_Intensity_CTRL.png"),     p_fig4_int_ctrl,     width = 5.5, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_Intensity_Combined.png"), p_fig4_int_combined, width = 15,  height = 5, dpi = 300)


# ── Figure 4 | Protein Intensity from QuantifiedProteins.tsv ─────────────────
PATH_DIVERSE_PRO <- "Data/DiversePtms/QuantifiedProteins.tsv"
PATH_LIMITED_PRO <- "Data/LimitedPtms/QuantifiedProteins.tsv"

# FlashLFQ names protein TSV columns Intensity_Default_N; rename them to the
# actual sample names by matching position with the peptide TSV header.
read_protein_tsv <- function(protein_path, peptide_path) {
  pro          <- read.csv(protein_path, sep = "\t", row.names = NULL, check.names = FALSE)
  pep_header   <- colnames(read.csv(peptide_path, sep = "\t", nrows = 0))
  pep_int_cols <- grep("^Intensity_", pep_header, value = TRUE)
  pro_int_idx  <- grep("^Intensity_Default_", colnames(pro))
  colnames(pro)[pro_int_idx] <- pep_int_cols
  pro %>%
    rename(ProteinGroup = `Protein Groups`, Gene = `Gene Name`) %>%
    filter(grepl("Homo sapiens", Organism))
}

mean_int_from_tsv <- function(pro_df, sample_names) {
  int_cols <- intersect(sample_names, colnames(pro_df))
  expr     <- apply(pro_df[, int_cols, drop = FALSE], 2, as.numeric)
  expr[expr == 0] <- NA
  expr <- log2(expr)
  data.frame(ProteinGroup = pro_df$ProteinGroup,
             Gene         = pro_df$Gene,
             MeanInt      = rowMeans(expr, na.rm = TRUE),
             stringsAsFactors = FALSE)
}

pro_tsv_div <- read_protein_tsv(PATH_DIVERSE_PRO, PATH_DIVERSE)
pro_tsv_lim <- read_protein_tsv(PATH_LIMITED_PRO, PATH_LIMITED)

build_int_tsv_df <- function(samples_div, samples_lim) {
  inner_join(
    mean_int_from_tsv(pro_tsv_div, samples_div) %>% select(ProteinGroup, Gene, Int_Diverse = MeanInt),
    mean_int_from_tsv(pro_tsv_lim, samples_lim)  %>% select(ProteinGroup,       Int_Limited = MeanInt),
    by = "ProteinGroup"
  ) %>%
    filter(is.finite(Int_Diverse), is.finite(Int_Limited)) %>%
    left_join(n_mod_strict[, c("Gene", "n_mod")], by = "Gene")
}

targets_div <- y.protein.filt$targets
targets_lim <- y.protein.lim$targets

int_tsv_all  <- build_int_tsv_df(targets_div$Sample,
                                 targets_lim$Sample)
int_tsv_als  <- build_int_tsv_df(targets_div$Sample[targets_div$Group == "ALS"],
                                 targets_lim$Sample[targets_lim$Group == "ALS"])
int_tsv_ctrl <- build_int_tsv_df(targets_div$Sample[targets_div$Group == "CTRL"],
                                 targets_lim$Sample[targets_lim$Group == "CTRL"])

p_fig4_tsv_all  <- make_int_plot(int_tsv_all,  "Protein Intensity (TSV): All Samples")
p_fig4_tsv_als  <- make_int_plot(int_tsv_als,  "Protein Intensity (TSV): ALS Cells")
p_fig4_tsv_ctrl <- make_int_plot(int_tsv_ctrl, "Protein Intensity (TSV): CTRL Cells")

p_fig4_tsv_combined <- plot_grid(p_fig4_tsv_all, p_fig4_tsv_als, p_fig4_tsv_ctrl,
                                 nrow = 1, labels = c("A", "B", "C"))

ggsave(file.path(FIG_DIR, "Fig4_IntensityTSV_All.png"),      p_fig4_tsv_all,      width = 5.5, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_IntensityTSV_ALS.png"),      p_fig4_tsv_als,      width = 5.5, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_IntensityTSV_CTRL.png"),     p_fig4_tsv_ctrl,     width = 5.5, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, "Fig4_IntensityTSV_Combined.png"), p_fig4_tsv_combined, width = 15,  height = 5, dpi = 300)
