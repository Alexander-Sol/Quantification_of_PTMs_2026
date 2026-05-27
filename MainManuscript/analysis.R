# analysis.R — Main Manuscript
# Differential abundance and PTM analyses that feed Figures 1–4.
# Run this script top-to-bottom before sourcing plotting.R.
# Output objects are used directly by plotting.R; no intermediate files are written.
#
# Slow intermediate objects (y.peptide, dpcfit, y.protein, y.peptide.imputed) are cached as
# uncompressed .rds files in Data/. Delete the cache files to force a full rerun.

library(tidyverse)
library(limpa)
library(limma)

setwd("C:/Users/Alex/Source/Repos/Quantification_of_PTMs_2026")
source("MainManuscript/CustomScripts.R")

# ── Data Paths ────────────────────────────────────────────────────────────────
PATH_DIVERSE     <- "Data/DiversePtms/QuantifiedPeptides.tsv"
PATH_LIMITED     <- "Data/LimitedPtms/QuantifiedPeptides.tsv"
PATH_DIVERSE_PSM <- "Data/DiversePtms/AllPSMs.psmtsv"
PATH_DIVERSE_PEP <- "Data/DiversePtms/AllPeptides.psmtsv"
PATH_LIMITED_PSM <- "Data/LimitedPtms/AllPSMs.psmtsv"
PATH_LIMITED_PEP <- "Data/LimitedPtms/AllPeptides.psmtsv"

PROP_OBS_MIN_ANALYSIS <- 0.02  # >= 0.12 ≈ observed in at least 10 samples; applied throughout Figures 2 and 3
PROP_OBS_MIN_RESULTS <- 0.12
logFC_cutoff <- 1.5
pval_cutoff  <- 0.05

# Cache paths
RDS_DIV_PEPTIDE <- "Data/cache_diverse_peptide.rds"
RDS_DIV_DPCFIT  <- "Data/cache_diverse_dpcfit.rds"
RDS_DIV_PROTEIN <- "Data/cache_diverse_protein.rds"
RDS_Z           <- "Data/cache_z.rds"

RDS_LIM_PEPTIDE <- "Data/cache_limited_peptide.rds"
RDS_LIM_DPCFIT  <- "Data/cache_limited_dpcfit.rds"
RDS_LIM_PROTEIN <- "Data/cache_limited_protein.rds"


# ── Figure 1 | PSM & Peptide Discovery: Diverse vs Limited PTMs ──────────────
# Reads from local .psmtsv files (gitignored). Script will error here if those
# files have not been copied into Data/DiversePtms/ and Data/LimitedPtms/.

filter_mm_results <- function(df, q_threshold = 1, pep_q_threshold = 0.01) {
  df %>%
    filter(PEP_QValue <= pep_q_threshold,
           QValue     <= q_threshold,
           !grepl("C|D", Decoy.Contaminant.Target))
}

extract_file_id <- function(file_name) {
  basename(file_name) %>% str_split("_") %>% sapply(function(x) x[2])
}

# Load and filter PSMs and peptide-level results
diverse_psm <- read.csv(PATH_DIVERSE_PSM, sep = "\t", row.names = NULL) %>% filter_mm_results()
limited_psm <- read.csv(PATH_LIMITED_PSM, sep = "\t", row.names = NULL) %>% filter_mm_results()
diverse_pep <- read.csv(PATH_DIVERSE_PEP, sep = "\t", row.names = NULL) %>% filter_mm_results()
limited_pep <- read.csv(PATH_LIMITED_PEP, sep = "\t", row.names = NULL) %>% filter_mm_results()

# Per-cell PSM counts ----
count_psms <- function(psm_df, type_label) {
  psm_df %>%
    group_by(File.Name) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(Type = type_label, File = extract_file_id(File.Name))
}

psm_counts_combined <- bind_rows(
  count_psms(diverse_psm, "Search 2: Diverse PTMs"),
  count_psms(limited_psm, "Search 1: Limited PTMs")
)

file_order_psm <- count_psms(diverse_psm, "") %>% arrange(count) %>% pull(File)
psm_counts_combined$File <- factor(psm_counts_combined$File, levels = file_order_psm)
psm_counts_combined <- arrange(psm_counts_combined, desc(Type))

# Restrict PSMs to sequences validated at the peptide level
diverse_psm <- diverse_psm[diverse_psm$Full.Sequence %in% diverse_pep$Full.Sequence, ]
limited_psm <- limited_psm[limited_psm$Full.Sequence %in% limited_pep$Full.Sequence, ]

# Per-cell unique peptide counts ----
count_unique_peps <- function(psm_df, pep_df, type_label) {
  psm_df %>%
    filter(Full.Sequence %in% pep_df$Full.Sequence) %>%
    group_by(File.Name) %>%
    distinct(Full.Sequence, .keep_all = TRUE) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(Type = type_label, File = extract_file_id(File.Name))
}

pep_counts_combined <- bind_rows(
  count_unique_peps(diverse_psm, diverse_pep, "Search 2: Diverse PTMs"),
  count_unique_peps(limited_psm, limited_pep, "Search 1: Limited PTMs")
)

file_order_pep <- count_unique_peps(diverse_psm, diverse_pep, "") %>% arrange(count) %>% pull(File)
pep_counts_combined$File <- factor(pep_counts_combined$File, levels = file_order_pep)
pep_counts_combined <- arrange(pep_counts_combined, desc(Type))

# Total counts (used for inset bars and manuscript numbers)
# PSM total = all PSMs summed across cells; arranged largest-first so the
# larger bar is drawn behind the smaller one in the identity-position inset.
psm_totals <- psm_counts_combined %>%
  group_by(Type) %>%
  summarise(Total = sum(count), .groups = "drop") %>%
  arrange(desc(Total))
# Peptide total = globally unique sequences (a peptide seen in N cells counted once)
pep_totals <- tibble::tribble(
  ~Type,                       ~Total,
  "Search 2: Diverse PTMs",   n_distinct(diverse_pep$Full.Sequence),
  "Search 1: Limited PTMs",   n_distinct(limited_pep$Full.Sequence)
)

# Modification type breakdown (diverse PSMs only) ----
mod_list <- GetModsTextWReplacement(diverse_psm$Full.Sequence, all_mods = TRUE)

mod_counts_df <- mod_list %>%
  unlist() %>%
  str_split(., ", ") %>%
  unlist() %>%
  table() %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  setNames(c("all_mods", "Freq")) %>%
  arrange(desc(Freq)) %>%
  filter(all_mods != "") %>%
  mutate(
    Category = ifelse(all_mods %in% biological_mods, "Biological Mod", "Other Mod"),
    all_mods_renamed = all_mods %>%
      str_replace("^Carbamidomethyl$",   "Carbamido-\nmethylation") %>%
      str_replace("Carboxymethyllysine", "Carboxy-\nmethyllysine")  %>%
      str_replace("Carboxymethylation",  "Carboxy-\nmethylation")   %>%
      str_replace("^Carbamyl$",          "Carbamylation")           %>%
      str_replace("^Ammonia$",           "Ammonia\nLoss")           %>%
      str_replace("^Water$",             "Water\nLoss")             %>%
      str_replace("^Iron$",              "Iron\nAdduct")            %>%
      str_replace("^Sodium$",            "Sodium\nAdduct")          %>%
      str_replace("^Calcium$",           "Calcium\nAdduct")
  )

bio_mods_df <- mod_counts_df %>%
  filter(Category == "Biological Mod") %>%
  head(12) %>%
  mutate(all_mods_renamed = factor(all_mods_renamed, levels = all_mods_renamed[order(Freq)]))

other_mods_df <- mod_counts_df %>%
  filter(Category == "Other Mod") %>%
  head(10) %>%
  mutate(all_mods_renamed = factor(all_mods_renamed, levels = all_mods_renamed[order(Freq)]))


# ── Figure 1 & 4 | Diverse-PTM Search (fasta database) ───────────────────────

diverse_cache_files <- c(RDS_DIV_PEPTIDE, RDS_DIV_DPCFIT, RDS_DIV_PROTEIN, RDS_Z)

if (all(file.exists(diverse_cache_files))) {
  message("Loading diverse-PTM cache...")
  y.peptide <- readRDS(RDS_DIV_PEPTIDE)
  dpcfit    <- readRDS(RDS_DIV_DPCFIT)
  y.protein <- readRDS(RDS_DIV_PROTEIN)
  y.peptide.imputed <- readRDS(RDS_Z)
} else {
  message("Diverse-PTM cache not found — running full pipeline...")

  df <- read.csv(PATH_DIVERSE, sep = "\t", row.names = NULL)
  df <- df[grepl("Homo sapiens", df$Organism), ]
  df <- CleanAmbiguousPeptides(df)

  intensity_cols <- grep("^Intensity_", colnames(df), value = TRUE)
  expr           <- GetExpressionMatrix(df, intensity_cols)
  genes          <- GetGeneTable(df)

  col_keep <- colSums(!is.na(expr)) > 2500
  expr     <- expr[, col_keep]
  targets  <- GetMetaDataTdpStratified(expr)
  expr     <- expr[, colnames(expr) %in% targets$Sample]
  keep     <- rowSums(!is.na(expr)) >= 3
  expr     <- expr[keep, ]
  genes    <- genes[keep, ]
  stopifnot(all(targets$Sample == colnames(expr)))

  y.peptide <- structure(list(E = expr, genes = genes, targets = targets), class = "EList")
  dpcfit    <- dpc(y.peptide)
  y.protein <- dpcQuant(y.peptide, "ProteinGroup", dpc = dpcfit, verbose = TRUE)
  y.protein <- AddProteinGeneNames(y.protein, y.peptide)
  y.peptide.imputed <- dpcImpute(y.peptide, "PeptideSequence", dpc = dpcfit, verbose = TRUE)

  saveRDS(y.peptide, RDS_DIV_PEPTIDE, compress = FALSE)
  saveRDS(dpcfit,    RDS_DIV_DPCFIT,  compress = FALSE)
  saveRDS(y.protein, RDS_DIV_PROTEIN, compress = FALSE)
  saveRDS(y.peptide.imputed, RDS_Z, compress = FALSE)
  message("Diverse-PTM cache saved.")
}

targets <- y.peptide$targets


# ── Figure 2 & 3 | PTM Occupancy Analysis (PropObs >= 0.12) ──────────────────
design <- model.matrix(formula(~0 + Group), data = targets)
colnames(design) <- gsub("Group", "", colnames(design))
#colnames(design)[ncol(design)] <- "Sex"

# SE-propagating normalization helpers (override the versions in CustomScripts.R)
NormalizeExpressionToProteins <- function(peptide_e, protein_e, mod_seq, pg) {
  zero <- list(expr = rep(0, ncol(peptide_e$E)), se = NULL)

  modified_row  <- which(peptide_e$genes$PeptideSequence == mod_seq)
  modified_expr <- peptide_e$E[modified_row, ]
  modified_se   <- peptide_e$other$standard.error[modified_row, ]

  if (is.null(pg)) return(zero)

  protein_row <- which(protein_e$genes$ProteinGroup == pg)
  if (length(protein_row) == 0) return(zero)

  protein_expr <- protein_e$E[protein_row, ] %>% as.vector()
  protein_se   <- protein_e$other$standard.error[protein_row, ] %>% as.vector()

  stopifnot(all(peptide_e$targets$Sample == y.protein$targets$Sample))

  #normed_expr <- median(protein_expr) + (modified_expr - protein_expr)
  normed_expr <- (modified_expr - protein_expr)
  combined_se <- sqrt(modified_se^2 + protein_se^2)

  return(list(expr = normed_expr, se = combined_se))
}

GetNormalizedEList <- function(peptides_of_interest, y.peptide,
                               y.protein = NULL, normalize_to_protein = FALSE) {
  mod.peptide <- y.peptide[y.peptide$genes$PeptideSequence %in% peptides_of_interest, ]

  if (!normalize_to_protein) {
    raw <- apply(mod.peptide$genes, MARGIN = 1, FUN = function(x)
      NormalizeExpressionToPeptides(y.peptide,
                                    x[["PeptideSequence"]],
                                    x[["BaseSequence"]],
                                    x[["ProteinGroup"]]))
    normed_elist <- if (is.list(raw)) do.call(rbind, raw) else t(raw)
    se_list      <- NULL
  } else if (!is.null(y.protein)) {
    filtered_protein <- y.protein[y.protein$genes$NPeptides >= 4, ]
    results <- lapply(seq_len(nrow(mod.peptide$genes)), function(i) {
      x <- mod.peptide$genes[i, , drop = FALSE]
      NormalizeExpressionToProteins(y.peptide, filtered_protein,
                                    x[["PeptideSequence"]], x[["ProteinGroup"]])
    })
    normed_elist <- do.call(rbind, lapply(results, `[[`, "expr"))
    se_list      <- lapply(results, `[[`, "se")
  }

  keep <- !apply(normed_elist, 1, function(x) all(x == 0))

  normed.peptide   <- mod.peptide[keep, ]
  normed.peptide$E <- normed_elist[keep, , drop = FALSE]

  if (!is.null(se_list)) {
    se_kept <- se_list[keep]
    if (all(!sapply(se_kept, is.null)))
      normed.peptide$other$standard.error <- do.call(rbind, se_kept)
  }

  return(normed.peptide)
}

# Normalize each modified peptide to its parent protein and propagate SE.
# PropObs filter applied before fitting so low-coverage peptides are excluded
# from both the model and the multiple-testing correction.
all.mod.normed.pro <- GetNormalizedEList(
  y.peptide.imputed$genes$PeptideSequence,
  y.peptide.imputed,
  y.protein            = y.protein,
  normalize_to_protein = TRUE
)
all.mod.normed.pro <- all.mod.normed.pro[all.mod.normed.pro$genes$PropObs >= PROP_OBS_MIN_ANALYSIS, ]

fit_occ    <- dpcDE(all.mod.normed.pro, block = targets$Donor, design = design, plot = FALSE)
fit_occ2   <- contrasts.fit(fit_occ, makeContrasts(ALS_vs_CTRL = "ALS - CTRL", levels = design))
eb_fit_occ <- eBayes(fit_occ2)

table_occ <- topTable(eb_fit_occ, coef = "ALS_vs_CTRL", number = Inf)
table_occ$AllMod   <- GetModsTextWReplacement(table_occ$PeptideSequence, all_mods = TRUE) %>% sapply(paste)
table_occ$Category <- ClassifyPeptideModsSpecific(table_occ$AllMod)
table_occ          <- AddPositionColumns(table_occ)

base_filtered <- table_occ %>% filter(PropObs >= PROP_OBS_MIN_RESULTS)


# Figure 3: significant peptides (adj.P.Val <= 0.05 and |logFC| >= 1.5)
de_splice <- base_filtered %>%
  filter(adj.P.Val <= pval_cutoff, abs(logFC) >= logFC_cutoff, PropObs > 0.12) %>%
  arrange(desc(abs(logFC)))


# ── Figure 4 | Protein-Level DA: Diverse PTMs ────────────────────────────────
y.protein.filt <- y.protein[y.protein$genes$NPeptides >= 2, ]

fit_pro    <- dpcDE(y.protein.filt, design = design, block = targets$Donor,
                    plot = FALSE, sample.weights = TRUE)
fit_pro2   <- contrasts.fit(fit_pro, makeContrasts(ALS_vs_CTRL = "ALS - CTRL", levels = design))
eb_fit_pro <- eBayes(fit_pro2)

table_protein_diverse <- topTable(eb_fit_pro, coef = "ALS_vs_CTRL", number = Inf)

label_de <- function(tbl) {
  tbl %>% mutate(DEStatus = case_when(
    adj.P.Val < pval_cutoff & logFC >=  logFC_cutoff ~ "UP",
    adj.P.Val < pval_cutoff & logFC <= -logFC_cutoff ~ "DOWN",
    TRUE ~ "NotDE"
  ))
}
table_protein_diverse <- label_de(table_protein_diverse)

de_proteins <- table_protein_diverse %>%
  filter(DEStatus != "NotDE") %>%
  select(Gene, ProteinGroup, logFC, AveExpr, adj.P.Val, DEStatus, NPeptides) %>%
  arrange(desc(abs(logFC)))

write.table(de_proteins, "Data/DE_proteins.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
message(sprintf("Wrote %d DE proteins to Data/DE_proteins.tsv", nrow(de_proteins)))


# ── Figure 4 | Protein-Level DA: Limited PTMs ────────────────────────────────
limited_cache_files <- c(RDS_LIM_PEPTIDE, RDS_LIM_DPCFIT, RDS_LIM_PROTEIN)

if (all(file.exists(limited_cache_files))) {
  message("Loading limited-PTM cache...")
  y.peptide.lim <- readRDS(RDS_LIM_PEPTIDE)
  dpcfit_lim    <- readRDS(RDS_LIM_DPCFIT)
  y.protein.lim <- readRDS(RDS_LIM_PROTEIN)
} else {
  message("Limited-PTM cache not found — running full pipeline...")

  df_lim <- read.csv(PATH_LIMITED, sep = "\t", row.names = NULL)
  df_lim <- df_lim[grepl("Homo sapiens", df_lim$Organism), ]
  df_lim <- CleanAmbiguousPeptides(df_lim)

  intensity_cols_lim <- grep("^Intensity_", colnames(df_lim), value = TRUE)
  expr_lim  <- GetExpressionMatrix(df_lim, intensity_cols_lim)
  genes_lim <- GetGeneTable(df_lim, expr_lim)

  col_keep_lim <- colSums(!is.na(expr_lim)) > 2500
  expr_lim     <- expr_lim[, col_keep_lim]
  targets_lim  <- GetMetaDataTdpStratified(expr_lim)
  expr_lim     <- expr_lim[, colnames(expr_lim) %in% targets_lim$Sample]
  keep_lim     <- rowSums(!is.na(expr_lim)) >= 3
  expr_lim     <- expr_lim[keep_lim, ]
  genes_lim    <- genes_lim[keep_lim, ]
  stopifnot(all(targets_lim$Sample == colnames(expr_lim)))

  y.peptide.lim <- structure(list(E = expr_lim, genes = genes_lim, targets = targets_lim), class = "EList")
  dpcfit_lim    <- dpc(y.peptide.lim)
  y.protein.lim <- dpcQuant(y.peptide.lim, "ProteinGroup", dpc = dpcfit_lim, verbose = FALSE)
  y.protein.lim <- AddProteinGeneNames(y.protein.lim, y.peptide.lim)
  y.protein.lim <- y.protein.lim[y.protein.lim$genes$NPeptides >= 2, ]

  saveRDS(y.peptide.lim, RDS_LIM_PEPTIDE, compress = FALSE)
  saveRDS(dpcfit_lim,    RDS_LIM_DPCFIT,  compress = FALSE)
  saveRDS(y.protein.lim, RDS_LIM_PROTEIN, compress = FALSE)
  message("Limited-PTM cache saved.")
}

targets_lim <- y.peptide.lim$targets

design_lim <- model.matrix(formula(~0 + Group), data = targets_lim)
colnames(design_lim) <- gsub("Group", "", colnames(design_lim))
#colnames(design_lim)[ncol(design_lim)] <- "Sex"

fit_lim    <- dpcDE(y.protein.lim, design = design_lim, block = targets_lim$Donor,
                    plot = FALSE, sample.weights = TRUE)
fit_lim2   <- contrasts.fit(fit_lim, makeContrasts(ALS_vs_CTRL = "ALS - CTRL", levels = design_lim))
eb_fit_lim <- eBayes(fit_lim2)

table_protein_limited <- topTable(eb_fit_lim, coef = "ALS_vs_CTRL", number = Inf)
table_protein_limited <- label_de(table_protein_limited)

# Merged table used by Figure 4 comparison plots
diverse_de_genes <- table_protein_diverse$Gene[table_protein_diverse$DEStatus != "NotDE"]
limited_de_genes <- table_protein_limited$Gene[table_protein_limited$DEStatus != "NotDE"]

merged_pro <- merge(
  table_protein_diverse,
  table_protein_limited[, c("Gene", "logFC", "adj.P.Val")],
  by = "Gene", suffixes = c("", "_NoMod")
) %>%
  mutate(
    LogFC_Diff = logFC - logFC_NoMod,
    ModDiff = case_when(
      Gene %in% setdiff(limited_de_genes, diverse_de_genes) ~ "NoModDA",
      Gene %in% setdiff(diverse_de_genes, limited_de_genes) ~ "ModDA",
      TRUE ~ "Shared"
    )
  ) %>%
  arrange(desc(ModDiff))

mod_discordant <- merged_pro %>%
  filter(ModDiff != "Shared") %>%
  select(Gene, ProteinGroup, logFC, adj.P.Val, DEStatus,
         logFC_NoMod, adj.P.Val_NoMod, LogFC_Diff, ModDiff, NPeptides) %>%
  arrange(ModDiff, desc(abs(LogFC_Diff)))

write.table(mod_discordant, "Data/ModDiscordant_proteins.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
message(sprintf("Wrote %d mod-discordant proteins to Data/ModDiscordant_proteins.tsv", nrow(mod_discordant)))

discordant_category_summary <- function(mod_discordant, table_occ, min_peptides = 20) {
  target_proteins <- mod_discordant %>%
    filter(NPeptides > min_peptides) %>%
    pull(ProteinGroup)

  table_occ %>%
    filter(ProteinGroup %in% target_proteins) %>%
    group_by(Gene, Category) %>%
    summarise(NPeptides = n(), .groups = "drop") %>%
    pivot_wider(names_from = Category, values_from = NPeptides, values_fill = 0) %>%
    arrange(Gene)
}

discordant_cats <- discordant_category_summary(mod_discordant, table_occ)

mod_discordant$p.diff <- -1*log10(mod_discordant$adj.P.Val) - -1*log10(mod_discordant$adj.P.Val_NoMod)
test <- merge(mod_discordant, discordant_cats, by = "Gene")
write.table(test, "Data/ModDiscordant_withCat_breakdown.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


protein_cats <- discordant_category_summary(table_protein_diverse, table_occ, min_peptides = 20)
protein_cats$pct_mod <-  (protein_cats$`Enzymatic Mod` + protein_cats$`Non-Enzymatic Mod` + protein_cats$Carboxymethylation) /
 (protein_cats$Unmodified + protein_cats$Carbamidomethylation +protein_cats$`Enzymatic Mod` + protein_cats$`Non-Enzymatic Mod` + protein_cats$Carboxymethylation)
