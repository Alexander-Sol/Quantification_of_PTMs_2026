# poster_ptm_diversity.R
# Adds PSM/peptide counts and UniProt/GPTMD breakdown to Data/PTM_groups_SI.tsv,
# then exports a formatted PowerPoint table for the poster.
#
# Filters: Homo sapiens, Decoy/Contaminant/Target == "T", PEP_QValue <= 0.01
#
# New TSV columns added:
#   Diverse_PSM_Count, Limited_PSM_Count  — PSMs matching each PTM group
#   UniProt_PSM_Count                     — subset of diverse PSMs from UniProt database
#   GPTMD_PSM_Count                       — subset from GPTMD / other sources
#   Diverse_Peptide_Count                 — distinct base sequences in diverse search
#
# Output:
#   Data/PTM_groups_SI.tsv          (updated in place)
#   MainManuscript/Figures/Poster_PTM_Table.pptx

# install.packages(c("officer", "flextable"))  # run once if needed
library(tidyverse)
library(officer)
library(flextable)

setwd("C:/Users/Alex/Source/Repos/Quantification_of_PTMs_2026")

FIG_DIR   <- "MainManuscript/Figures"
CACHE_DIR <- "MainManuscript/cache"
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

diverse_psm_tsv <- "Data/DiversePtms/AllPSMs.psmtsv"
limited_psm_tsv <- "Data/LimitedPtms/AllPSMs.psmtsv"
ptm_groups_tsv  <- "Data/PTM_groups_SI.tsv"

cache_valid <- function(cache_path, source_paths) {
  file.exists(cache_path) &&
    file.mtime(cache_path) > max(file.mtime(source_paths))
}

# ── Read diverse PSM file: Mods + Full Sequence + Base Sequence ───────────────
read_diverse_psm <- function(psm_path) {
  message(sprintf("Reading %s ...", basename(psm_path)))
  needed_cols <- c("Organism Name", "Mods", "Full Sequence", "Base Sequence",
                   "Decoy/Contaminant/Target", "PEP_QValue")

  if (requireNamespace("data.table", quietly = TRUE)) {
    psm <- data.table::fread(psm_path, sep = "\t", select = needed_cols,
                              data.table = FALSE, showProgress = FALSE)
  } else {
    psm <- readr::read_tsv(psm_path, col_select = any_of(needed_cols),
                            show_col_types = FALSE)
  }

  colnames(psm) <- make.names(colnames(psm))

  psm <- psm[
    grepl("Homo sapiens", psm$Organism.Name, fixed = TRUE) &
    psm$Decoy.Contaminant.Target == "T" &
    !is.na(psm$PEP_QValue) & psm$PEP_QValue <= 0.01,
    c("Mods", "Full.Sequence", "Base.Sequence")
  ]

  message(sprintf("  → %d qualifying PSMs", nrow(psm)))
  psm
}

# ── Read limited PSM file: Mods only ─────────────────────────────────────────
read_psm_mods <- function(psm_path) {
  message(sprintf("Reading %s ...", basename(psm_path)))
  needed_cols <- c("Organism Name", "Mods", "Decoy/Contaminant/Target", "PEP_QValue")

  if (requireNamespace("data.table", quietly = TRUE)) {
    psm <- data.table::fread(psm_path, sep = "\t", select = needed_cols,
                              data.table = FALSE, showProgress = FALSE)
  } else {
    psm <- readr::read_tsv(psm_path, col_select = any_of(needed_cols),
                            show_col_types = FALSE)
  }

  colnames(psm) <- make.names(colnames(psm))

  psm <- psm[
    grepl("Homo sapiens", psm$Organism.Name, fixed = TRUE) &
    psm$Decoy.Contaminant.Target == "T" &
    !is.na(psm$PEP_QValue) & psm$PEP_QValue <= 0.01,
  ]

  message(sprintf("  → %d qualifying PSMs", nrow(psm)))
  psm$Mods
}

# ── Caches ────────────────────────────────────────────────────────────────────
# diverse_psm_data.rds stores a 3-column data frame (Mods, Full.Sequence, Base.Sequence)
div_cache <- file.path(CACHE_DIR, "diverse_psm_data.rds")
lim_cache <- file.path(CACHE_DIR, "limited_psm_mods.rds")

if (cache_valid(div_cache, diverse_psm_tsv)) {
  message("Loading diverse PSM data from cache...")
  diverse_psm_df <- readRDS(div_cache)
} else {
  diverse_psm_df <- read_diverse_psm(diverse_psm_tsv)
  saveRDS(diverse_psm_df, div_cache)
}
diverse_mods <- diverse_psm_df$Mods

if (cache_valid(lim_cache, limited_psm_tsv)) {
  message("Loading limited PSM mods from cache...")
  limited_mods <- readRDS(lim_cache)
} else {
  limited_mods <- read_psm_mods(limited_psm_tsv)
  saveRDS(limited_mods, lim_cache)
}

# ── Counting functions ────────────────────────────────────────────────────────

# Build escaped regex pattern from Merged_From mod strings
mod_pattern <- function(mod_strings) {
  mod_strings <- mod_strings[nchar(mod_strings) > 0]
  if (length(mod_strings) == 0) return(NULL)
  escaped <- gsub("([.+*?^${}()|\\[\\]\\\\])", "\\\\\\1", mod_strings, perl = TRUE)
  paste(escaped, collapse = "|")
}

# PSMs where Mods column matches any string in this group
count_psms <- function(mods_col, mod_strings) {
  pat <- mod_pattern(mod_strings)
  if (is.null(pat)) return(0L)
  sum(grepl(pat, mods_col, perl = TRUE), na.rm = TRUE)
}

# Diverse PSMs where Full.Sequence contains [UniProt:matched_mod...]
count_uniprot_psms <- function(mods_col, fullseq_col, mod_strings) {
  pat <- mod_pattern(mod_strings)
  if (is.null(pat)) return(0L)
  matching     <- grepl(pat, mods_col, perl = TRUE) & !is.na(mods_col)
  uniprot_pat  <- paste0("\\[UniProt:(?:", pat, ")")
  sum(matching & grepl(uniprot_pat, fullseq_col, perl = TRUE), na.rm = TRUE)
}

# Distinct base sequences among matching diverse PSMs
count_peptides <- function(baseseq_col, mods_col, mod_strings) {
  pat <- mod_pattern(mod_strings)
  if (is.null(pat)) return(0L)
  matching <- grepl(pat, mods_col, perl = TRUE) & !is.na(mods_col)
  length(unique(baseseq_col[matching]))
}

# ── Load PTM groups, compute all columns, save ───────────────────────────────
ptm_groups <- read.csv(ptm_groups_tsv, sep = "\t", stringsAsFactors = FALSE,
                        check.names = TRUE)

mod_lists <- strsplit(ptm_groups$Merged_From, " \\| ")

ptm_groups$Diverse_PSM_Count    <- sapply(mod_lists, count_psms,
                                           mods_col = diverse_mods)
ptm_groups$Limited_PSM_Count    <- sapply(mod_lists, count_psms,
                                           mods_col = limited_mods)
ptm_groups$UniProt_PSM_Count    <- sapply(mod_lists, count_uniprot_psms,
                                           mods_col    = diverse_mods,
                                           fullseq_col = diverse_psm_df$Full.Sequence)
ptm_groups$GPTMD_PSM_Count      <- ptm_groups$Diverse_PSM_Count -
                                    ptm_groups$UniProt_PSM_Count
ptm_groups$Diverse_Peptide_Count <- sapply(mod_lists, count_peptides,
                                            baseseq_col = diverse_psm_df$Base.Sequence,
                                            mods_col    = diverse_mods)

write.table(ptm_groups, ptm_groups_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
message(sprintf("Updated %s", ptm_groups_tsv))
message(sprintf("  Diverse: %s PSMs, %s unique peptides",
                format(sum(ptm_groups$Diverse_PSM_Count),    big.mark = ","),
                format(sum(ptm_groups$Diverse_Peptide_Count), big.mark = ",")))
message(sprintf("  UniProt: %s PSMs  |  GPTMD: %s PSMs",
                format(sum(ptm_groups$UniProt_PSM_Count), big.mark = ","),
                format(sum(ptm_groups$GPTMD_PSM_Count),  big.mark = ",")))
message(sprintf("  Limited: %s PSMs",
                format(sum(ptm_groups$Limited_PSM_Count), big.mark = ",")))


# ══════════════════════════════════════════════════════════════════════════════
# POWERPOINT TABLE
# ══════════════════════════════════════════════════════════════════════════════

category_order <- c(
  "Phosphorylation", "Methylation", "Acetylation",
  "Lysine acylation", "Lysine (Other)",
  "Citrullination", "Glycosylation",
  "Hydroxylation", "Deamidation", "Nitrogen loss", "Nitrosylation",
  "Carboxylation", "Other biological",
  "Reagent artifact", "In-source / oxidative artifact", "Metal adduct"
)

table_data <- ptm_groups %>%
  mutate(Category = factor(Category, levels = category_order)) %>%
  arrange(Category, desc(Diverse_PSM_Count)) %>%
  transmute(
    Category,
    Modification,
    Residues       = Residues_Observed,
    Proteins       = as.integer(Total_Count),
    diverse_pep    = as.integer(Diverse_Peptide_Count),
    diverse_psms   = as.integer(Diverse_PSM_Count),
    uniprot_psms   = as.integer(UniProt_PSM_Count),
    limited_psms   = ifelse(
      Limited_PSM_Count == 0, "—",
      formatC(Limited_PSM_Count, format = "d", big.mark = ",")
    )
  )

# Alternating background per category block
cat_sequence <- levels(droplevels(table_data$Category))
cat_bg       <- setNames(
  rep(c("white", "#E1E5E7"), length.out = length(cat_sequence)),
  cat_sequence
)

# Build flextable ──────────────────────────────────────────────────────────────
ft <- flextable(table_data) %>%
  set_header_labels(
    diverse_pep  = "Diverse Peptides",
    diverse_psms = "Diverse PSMs",
    uniprot_psms = "UniProt PSMs",
    limited_psms = "Limited PSMs"
  ) %>%
  merge_v(j = "Category") %>%
  valign(j = "Category", valign = "top") %>%
  colformat_int(j = c("Proteins", "diverse_pep", "diverse_psms", "uniprot_psms"),
                big.mark = ",") %>%
  align(j = c("Proteins", "diverse_pep", "diverse_psms", "uniprot_psms", "limited_psms"),
        align = "right", part = "all") %>%
  align(j = c("Category", "Modification", "Residues"),
        align = "left", part = "all") %>%
  bold(part = "header") %>%
  bg(part = "header", bg = "#9B0000") %>%
  color(part = "header", color = "white") %>%
  bold(j = "Category") %>%
  # Column widths (total ~12.9")
  width(j = "Category",     width = 1.40) %>%
  width(j = "Modification", width = 1.70) %>%
  width(j = "Residues",     width = 1.70) %>%
  width(j = "Proteins",     width = 0.80) %>%
  width(j = "diverse_pep",  width = 0.90) %>%
  width(j = "diverse_psms", width = 0.90) %>%
  width(j = "uniprot_psms", width = 0.90) %>%
  width(j = "limited_psms", width = 0.85) %>%
  font(fontname = "Calibri", part = "all") %>%
  fontsize(size = 9,  part = "body") %>%
  fontsize(size = 10, part = "header") %>%
  height_all(height = 0.145) %>%
  padding(padding.top = 2, padding.bottom = 2,
          padding.left = 4, padding.right = 4, part = "all") %>%
  border_outer(border = fp_border(color = "#aaaaaa", width = 0.75)) %>%
  border_inner_h(border = fp_border(color = "#dddddd", width = 0.5)) %>%
  border_inner_v(border = fp_border(color = "#cccccc", width = 0.5))

# Alternating background by category block
for (cat_name in names(cat_bg)) {
  rows <- which(as.character(table_data$Category) == cat_name)
  if (length(rows) > 0 && cat_bg[[cat_name]] != "white")
    ft <- bg(ft, i = rows, bg = cat_bg[[cat_name]])
}

# Bold non-zero limited PSM values
lim_rows <- which(ptm_groups$Limited_PSM_Count[
  order(factor(ptm_groups$Category, levels = category_order),
        -ptm_groups$Diverse_PSM_Count)
] > 0)
if (length(lim_rows) > 0)
  ft <- bold(ft, i = lim_rows, j = "limited_psms")

prs <- read_pptx() %>%
  add_slide(layout = "Blank", master = "Office Theme") %>%
  ph_with(ft, location = ph_location(left = 0.2, top = 0.15,
                                      width = 12.9, height = 7.5))

print(prs, target = file.path(FIG_DIR, "Poster_PTM_Table.pptx"))
message("Saved: Poster_PTM_Table.pptx")
