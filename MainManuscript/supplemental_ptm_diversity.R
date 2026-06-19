# supplemental_ptm_diversity.R
# Exports the full PTM diversity table as a supplemental .xlsx data file.
#
# This is the unabridged companion to poster_ptm_diversity.R. Where the poster
# table drops categories and collapses the minor lysine acylations into a single
# "Other" row, this file keeps every modification as its own row and adds the
# UniProt-vs-GPTMD breakdown of where each PTM's PSMs were identified.
#
# This script is self-contained: every count is generated fresh from the raw
# search output, with no dependency on poster_ptm_diversity.R or any data it
# produces. It maintains its own caches (cache/supp_*.rds) so reruns are fast;
# those caches are rebuilt automatically whenever a source .psmtsv changes.
#   Data/DiversePtms/AllPSMs.psmtsv   (UniProt + GPTMD diverse search)
#   Data/LimitedPtms/AllPSMs.psmtsv   (limited search)
# The PTM group definitions (category, modification, residues, merged-from mod
# strings, notes) are defined inline below, so the script reads no external table.
#
# Filters: Homo sapiens, Decoy/Contaminant/Target == "T", PEP_QValue <= 0.01
#
# Output:
#   MainManuscript/SupplementaryData_PTM_Diversity.xlsx

# install.packages("writexl")  # run once if needed
library(tidyverse)
library(writexl)

# Paths are relative to the repository root — set the working directory there
# before running (opening Quantification_of_PTMs.Rproj in RStudio does this).

CACHE_DIR <- "MainManuscript/cache"
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

diverse_psm_tsv <- "Data/DiversePtms/AllPSMs.psmtsv"
limited_psm_tsv <- "Data/LimitedPtms/AllPSMs.psmtsv"
out_xlsx        <- "MainManuscript/SupplementaryData_PTM_Diversity.xlsx"

cache_valid <- function(cache_path, source_paths) {
  file.exists(cache_path) &&
    file.mtime(cache_path) > max(file.mtime(source_paths))
}

# ── Read diverse PSM file: Mods + Full Sequence + Base Sequence ───────────────
read_diverse_psm <- function(psm_path) {
  message(sprintf("Reading %s ...", basename(psm_path)))

  # Sniff header to find exact protein accession column name
  header   <- colnames(data.table::fread(psm_path, nrows = 0, sep = "\t"))
  acc_col  <- grep("^accession$|protein.*accession", header, ignore.case = TRUE, value = TRUE)[1]
  if (is.na(acc_col))
    stop(sprintf("No protein accession column found in %s.\nAvailable: %s",
                 basename(psm_path), paste(header, collapse = ", ")))

  needed_cols <- c("Organism Name", "Mods", "Full Sequence", "Base Sequence",
                   acc_col, "Decoy/Contaminant/Target", "PEP_QValue")

  if (requireNamespace("data.table", quietly = TRUE)) {
    psm <- data.table::fread(psm_path, sep = "\t", select = needed_cols,
                              data.table = FALSE, showProgress = FALSE)
  } else {
    psm <- readr::read_tsv(psm_path, col_select = all_of(needed_cols),
                            show_col_types = FALSE)
  }

  colnames(psm) <- make.names(colnames(psm))
  names(psm)[names(psm) == make.names(acc_col)] <- "Protein.Accession"

  psm <- psm[
    grepl("Homo sapiens", psm$Organism.Name, fixed = TRUE) &
    psm$Decoy.Contaminant.Target == "T" &
    !is.na(psm$PEP_QValue) & psm$PEP_QValue <= 0.01,
    c("Mods", "Full.Sequence", "Base.Sequence", "Protein.Accession")
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

# ── Read PSM data, caching this script's own parse of the raw files ───────────
div_cache <- file.path(CACHE_DIR, "supp_diverse_psm.rds")
lim_cache <- file.path(CACHE_DIR, "supp_limited_mods.rds")

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

# Unique proteins (by accession) among matching diverse PSMs
count_proteins <- function(accession_col, mods_col, mod_strings) {
  pat <- mod_pattern(mod_strings)
  if (is.null(pat)) return(0L)
  matching <- grepl(pat, mods_col, perl = TRUE) & !is.na(mods_col)
  length(unique(accession_col[matching]))
}

# Distinct base sequences among matching diverse PSMs
count_peptides <- function(baseseq_col, mods_col, mod_strings) {
  pat <- mod_pattern(mod_strings)
  if (is.null(pat)) return(0L)
  matching <- grepl(pat, mods_col, perl = TRUE) & !is.na(mods_col)
  length(unique(baseseq_col[matching]))
}

# ── PTM group definitions ─────────────────────────────────────────────────────
# Each row defines one modification group: the residues it was observed on, the
# exact mod strings merged into it (matched against the PSM "Mods" column,
# pipe-separated), and a note. Counts are computed fresh from the PSM data below;
# none are stored here.
ptm_groups <- tibble::tribble(
  ~Category, ~Modification, ~Residues_Observed, ~Merged_From, ~Notes,

  "Phosphorylation", "Phosphorylation", "S; T; Y; H",
  "Phosphoserine on S | Phosphorylation on S | Phosphothreonine on T | Phosphorylation on T | Phosphotyrosine on Y | Phosphorylation on Y | Phosphohistidine on H", "",

  "Methylation", "Monomethylation", "K; R; H; Q",
  "Methylation on K | N6-methyllysine on K | Methylation on R | Omega-N-methylarginine on R | Tele-methylhistidine on H | Pros-methylhistidine on H | N5-methylglutamine on Q", "",

  "Methylation", "Dimethylation", "K; R",
  "N6,N6-dimethyllysine on K | Dimethylation on R | Asymmetric dimethylarginine on R | Symmetric dimethylarginine on R | Dimethylated arginine on R",
  "ADMA and SDMA are biologically distinct but isobaric in standard LC-MS/MS",

  "Methylation", "Trimethylation", "K",
  "Trimethylation on K | N6,N6,N6-trimethyllysine on K", "",

  "Acetylation", "Lysine acetylation", "K",
  "Acetylation on K | N6-acetyllysine on K",
  "Modification of the epsilon-amino group of lysine (epsilon-amino group of K side chain)",

  "Acetylation", "N-terminal protein acetylation", "Protein N-terminus",
  "N-acetylalanine on A | N-acetylserine on S | N-acetylmethionine on M | N-acetylthreonine on T | N-acetylvaline on V | N-acetylglycine on G | N-acetylcysteine on C | N-acetylaspartate on D | N-acetylglutamate on E | Acetylation on X",
  "Modification of the protein N-terminal alpha-amino group; residue refers to identity of N-terminal residue after Met processing",

  "Lysine acylation", "Formylation", "K",
  "Formylation on K", "",

  "Lysine acylation", "Succinylation", "K",
  "Succinylation on K | N6-succinyllysine on K", "",

  "Lysine acylation", "Glutarylation", "K",
  "Glutarylation on K | N6-glutaryllysine on K", "",

  "Lysine acylation", "Malonylation", "K",
  "Malonylation on K | N6-malonyllysine on K", "",

  "Lysine acylation", "Crotonylation", "K",
  "Crotonylation on K | N6-crotonyllysine on K", "",

  "Lysine acylation", "2-Hydroxyisobutyrylation", "K",
  "Hydroxybutyrylation on K | N6-(2-hydroxyisobutyryl)lysine on K", "",

  "Lysine acylation", "Butyrylation", "K",
  "Butyrylation on K | N6-butyryllysine on K", "",

  "Lysine (Other)", "Ubiquitination (GG remnant)", "K",
  "GG (Ubiquitination Site) on K",
  "Diglycine tag on lysine after tryptic digestion of ubiquitinated proteins",

  "Lysine (Other)", "Carboxymethyllysine", "K",
  "Carboxymethylation on K",
  "Carboxymethylation of the lysine epsilon-amino group by iodoacetic acid",

  "Hydroxylation", "Hydroxylation", "P; N; K; H",
  "Hydroxylation on P | 4-hydroxyproline on P | Hydroxyproline on P | 3-hydroxyproline on P | Hydroxylation on N | (3S)-3-hydroxyasparagine on N | Hydroxylation on K | 5-hydroxylysine on K | (3S)-3-hydroxyhistidine on H", "",

  "Deamidation", "Deamidation", "N; Q",
  "Deamidation on N | Deamidated asparagine on N | Deamidation on Q | Deamidated glutamine on Q", "",

  "Citrullination", "Citrullination", "R",
  "Citrullination on R | Citrulline on R",
  "Enzymatic deamidation of arginine to citrulline by peptidylarginine deiminases (PADs)",

  "Nitrogen loss", "Ammonia loss", "N; C",
  "Ammonia loss on N | Ammonia loss on C", "",

  "Nitrogen loss", "Pyrrolidone carboxylic acid", "Q",
  "Pyrrolidone carboxylic acid on Q",
  "N-terminal glutamine cyclization to pyroglutamate",

  "Glycosylation", "HexNAc", "S; T; N",
  "HexNAc on S | HexNAc on T | HexNAc on Nxs | HexNAc on Nxt",
  "O-GlcNAc on S/T; N-linked HexNAc on N in NxS and NxT sequons",

  "Nitrosylation", "Nitrosylation (+NO)", "C; Y",
  "Nitrosylation on C | S-nitrosocysteine on C | Nitrosylation on Y", "",

  "Carboxylation", "Carboxylation (+CO2)", "D; E; K",
  "Carboxylation on D | Carboxylation on E | 4-carboxyglutamate on E | Carboxylation on K",
  "Vitamin K-dependent gamma-carboxylation on E (Gla residue); mechanism on D/K may differ",

  "Other biological", "5-glutamyl glycerylphosphorylethanolamine", "E",
  "5-glutamyl glycerylphosphorylethanolamine on E", "",

  "Other biological", "Diphthamide", "H",
  "Diphthamide on H",
  "Unique two-step modification of a conserved His in elongation factor EF-2; target of diphtheria toxin ADP-ribosyltransferase",

  "Reagent artifact", "Carbamidomethylation (IAA)", "D; E; H; K; Y; X; U",
  "Carbamidomethyl on D | Carbamidomethyl on E | Carbamidomethyl on H | Carbamidomethyl on K | Carbamidomethyl on Y | Carbamidomethyl on X | Carbamidomethyl on U",
  "Fixed alkylation of Cys excluded (common to both searches); counts reflect off-target labeling on non-Cys residues only",

  "Reagent artifact", "Carboxymethylation (IAc)", "Peptide N-terminus; C; W",
  "Carboxymethylation on C | Carboxymethylation on W | Carboxymethylation on X",
  "Carboxymethylation of Cys/Trp side chains and peptide N-terminus by iodoacetic acid",

  "Reagent artifact", "Carbamylation", "M; K; R; C; X",
  "Carbamyl on M | Carbamyl on K | Carbamyl on R | Carbamyl on C | Carbamyl on X",
  "Artifact from urea in denaturation buffer or in-solution cyanate",

  "In-source / oxidative artifact", "Oxidation", "M",
  "Oxidation on M | Methionine (R)-sulfoxide on M | Methionine sulfoxide on M",
  "Common sample-handling oxidation artifact; methionine sulfoxide. Includes the UniProt-annotated forms Methionine (R)-sulfoxide and Methionine sulfoxide",

  "In-source / oxidative artifact", "Water loss", "E",
  "Water Loss on E",
  "In-source neutral loss or dehydration artifact",

  "Metal adduct", "Sodium adduct", "D; E",
  "Sodium on D | Sodium on E", "Metal contamination adduct on acidic residues",

  "Metal adduct", "Potassium adduct", "D; E",
  "Potassium on D | Potassium on E", "Metal contamination adduct on acidic residues",

  "Metal adduct", "Fe(III) adduct", "D; E",
  "Fe[III] on D | Fe[III] on E", "Metal contamination adduct on acidic residues",

  "Metal adduct", "Fe(II) adduct", "D; E",
  "Fe[II] on D | Fe[II] on E", "Metal contamination adduct on acidic residues",

  "Metal adduct", "Calcium adduct", "D; E",
  "Calcium on D | Calcium on E", "Metal contamination adduct on acidic residues",

  "Metal adduct", "Magnesium adduct", "D; E",
  "Magnesium on D | Magnesium on E", "Metal contamination adduct on acidic residues",

  "Metal adduct", "Zinc adduct", "D; E",
  "Zinc on D | Zinc on E", "Metal contamination adduct on acidic residues",

  "Metal adduct", "Copper(I) adduct", "D; E",
  "Cu[I] on D | Cu[I] on E",
  "Metal contamination adduct on acidic residues; previously omitted as rare"
)

mod_lists <- strsplit(ptm_groups$Merged_From, " \\| ")

ptm_groups$Diverse_PSM_Count     <- sapply(mod_lists, count_psms,
                                           mods_col = diverse_mods)
ptm_groups$Limited_PSM_Count     <- sapply(mod_lists, count_psms,
                                           mods_col = limited_mods)
ptm_groups$UniProt_PSM_Count     <- sapply(mod_lists, count_uniprot_psms,
                                           mods_col    = diverse_mods,
                                           fullseq_col = diverse_psm_df$Full.Sequence)
ptm_groups$GPTMD_PSM_Count       <- ptm_groups$Diverse_PSM_Count -
                                    ptm_groups$UniProt_PSM_Count
ptm_groups$Diverse_Peptide_Count <- sapply(mod_lists, count_peptides,
                                           baseseq_col = diverse_psm_df$Base.Sequence,
                                           mods_col    = diverse_mods)
ptm_groups$Diverse_Protein_Count <- sapply(mod_lists, count_proteins,
                                           accession_col = diverse_psm_df$Protein.Accession,
                                           mods_col      = diverse_mods)

# ── Category ordering ──────────────────────────────────────────────────────────
# Mirror the poster's biological-then-artifact ordering, then append any category
# present in the data but not listed (so nothing is silently dropped).
category_order <- c(
  "Phosphorylation", "Methylation", "Acetylation",
  "Lysine acylation", "Lysine (Other)",
  "Citrullination", "Glycosylation",
  "Hydroxylation", "Deamidation", "Nitrogen loss", "Nitrosylation",
  "Carboxylation", "Other biological",
  "Reagent artifact", "In-source / oxidative artifact", "Metal adduct"
)
category_order <- c(category_order,
                    setdiff(unique(ptm_groups$Category), category_order))

# ── Assemble the supplemental table: every row, full count breakdown ───────────
supp_table <- ptm_groups %>%
  mutate(Category = factor(Category, levels = category_order)) %>%
  arrange(Category, desc(Diverse_PSM_Count)) %>%
  transmute(
    Category                       = as.character(Category),
    Modification,
    `Residues observed`            = Residues_Observed,
    Proteins                       = as.integer(Diverse_Protein_Count),
    Peptidoforms                   = as.integer(Diverse_Peptide_Count),
    `PSMs (diverse search)`        = as.integer(Diverse_PSM_Count),
    `PSMs identified from UniProt` = as.integer(UniProt_PSM_Count),
    `PSMs discovered by GPTMD`     = as.integer(GPTMD_PSM_Count),
    `PSMs (limited search)`        = as.integer(Limited_PSM_Count),
    `Source modifications merged`  = Merged_From,
    Notes
  )

write_xlsx(supp_table, out_xlsx)

message(sprintf("Saved: %s  (%d modifications across %d categories)",
                out_xlsx, nrow(supp_table), dplyr::n_distinct(supp_table$Category)))
message(sprintf("  Diverse PSMs: %s  |  UniProt: %s  |  GPTMD: %s",
                format(sum(supp_table$`PSMs (diverse search)`),        big.mark = ","),
                format(sum(supp_table$`PSMs identified from UniProt`), big.mark = ","),
                format(sum(supp_table$`PSMs discovered by GPTMD`),     big.mark = ",")))
