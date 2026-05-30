# spectrum_viewer.R
# Annotated MS2 spectrum viewer.
# Reads matched-ion annotations from AllPeptides.psmtsv and experimental
# peaks from the corresponding mzML file, then plots with ggplot2.
# Requires analysis.R to have been sourced (for CustomScripts.R path resolution).

library(Spectra)
library(tidyverse)

MZML_DIR <- "D:/Kelly_ALS_motor_nueron_dataset/TDP_Stratified/MM1p1p0_CalAvgGPTMDSearch/Task2-AveragingTask"
PSM_PATH <- "Data/DiversePtms/AllPeptides.psmtsv"

if (exists("diverse_pep")) {
  psm_table <- diverse_pep
} else {
  psm_table <- read.csv(PSM_PATH, sep = "\t", check.names = FALSE, row.names = NULL)
}
names(psm_table) <- make.names(names(psm_table))


# ── Helpers ───────────────────────────────────────────────────────────────────

# Parse "[b1+1:129.1, b2+1:244.1, ...];[y1+1:147.1, ...]" style fields.
# Returns a data frame with ion label, type, and numeric value columns.
.parse_ion_field <- function(value_str) {
  pairs <- value_str |>
    gsub(pattern = "[\\[\\]]", replacement = "", perl = TRUE) |>
    gsub(pattern = ";", replacement = ",") |>
    strsplit(",") |>
    unlist() |>
    trimws()
  pairs <- pairs[nchar(pairs) > 0]
  data.frame(
    ion      = sub(":.*", "", pairs),
    value    = as.numeric(sub(".*:", "", pairs)),
    stringsAsFactors = FALSE
  )
}

# Load experimental peaks for a single scan from an mzML file.
# MetaMorpheus scan numbers correspond to acquisitionNum in Spectra.
.load_peaks <- function(file_name, scan_num, mzml_dir = MZML_DIR) {
  path <- file.path(mzml_dir, paste0(file_name, ".mzML"))
  sps  <- Spectra(path)
  hit  <- which(sps$acquisitionNum == scan_num)
  if (length(hit) == 0) stop("Scan ", scan_num, " not found in ", basename(path))
  peaksData(sps[hit[1]])[[1]]
}


# ── Sequence ion diagram ──────────────────────────────────────────────────────

# Parse a MetaMorpheus Full.Sequence string into a per-residue data frame.
# Modifications follow their residue in brackets: "SM[Ox]R" → M at pos 2 has mod.
.parse_full_sequence <- function(full_seq) {
  chars    <- strsplit(full_seq, "")[[1]]
  residues <- character(0)
  mods     <- character(0)
  i <- 1
  while (i <= length(chars)) {
    if (chars[i] == "[") {
      j <- i + 1; depth <- 1
      while (j <= length(chars) && depth > 0) {
        if (chars[j] == "[") depth <- depth + 1
        if (chars[j] == "]") depth <- depth - 1
        j <- j + 1
      }
      mods[length(mods)] <- paste(chars[(i + 1):(j - 2)], collapse = "")
      i <- j
    } else {
      residues <- c(residues, chars[i])
      mods     <- c(mods, NA_character_)
      i <- i + 1
    }
  }
  data.frame(position = seq_along(residues), residue = residues, mod = mods,
             stringsAsFactors = FALSE)
}

# Extract a short display name from a raw mod string.
# "UniProt:Phosphoserine on S" → "Phosphoserine"
# "Common Variable:Oxidation on M" → "Oxidation"
.mod_label <- function(mod_str) {
  if (is.na(mod_str)) return(NA_character_)
  trimws(sub(" on [A-Z].*$", "", sub("^[^:]+:", "", mod_str)))
}

# Draws the peptide sequence with:
#   - colored circles behind modified residues (fill = modification type)
#   - corner marks at each detected b/y cut site
#   - b ions on bottom (hash left), y ions on top (hash right)
#   - modification legend along the bottom
.plot_sequence_ions <- function(full_seq, matched_ions) {
  pep_df           <- .parse_full_sequence(full_seq)
  n                <- nrow(pep_df)
  pep_df$mod_label <- vapply(pep_df$mod, .mod_label, character(1))

  mod_types <- unique(pep_df$mod_label[!is.na(pep_df$mod_label)])
  mod_pal   <- if (length(mod_types) > 0)
    setNames(scales::hue_pal()(length(mod_types)), mod_types) else character(0)

  mod_df <- pep_df[!is.na(pep_df$mod_label), ]

  ion_num <- function(ions)
    as.integer(regmatches(ions$ion, regexpr("\\d+", ions$ion)))

  b_ions <- matched_ions[matched_ions$ion_type == "b", ]
  y_ions <- matched_ions[matched_ions$ion_type == "y", ]

  make_corners <- function(cx, side) {
    if (length(cx) == 0)
      return(data.frame(x = numeric(0), xend = numeric(0),
                        y = numeric(0), yend = numeric(0)))
    tick <- 0.55; hash <- 0.35
    if (side == "bottom")
      rbind(data.frame(x = cx, xend = cx,        y = 0,     yend = -tick),
            data.frame(x = cx, xend = cx - hash, y = -tick, yend = -tick))
    else
      rbind(data.frame(x = cx, xend = cx,        y = 0,     yend =  tick),
            data.frame(x = cx, xend = cx + hash, y =  tick, yend =  tick))
  }

  b_cx   <- if (nrow(b_ions) > 0) ion_num(b_ions) + 0.5     else numeric(0)
  y_cx   <- if (nrow(y_ions) > 0) n - ion_num(y_ions) + 0.5 else numeric(0)
  b_segs <- make_corners(b_cx, "bottom")
  y_segs <- make_corners(y_cx, "top")

  # Build in a single chain so scale_fill_manual binds to its geom correctly.
  # shape 21 needs an explicit color (stroke) to render; "white" keeps it clean.
  ggplot() +
    geom_segment(aes(x = 0.5, xend = n + 0.5, y = 0, yend = 0),
                 color = "grey50", linewidth = 0.5) +
    (if (nrow(mod_df) > 0)
      geom_point(data = mod_df, aes(x = position, y = 0, fill = mod_label),
                 shape = 21, size = 9, color = "white", stroke = 0.4, alpha = 0.9)
     else NULL) +
    (if (length(mod_pal) > 0)
      scale_fill_manual(values = mod_pal, name = NULL,
                        guide = guide_legend(nrow = 1,
                                             override.aes = list(size = 4)))
     else NULL) +
    geom_text(data = pep_df, aes(x = position, y = 0, label = residue),
              size = 3.8, fontface = "bold", family = "mono") +
    geom_segment(data = b_segs, aes(x = x, xend = xend, y = y, yend = yend),
                 color = "#2166ac", linewidth = 1.2) +
    geom_segment(data = y_segs, aes(x = x, xend = xend, y = y, yend = yend),
                 color = "#d6604d", linewidth = 1.2) +
    annotate("text", x = 0.2, y = -0.45, label = "b",
             color = "#2166ac", size = 3.5, fontface = "bold", hjust = 1) +
    annotate("text", x = 0.2, y =  0.45, label = "y",
             color = "#d6604d", size = 3.5, fontface = "bold", hjust = 1) +
    scale_x_continuous(limits = c(-0.2, n + 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1, 1)) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.text     = element_text(size = 8),
          legend.key.size = unit(0.9, "lines"),
          legend.margin   = margin(t = 2, unit = "pt"))
}


# ── Main plotting function ────────────────────────────────────────────────────

#' Plot an annotated MS2 spectrum with sequence ion diagram for a peptide.
#'
#' @param full_sequence  Full sequence string as it appears in the psmtsv
#'                       "Full Sequence" column (e.g. with mod brackets).
#' @param psm            PSM table (defaults to psm_table loaded above).
#' @param mzml_dir       Directory containing the mzML files.
#' @param top_n_label    Maximum number of matched peaks to label (ranked by
#'                       intensity); prevents overcrowding on dense spectra.
plot_annotated_psm <- function(full_sequence,
                               psm         = psm_table,
                               mzml_dir    = MZML_DIR,
                               top_n_label = 20) {
  rows <- psm[psm[["Full.Sequence"]] == full_sequence, ]
  if (nrow(rows) == 0) stop("No PSM found for: ", full_sequence)
  row  <- rows[which.max(rows$Score), ]

  # Experimental peaks
  pk     <- .load_peaks(row[["File.Name"]], row[["Scan.Number"]], mzml_dir)
  exp_df <- data.frame(mz = pk[, "mz"], intensity = pk[, "intensity"]) |>
    mutate(rel_int = intensity / max(intensity) * 100)

  # Matched ions: merge m/z and intensity tables on ion label
  matched_mz  <- .parse_ion_field(row[["Matched.Ion.Mass.To.Charge.Ratios"]])
  matched_int <- .parse_ion_field(row[["Matched.Ion.Intensities"]])
  matched <- left_join(matched_mz, matched_int, by = "ion",
                       suffix = c("_mz", "_int")) |>
    mutate(ion_type = case_when(
      startsWith(ion, "b") ~ "b",
      startsWith(ion, "y") ~ "y",
      TRUE                 ~ "other"
    ))

  # Annotate each experimental peak with the closest matched ion (≤20 ppm)
  exp_df$ion      <- NA_character_
  exp_df$ion_type <- NA_character_
  for (i in seq_len(nrow(matched))) {
    tol <- matched$value_mz[i] * 20 / 1e6
    j   <- which.min(abs(exp_df$mz - matched$value_mz[i]))
    if (abs(exp_df$mz[j] - matched$value_mz[i]) <= tol) {
      exp_df$ion[j]      <- matched$ion[i]
      exp_df$ion_type[j] <- matched$ion_type[i]
    }
  }

  label_df   <- exp_df |> filter(!is.na(ion)) |> slice_max(intensity, n = top_n_label)
  ion_colors <- c("b" = "#2166ac", "y" = "#d6604d", "other" = "#4dac26")

  p_spectrum <- ggplot(exp_df, aes(x = mz, xend = mz, y = 0, yend = rel_int)) +
    geom_segment(color = "grey75", linewidth = 0.35) +
    geom_segment(data = filter(exp_df, !is.na(ion)),
                 aes(color = ion_type), linewidth = 0.7) +
    geom_text(data = label_df,
              aes(y = rel_int + 1.5, label = ion, color = ion_type),
              size = 2.4, angle = 90, hjust = 0, vjust = 0.5,
              show.legend = FALSE) +
    scale_color_manual(values = ion_colors, name = "Ion series") +
    scale_y_continuous(expand = expansion(add = c(0, 3)), limits = c(0, 100)) +
    labs(
      title    = row[["Full.Sequence"]],
      subtitle = sprintf("File: %s  |  Scan: %d  |  Score: %.2f  |  Q-value: %.5f",
                         row[["File.Name"]], row[["Scan.Number"]],
                         row$Score, row$PEP_Qvalue),
      x = "m/z",
      y = "Relative Intensity (%)"
    ) +
    theme_bw() +
    theme(
      plot.title    = element_text(size = 9,  family = "mono"),
      plot.subtitle = element_text(size = 7.5, color = "grey40"),
      axis.title    = element_text(size = 10),
      axis.text     = element_text(size = 9),
      legend.position = "right"
    )

  p_seq <- .plot_sequence_ions(row[["Full.Sequence"]], matched)

  plot_grid(p_spectrum, p_seq, ncol = 1, rel_heights = c(5, 1.4))
}


# ── Example usage ─────────────────────────────────────────────────────────────

plot_annotated_psm("GKS[UniProt:Phosphoserine on S]PVPKS[UniProt:Phosphoserine on S]PVEEK")
