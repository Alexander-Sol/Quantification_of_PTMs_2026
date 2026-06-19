# Quantification of PTMs in Single-Cell Proteomics

Code companion to **_Improved discovery and quantification of post-translational modifications in
single-cell proteomics_** — an R-based re-analysis of single-cell proteomics data from ALS
(amyotrophic lateral sclerosis) motor neurons
([original dataset](https://doi.org/10.1016/j.celrep.2023.113636)).

📖 **Walkthrough / pipeline tutorial:** <https://alexander-sol.github.io/Quantification_of_PTMs_2026/>

The site walks through the entire workflow — from raw spectra to publication figures — so the
analysis can be reproduced or adapted to a new dataset.

## Pipeline

| Step | Tool | Purpose |
|------|------|---------|
| 1 | **MetaMorpheus** | Peptide/protein identification + PTM discovery (GPTMD) |
| 2 | **FlashLFQ + PIP-ECHO** | Label-free quantification with FDR-controlled match-between-runs |
| 3 | **limpa** | Dropout-curve modeling, protein summarization, probabilistic imputation |
| 4 | **R (limma / tidyverse)** | Differential abundance, PTM-occupancy analysis, figures |

Two parallel searches are compared throughout:

- **Diverse PTMs** — GPTMD + search against the UniProt human XML, allowing a broad PTM panel.
- **Limited PTMs** — GPTMD + search against the UniProt human FASTA, restricted to three variable
  modifications (deamidation, pyroglutamate, oxidation of methionine).

## Repository structure

```
├── MainManuscript/
│   ├── CustomScripts.R              # Data loading, cleaning, PTM parsing, metadata helpers
│   ├── analysis.R                   # Differential abundance + PTM analyses (run first)
│   ├── plotting.R                   # Publication figures (run after analysis.R)
│   ├── ptm_occupancy.R              # Alternative occupancy approach (protein-normalized)
│   ├── supplemental_ptm_diversity.R # Exports the full PTM-diversity supplemental table
│   ├── cache/                       # Cached intermediate objects (speed-ups)
│   └── extras/                      # Utilities NOT part of the core pipeline (see note below)
├── Supplemental/
│   ├── analysis.R                   # Supplementary analyses (needs MainManuscript objects)
│   └── plotting.R                   # Supplementary figures
├── Data/
│   ├── DiversePtms/ , LimitedPtms/  # FlashLFQ QuantifiedPeptides/Proteins (LFS)
│   ├── SequencePositionTable.tsv    # Peptide → sequence-position lookup
│   └── cache_*.rds                  # Cached EList / dpc objects (LFS)
├── site/                            # Quarto source for the walkthrough website
└── Quantification_of_PTMs.Rproj     # Open in RStudio to set the working directory
```

`MainManuscript/extras/` contains real but non-core utilities (`peptide_explorer.R`,
`spectrum_viewer.R`, `ProteinHclustExport.R`). They are **not** part of the reproducible pipeline:
they require either local raw spectra (e.g. `.mzML` files) or a separate code repository, and are
kept for reference only.

## Dependencies

- R ≥ 4.3 with **limpa** (Bioconductor), **limma**, and **tidyverse**
- Plotting also uses **cowplot**, **scales**, **ggrepel**, and **eulerr**
- `supplemental_ptm_diversity.R` additionally uses **writexl**

## Reproducing the analysis

1. **Clone with [Git LFS](https://git-lfs.com)** installed (the large `.tsv`/`.rds` files are tracked
   via LFS):
   ```bash
   git lfs install
   git clone https://github.com/Alexander-Sol/Quantification_of_PTMs_2026.git
   ```
2. **Open `Quantification_of_PTMs.Rproj` in RStudio.** This sets the working directory to the repo
   root, which every script assumes. (If running outside RStudio, `setwd()` to the repo root first.)
3. **Run the analysis:**
   ```r
   source("MainManuscript/analysis.R")   # builds analysis objects (uses Data/ cache if present)
   source("MainManuscript/plotting.R")   # writes figures to MainManuscript/Figures/
   ```

Steps 1–3 of the pipeline (MetaMorpheus → FlashLFQ → limpa) operate on raw spectra that are **not**
redistributed here. The FlashLFQ quantification tables that the R analysis consumes are included
under `Data/`. To regenerate them from raw data, follow the
[walkthrough](https://alexander-sol.github.io/Quantification_of_PTMs_2026/).

> **Note:** Figure 1 and the modification breakdowns read the raw MetaMorpheus `AllPSMs.psmtsv` /
> `AllPeptides.psmtsv` files (git-ignored). Copy your local search results into
> `Data/DiversePtms/` and `Data/LimitedPtms/` to regenerate those panels.

## The walkthrough website

The site under `site/` is a [Quarto](https://quarto.org) project. A GitHub Actions workflow
(`.github/workflows/publish.yml`) renders it and publishes to the `gh-pages` branch on every push.
To enable hosting once, set **Settings → Pages → Source** to the `gh-pages` branch. To preview
locally:

```bash
quarto preview site
```
