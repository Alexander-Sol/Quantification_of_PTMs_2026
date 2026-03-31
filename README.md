# Quantification_of_PTMs_2026

Code companion for: **Improved discovery and quantification of post-translational modifications in single-cell proteomics**

R-based re-analysis of single-cell proteomics data from ALS (Amyotrophic Lateral Sclerosis) motor neurons, original publication [here](https://doi.org/10.1016/j.celrep.2023.113636)

## Dependencies

- **limpa**: Dropout curve modeling and differential expression (`dpc()`, `dpcQuant()`, `dpcDE()`)
- **limma**: Linear models for differential abundance
- **tidyverse / ggplot2**: Data manipulation and visualization

## Repository Structure

```
├── Supplemental/
│   └── CustomScripts.R           # Data loading, cleaning, PTM parsing, and metadata functions
├── Tdp43_Stratified/
│   ├── PublicationWorkflow.R     # Primary analysis: PTM analysis, protein DA, figures 2–4
│   ├── LimmaProteinWorkflow.R    # Limma-based protein-level differential abundance
│   ├── LimpaProteinWorkflow.R    # Limpa-based protein-level differential abundance
│   ├── SequencePositionTable.tsv # Reference table mapping peptides to sequence positions
│   └── QuantData/
│       ├── DiversePtms/
│       │   ├── QuantifiedPeptides.tsv   # FlashLFQ+PIP-ECHO output (diverse PTM search); LFS
│       │   └── QuantifiedProteins.tsv  # Protein-level quantification
│       └── LimitedPtms/
│           ├── QuantifiedPeptides.tsv   # FlashLFQ+PIP-ECHO output (limited PTM search); LFS
│           └── QuantifiedProteins.tsv  # Protein-level quantification
└── Stratified_PSMs_mod_vs_nomod.R    # Figure 1: PSM and peptide counts across search strategies
```

## Data

Input data (raw `.psmtsv` /  files from MetaMorpheus) resides on the Smith lab network, will be added to git LFS shortly.
`\\bison.chem.wisc.edu\share\Projects\Kelly_ALS_motor_nueron_dataset\`

The `QuantData/` files in this repo are the pre-processed FlashLFQ outputs used as direct inputs to the R analysis scripts. The large `QuantifiedPeptides.tsv` files are stored via Git LFS.

## Usage

Open the project in RStudio. Scripts source helper functions relatively:

```r
source("../Supplemental/CustomScripts.R")
```

Run `PublicationWorkflow.R` for the primary PTM and protein-level analysis. Run `Stratified_PSMs_mod_vs_nomod.R` for Figure 1 (update the `psmtsv` paths to point to your local copies of the MetaMorpheus search results).

## Dataset

~100 single motor neurons stratified by TDP-43 pathology severity (CTL / NON / MLD / MOD / SEV), from 5 post-mortem donors. Two parallel searches:
- **DiversePtms**: GPTMD+Search using the UniProt human XML allowing a broad set of PTMs
- **LimitedPtms**: GPTMD+Search  using the UniProt human fasta restricted to 3 variable PTMs (Deamidation, pyroglutamate, oxidation of methionine) 
