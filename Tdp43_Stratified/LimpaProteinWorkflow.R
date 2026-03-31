# PTM Occupancy for TDP stratified data - PTMOccupancy_TDP43_Stratified_3
library(tidyverse)
library(limpa)
library(cowplot)
library(scales)
library(ggsignif)
library(ggrepel)
library(dplyr)

source("../Supplemental/CustomScripts.R")
2
# Data loading, cleaning ----
# Load the data
# path <- r"(D:\Kelly_ALS_motor_nueron_dataset\MM1p1p0_Stratified_PubModsOnly\FlashLFQ_Norm_MbrFdr0p05\QuantifiedPeptides.tsv)" # Unmodified
path <- r"(D:\Kelly_ALS_motor_nueron_dataset\MM1p1p4_GPTMD_Search_Carboxymethyl_Carbamido_2\FlashLFQ_Mbr0p05\QuantifiedPeptides.tsv)" 
df <- read.csv(path, sep = '\t', row.names = NULL)
df <- df[grepl("Homo sapiens", df$Organism),] # Remove contaminant proteins

# Function to remove non-human attributions from ambiguous peptides
df <- CleanAmbiguousPeptides(df)
quant_pep_df <- df

# Identify intensity columns
intensity_cols <- grep("^Intensity_", colnames(df), value = TRUE)

# Get expression matrix  and annotation table
expr <- GetExpressionMatrix(df, intensity_cols)
genes <- GetGeneTable(df)

# Filter out columns with fewer than 2500 observed peptides ----
col_keep <- colSums(!is.na(expr)) > 2500
expr <- expr[,col_keep]

# Get  metadata table ----
targets <- GetMetaDataTdpStratified(expr)

# Remove samples without corresponding metadata
expr <- expr[,colnames(expr) %in% targets$Sample]

# Filter out rows that contain ONLY missing values
keep <- rowSums(!is.na(expr)) >= 3
expr <- expr[keep,]
genes <- genes[keep,]

# Ensure sample order matches expression matrix
stopifnot(all(targets$Sample == colnames(expr)))

# Create Elist objects ----

# Create peptide EList
y.peptide <- structure(
  list(E = expr, genes = genes, targets = targets),
  class = "EList"
)

# Dropout curve modeling
dpcfit <- dpc(y.peptide)

# Protein summarization
y.protein <- dpcQuant(y.peptide, "ProteinGroup", dpc = dpcfit, verbose = TRUE)

# Create the design for differential abundance testing ----
design <- model.matrix(formula(~0 + Group + PMI + Sex), data = targets)
colnames(design) <- gsub("Group", "", colnames(design))
colnames(design)[ncol(design)] <- "Sex"

# Protein level differential abundance ----
y.protein.no.OHW <- y.protein[y.protein$genes$NPeptides >= 2, ]

fit <- dpcDE(y.protein.no.OHW, design = design, block = targets$Donor, plot = TRUE, sample.weights = T)
contrast.matrix <- makeContrasts(ALS_vs_CTRL = "ALS - CTRL", levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
eb_fit <- eBayes(fit2)

# Extract results
table_protein <- topTable(eb_fit, coef = "ALS_vs_CTRL", number = Inf)
logFC_cutoff <- 1.5
pval_cutoff <- 0.05
table_protein$DEStatus <- ifelse(
  table_protein$adj.P.Val < pval_cutoff & table_protein$logFC >= logFC_cutoff, "UP",
  ifelse(table_protein$adj.P.Val < pval_cutoff & table_protein$logFC <= -logFC_cutoff, "DOWN", "NotDE")
)
de_table_protein <- table_protein[table_protein$DEStatus %in% c("UP", "DOWN"), ]
table(table_protein$DEStatus)

# Merge, add additional info
unique_df <- quant_pep_df[!duplicated(quant_pep_df$`Protein.Groups`), ]
merged_table_protein <- merge(unique_df[c("Gene.Names","Protein.Groups")], table_protein, by.x = "Protein.Groups",by.y = "ProteinGroup", all = FALSE)

# Clean the gene names
merged_table_protein$Gene.Names <- CleanGeneNames(merged_table_protein$Gene.Names)

# Write the output
write_tsv(merged_table_protein, file = r"(C:\Users\Alex\Source\Repos\AlsMotorNeuronAnalysis\FullWorkflow_ModifiedProteinResults.tsv)")

# Create volcano plot  ----
text_size = 1.1

p <- ggplot(merged_table_protein, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = DEStatus), alpha = 0.6, size = 2, shape = 19) +
  scale_color_manual(values = c("NotDE" = "grey", "UP" = "red", "DOWN" = "red")) +
  theme_minimal() +
  theme(legend.position = "none") +  # 👈 This removes the legend
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  ggtitle("Differentially Abundant Proteins") +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(plot.title   = element_text(size = 16 * text_size),   # title
        axis.title.x = element_text(size = 14 * text_size),                  # x-axis label
        axis.title.y = element_text(size = 14 * text_size),                  # y-axis label
        axis.text.x  = element_text(size = 14 * text_size),                  # x tick labels
        axis.text.y  = element_text(size = 12 * text_size),                   # y tick labels
        legend.position = "none")

genes_of_interest <- c("PADI2", "PRKCA", "FYN", "MBP", "GFAP", "CALM1", "MAP2K4", "CAMK1D", "SIRT2")

highlight <- subset(merged_table_protein, Gene.Names %in% genes_of_interest)
# Add highlighted points
p <- p +
  geom_point(data = highlight, aes(x = logFC, y = -log10(adj.P.Val)),
             color = "darkred", size = 4, shape = 19) +
  geom_text_repel(data = highlight,
                  aes(x = logFC, y = -log10(adj.P.Val), label = Gene.Names),
                  color = "black", size = 5, fontface = "bold")
p


# Re-run the protein level analysis  ----------------

# mod_pros <- read.csv(r"(C:\Users\Alex\Source\Repos\AlsMotorNeuronAnalysis\Tdp43_Stratified\Tdp43_ProteinDA_PMI_Sex_noMods.tsv)",
#                      sep = '\t', row.names = NULL)
# mod_pros$Gene.Names <- CleanGeneNames(mod_pros$Gene.Names)
# mod_pro_genes <- mod_pros$Gene.Names[mod_pros$DEStatus != "NotDE"]

# Read in the results from non-modified DAA
nomod_pros <- read.csv(r"(C:\Users\Alex\Source\Repos\AlsMotorNeuronAnalysis\Tdp43_Stratified\Tdp43_ProteinDA_PMI_Sex_noMods.tsv)",
                     sep = '\t', row.names = NULL)
nomod_genes <- nomod_pros$Gene.Names[nomod_pros$DEStatus != "NotDE"]
mod_genes <- merged_table_protein$Gene.Names[merged_table_protein$DEStatus != "NotDE"]

# Find the differences between the analyses
no_mod_da <- setdiff(nomod_genes, mod_genes)
mod_da <- setdiff(mod_genes, nomod_genes)

# Join the merged_table_protein with the nomod_pros based on Gene.Names, then find the difference in logFC between the two groups
merged_table_protein <- merge(merged_table_protein,
                               nomod_pros[, c("Gene.Names", "logFC", "adj.P.Val")],
                               by = "Gene.Names",
                               suffixes = c("", "_NoMod"))
merged_table_protein$LogFC_Diff <- merged_table_protein$logFC - merged_table_protein$logFC_NoMod

merged_table_protein$ModDiff <- ifelse(merged_table_protein$Gene.Names %in% no_mod_da, "NoModDA", ifelse(
  merged_table_protein$Gene.Names %in% mod_da, "ModDA",  "Shared"))
merged_table_protein <- merged_table_protein[order(merged_table_protein$ModDiff, decreasing = T), ]

# Find proteins that are DA in both but have opposite directionality
merged_table_protein$Directionality <- ifelse(
  ((merged_table_protein$logFC > 0 & merged_table_protein$logFC_NoMod < 0) |
    (merged_table_protein$logFC < 0 & merged_table_protein$logFC_NoMod > 0)) & 
    merged_table_protein$adj.P.Val < 0.05  &
    merged_table_protein$adj.P.Val_NoMod < 0.05, 
  "Opposite", "Same"
)

# There's only one Titin, and we only observe two peptides. Not gonna make a big 
# deal about that one

merged_table_protein <- merged_table_protein[order(merged_table_protein$ModDiff, decreasing=T),]
p <- ggplot(merged_table_protein, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = ModDiff), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("Shared" = "grey80", "ModDA" = "forestgreen", "NoModDA" = "darkred")) +
  theme_minimal() +
  theme(legend.position = "none") +  # 👈 This removes the legend
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  ggtitle("Changes in Differential Abundance when\nModifications are Considered") +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(plot.title   = element_text(size = 16 * text_size),   # title
        axis.title.x = element_text(size = 14 * text_size),                  # x-axis label
        axis.title.y = element_text(size = 14 * text_size),                  # y-axis label
        axis.text.x  = element_text(size = 14 * text_size),                  # x tick labels
        axis.text.y  = element_text(size = 12 * text_size),                   # y tick labels
        legend.position = "none")

genes_of_interest <- c("NEFM",  "DRG1", "CLASP1", "EIF3J", "MARCKSL1", "EEF1A1", "USP9X")
highlight <- subset(merged_table_protein, Gene.Names %in% genes_of_interest)
# Add highlighted points
p <- p +
  geom_text_repel(data = highlight,
                  aes(x = logFC, y = -log10(adj.P.Val), label = Gene.Names),
                  color = "black", size = 5, fontface = "bold", box.padding = 0.5)
p


# Create a plot showing the comparison of logFC with and without modified peptides ----
p <- ggplot(merged_table_protein, aes(x = logFC, y = logFC_NoMod)) +
  geom_point(aes(color = ModDiff), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c("Shared" = "grey80", "ModDA" = "forestgreen", "NoModDA" = "darkred")) +
  geom_hline(yintercept = 0, color = "grey20", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "grey20", linewidth = 0.5) +
  ylab("Log2 Fold Change (No Mods)") +
  xlab("Log2 Fold Change (With Mods)") +
  ggtitle("Comparison of Protein Log2FC With\nand Without Modified Peptides") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(plot.title   = element_text(size = 16 * text_size),   # title
        axis.title.x = element_text(size = 14 * text_size),                  # x-axis label
        axis.title.y = element_text(size = 14 * text_size),                  # y-axis label
        axis.text.x  = element_text(size = 14 * text_size),                  # x tick labels
        axis.text.y  = element_text(size = 12 * text_size),                   # y tick labels
        legend.position = "none")

genes_of_interest <- c("NEFM",  "DRG1", "CLASP1", "EIF3J", "MARCKSL1", "EEF1A1", "USP9X")
highlight <- subset(merged_table_protein, Gene.Names %in% genes_of_interest)
# Add highlighted points
p <- p +
  geom_text_repel(data = highlight,
                  aes(x = logFC, y = logFC_NoMod, label = Gene.Names),
                  color = "black", size = 5, fontface = "bold",
                  box.padding = 1.5,
                  point.padding = 0.5,
                  min.segment.length = 0,
                  segment.color = "black")
p



