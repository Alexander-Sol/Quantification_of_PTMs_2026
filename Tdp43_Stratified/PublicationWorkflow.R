library(tidyverse)
library(limpa)
library(cowplot)
library(scales)
library(ggsignif)
library(ggrepel)

source("Supplemental/CustomScripts.R")

# Data loading, cleaning ----
# path <- r"(D:\Kelly_ALS_motor_nueron_dataset\MM1p1p0_Stratified_PubModsOnly\FlashLFQ_Norm_MbrFdr0p05\QuantifiedPeptides.tsv)" # Unmodified
#path <- r"(D:\Kelly_ALS_motor_nueron_dataset\MM1p1p4_GPTMD_Search_Carboxymethyl_Carbamido_2\FlashLFQ_Mbr0p05\QuantifiedPeptides.tsv)" 

path <- r"(Tdp43_Stratified/QuantData/DiversePtms/QuantifiedPeptides.tsv)" # Modified and unmodified

df <- read.csv(path, sep = '\t', row.names = NULL)
df <- df[grepl("Homo sapiens", df$Organism),] # Remove contaminant proteins

# Function to remove non-human attributions from ambiguous peptides
df <- CleanAmbiguousPeptides(df)

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

#Dropout curve modeling
dpcfit <- dpc(y.peptide)

# Perform protein quant using the dpc curve
y.protein <- dpcQuant(y.peptide, "ProteinGroup", dpc = dpcfit, verbose = TRUE)
y.protein <- AddProteinGeneNames(y.protein, y.peptide)
one_hit_wonders <- y.protein$genes$ProteinGroup[y.protein$genes$NPeptides == 1]

# PTM Analysis using diffSplice ------------------------------------------------
z <- dpcQuantByRow(y.peptide, "PeptideSequence", dpc = dpcfit, verbose = TRUE)

# Create the model design
design <- model.matrix(formula(~0 + Group + PMI + Sex), data = targets)
colnames(design) <- gsub("Group", "", colnames(design))
colnames(design)[ncol(design)] <- "Sex"

fit <- dpcDE(z, design, plot = FALSE)
contrast.matrix <- makeContrasts(ALS_vs_CTRL = "ALS - CTRL", levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
eb_fit <- eBayes(fit2)

du <- diffSplice(eb_fit, geneid="Gene", exonid="Sequence")
table_splice <- topSplice(du, coef = "ALS_vs_CTRL", number = Inf, test = "t")
table_splice_proteins <- topSplice(du, coef = "ALS_vs_CTRL", number = Inf)

# Exploring table_splice and optimizing filtering ------------------------------
table_splice$AllMod <- GetModsTextWReplacement(table_splice$PeptideSequence, all_mods = T) %>% sapply(., paste)
table_splice$Category <- ClassifyPeptideModsSpecific(table_splice$AllMod)

# Calculate the filtering stat
table_splice$TestStat <- (-1*log10(table_splice$FDR) + log10(0.05)) * abs(table_splice$logFC)

# annotate the position of the peptides within the protein
table_splice <- AddPositionColumns(table_splice)

# Levels of interest for filtering: 
b = 27 # Top 100
b = 12.725 # Top 250
b = 7.2 # Top 500
b = 3.78 # Top 1000
b = 1.1825 # Top 2500

# Take the top 1000 peptides filtered set
de_splice <- table_splice %>%
  filter(PropObs >= 0.09, TestStat > 3.78) %>%
  arrange(TestStat)

# Opposing Peptide analysis (Possible unneccesary) ----

# For unmodified peptides in de_splice, see if they can be explained by
# an opposing peptides
de_splice_unmod <- de_splice %>%
  filter(Category == "Unmodified")

de_splice_unmod$Opposing <- CheckForOpposingPeptidesTwoTable(de_splice_unmod, table_splice)

# For a random selection of unmodified proteins, how often is this true?
random_opposing <- CheckForOpposingPeptidesTwoTable(table_splice[table_splice$Category == "Unmodified", ] %>% 
                                   .[sample(nrow(.), 100), ], table_splice)

CheckForOpposingPeptidesTwoTable(table_splice[table_splice$Category == "Unmodified", ] %>% 
                                   .[sample(nrow(.), 100), ], table_splice) %>%
  table()


# ── Publication Plots: Peptide-Level Analysis ─────────────────────────────----
text_size <- 1.1

# Stringency levels: PropObs > 0.09 always applied first, then TestStat > b
# B values correspond to approximate top-N thresholds (see lines 89-95)
stringency_defs <- list(
  list(label = "All Peptides\n(n=30,059)",  b = -Inf),
  list(label = "Top 2,500\nPeptides",      b = 1.1825 ),
  list(label = "Top 1,000\nPeptides",      b = 3.78),
  list(label = "Top 250\nPeptides",       b = 12.725) #, list(label = "Top\n100",       b = 27)
)

base_filtered <- table_splice %>% filter(PropObs > 0.09)

apply_stringency <- function(b) {
  if (is.infinite(b) && b < 0) base_filtered else base_filtered %>% filter(TestStat > b)
}

category_levels_ordered <- c("Enzymatic Mod", "Carboxymethylation", "Non-Enzymatic Mod",
                             "Carbamidomethylation", "Unmodified")

# Build category frequency dataframe across all stringency levels ----
bar_df <- map_dfr(stringency_defs, function(lvl) {
  filt     <- apply_stringency(lvl$b)
  n_total  <- nrow(filt)
  #full_label <- paste0(lvl$label, "\n(n=", n_total, ")")
  full_label <- lvl$label # BE CAREFUL - the n for all peptides is hardcoded above
  filt %>%
    group_by(Category) %>%
    summarise(Count = n(), .groups = "drop") %>%
    complete(Category = factor(category_levels_ordered, levels = category_levels_ordered),
             fill = list(Count = 0)) %>%
    mutate(
      percent = 100 * Count / n_total,
      Group   = full_label,
      b       = lvl$b
    )
})

bar_df$Category <- factor(bar_df$Category, levels = category_levels_ordered)

# Preserve display order left-to-right (All → Top 100)
group_order <- map_chr(stringency_defs, function(lvl) {
  n <- nrow(apply_stringency(lvl$b))
  #paste0(lvl$label, "\n(n=", n, ")")
  lvl$label
})
bar_df$Group <- factor(bar_df$Group, levels = group_order)

small_threshold <- 5  # percent; segments below this skip in-bar labels

p_bar_stringency <- ggplot(bar_df, aes(x = Group, y = percent, fill = Category)) +
  geom_col(position = "stack", alpha = 0.8, color = "black", linewidth = 0.4) +
  geom_text(
    data = filter(bar_df, percent >= small_threshold),
    aes(label = sprintf("%d%%", round(percent))),
    position = position_stack(vjust = 0.5),
    size = 4, color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Enzymatic Mod"            = "#d73027",
      "Carboxymethylation"   = "#f46d43",
      "Non-Enzymatic Mod"        = "#fdae61",
      "Carbamidomethylation" = "#fee090",
      "Unmodified"           = "#4575b4"
    ),
    breaks = category_levels_ordered
  ) +
  coord_cartesian(clip = "off") +
  labs(y = "Percent of Peptides") +
  theme_minimal() +
  theme(
    axis.title.x  = element_blank(),
    axis.title.y  = element_text(size = 12 * text_size),
    axis.text.x   = element_text(size = 10 * text_size),
    axis.text.y   = element_text(size = 11 * text_size),
    legend.position = "right",
    legend.text   = element_text(size = 11 * text_size),
    legend.title  = element_blank(),
    plot.margin   = margin(5, 5, 5, 5)
  )

print(p_bar_stringency)
ggsave("ModBarBreakdown_Stringency.png", p_bar_stringency, width = 7, height = 5.5, dpi = 300)


# Volcano plot ----
b_top1000 <- 3.78  # TestStat threshold for top ~1000 peptides

volcano_df <- base_filtered %>%
  mutate(
    InTop1000 = TestStat > b_top1000,
    Highlight = case_when(
      InTop1000 & Gene %in% c("NEFH", "NEFM", "NEFL")                          ~ "Neurofilaments",
      InTop1000 & Gene == "GFAP"                                                 ~ "GFAP",
      InTop1000 & Gene %in% c("YWHAB", "YWHAG", "YWHAE", "YWHAH", "YWHAQ", "YWHAZ") ~ "14-3-3",
      TRUE                                                                        ~ NA_character_
    )
  )

highlight_df <- volcano_df %>% filter(!is.na(Highlight)) %>%
  mutate(Highlight = factor(Highlight, levels = c("GFAP", "14-3-3", "Neurofilaments")))

mod_colors <- c(
  "Unmodified"           = "#4575b4",
  "Carbamidomethylation" = "#fee090",
  "Non-Enzymatic Mod"            = "#fdae61",
  "Carboxymethylation"   = "#f46d43",
  "Enzymatic Mod"       = "#d73027"
)

# Hyperbolic boundary: TestStat = b_top1000
# Rearranged: -log10(FDR) = b / |logFC| - log10(0.05)
#                         = b / |logFC| + 1.301
x_seq <- c(seq(-3.4, -0.03, length.out = 500), seq(0.03, 3.4, length.out = 500))
boundary_df <- data.frame(logFC = x_seq) %>%
  mutate(y = b_top1000 / abs(logFC) - log10(0.05))

# Compress -log10(FDR) values above 20 into [20, 23]
y_linear_max <- 20
y_compressed_max <- 23

compress_y <- function(y_raw) {
  above <- y_raw > y_linear_max
  y_raw_max <- max(y_raw)
  y_out <- y_raw
  # Log-scale within compressed region so points spread evenly rather than
  # clustering near y_linear_max
  excess     <- y_raw[above] - y_linear_max
  excess_max <- y_raw_max    - y_linear_max
  y_out[above] <- y_linear_max +
    (y_compressed_max - y_linear_max) * log1p(excess) / log1p(excess_max)
  y_out
}

volcano_df <- volcano_df %>%
  mutate(y_display = compress_y(-log10(FDR)))

highlight_df <- highlight_df %>%
  mutate(y_display = compress_y(-log10(FDR)))

# Clamp boundary curve to the compressed y range so it doesn't expand the axis
boundary_df <- boundary_df %>% mutate(y = pmin(y, y_compressed_max+1))

# Axis breaks: linear region uses normal values; compressed region gets a label
# showing the actual max
y_raw_max_label <- round(max(-log10(volcano_df$FDR)))
y_breaks <- c(0, 5, 10, 15, 20, y_compressed_max)
y_labels <- c("0", "5", "10", "15", "20", as.character(y_raw_max_label))

p_volcano <- ggplot(volcano_df, aes(x = logFC, y = y_display)) +
  # Shaded band indicating compressed region
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = y_linear_max, ymax = y_compressed_max+0.25,
           fill = "grey92", alpha = 0.6) +
  # Background (below threshold): grey
  geom_point(data = filter(volcano_df, !InTop1000, is.na(Highlight)),
             color = "grey75", alpha = 0.4, size = 1.5, shape = 19) +
  # Top 1000, non-highlighted
  geom_point(data = filter(volcano_df, InTop1000, is.na(Highlight)),
             aes(color = Category), alpha = 0.8, size = 1.8, shape = 19) +
  # Highlighted shapes
  geom_point(data = filter(highlight_df, Highlight != "Neurofilaments"),
             aes(fill = Category, shape = Highlight),
             size = 2.3, alpha = 0.9, stroke = 0.6, color = "black") +
  geom_point(data = filter(highlight_df, Highlight == "Neurofilaments"),
             aes(fill = Category, shape = Highlight),
             size = 3.0, alpha = 0.9, stroke = 0.6, color = "black") +
  # Hyperbolic threshold boundary
  geom_line(data = filter(boundary_df, logFC < 0), aes(x = logFC, y = y),
            color = "grey30", linewidth = 0.5, linetype = "solid", inherit.aes = FALSE) +
  geom_line(data = filter(boundary_df, logFC > 0), aes(x = logFC, y = y),
            color = "grey30", linewidth = 0.5, linetype = "solid", inherit.aes = FALSE) +
  # Line marking start of compressed region
  geom_hline(yintercept = y_linear_max, linetype = "dotted", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = mod_colors, breaks = names(mod_colors)) +
  scale_fill_manual(values = mod_colors, guide = "none") +
  scale_shape_manual(values = c("GFAP" = 24, "14-3-3" = 22, "Neurofilaments" = 23),
                     name = NULL) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  labs(x = "Log2 Fold Change", y = "-Log10 FDR",
       title = "Differing Peptidoform Fractions") +
  theme_minimal() +
  theme(
    plot.title             = element_text(size = 14 * text_size),
    axis.title             = element_text(size = 12 * text_size),
    axis.text              = element_text(size = 11 * text_size),
    legend.position        = "inside",
    legend.position.inside = c(0.125, 0.75),
    legend.text            = element_text(size = 11 * text_size),
    legend.background      = element_rect(fill = alpha("white", 0.5), color = NA),
    legend.title           = element_blank(),
    legend.box             = "vertical",
    legend.margin          = margin(0, 0, 0, 0)
  ) +
  guides(color = guide_legend(order = 1, reverse = TRUE, override.aes = list(size = 3.5)),
         shape = guide_legend(order = 2, override.aes = list(size = 3.5))) +
  ylim(0, 23.5)

print(p_volcano)
ggsave("ModVolcano.png", p_volcano, width = 8.5, height = 6, dpi = 300)


# Mod bar charts (applied to de_splice — strict filter set) ----

# Non-Enzymatic bar chart
other_mod_df <- table(
  GetModsTextWReplacement(de_splice$PeptideSequence[de_splice$Category == "Non-Enzymatic Mod"]) %>%
    unlist() %>% str_split(., ", ") %>% unlist()
) %>% as.data.frame()
other_mod_df <- other_mod_df[!grepl("Carboxymeth|Carbamidomethy", other_mod_df$Var1), ]
other_mod_df <- other_mod_df[order(other_mod_df$Freq, decreasing = TRUE), ]
other_mod_df$mod <- factor(other_mod_df$Var1, levels = other_mod_df$Var1)
levels(other_mod_df$mod)[levels(other_mod_df$mod) == "Carbamidomethyl"]  <- "Carbamidomethylation"
levels(other_mod_df$mod)[levels(other_mod_df$mod) == "Carbamyl"]         <- "Carbamylation"
levels(other_mod_df$mod)[levels(other_mod_df$mod) == "Ammonia"]          <- "Ammonia Loss"
levels(other_mod_df$mod)[levels(other_mod_df$mod) == "Water"]            <- "Water Loss"
levels(other_mod_df$mod)[levels(other_mod_df$mod) == "Ubiquitination"]   <- "Dialkylation*"

p_other_mod <- ggplot(other_mod_df, aes(x = reorder(mod, Freq), y = Freq)) +
  geom_bar(stat = "identity", color = "black", fill = "#fdae61", alpha = 0.9, linewidth = 0.5) +
  coord_flip() +
  labs(title = "Non-Enzymatic Modifications in\nDifferentially Modified Peptides",
       x = "", y = "Instances of Modification Observed") +
  theme_minimal() +
  theme(
    plot.title   = element_text(size = 14 * text_size, hjust = 0,
                                margin = margin(b = 10, l = 16, unit = "pt")),
    axis.title.x = element_text(size = 12 * text_size),
    axis.text.x  = element_text(size = 11 * text_size),
    axis.text.y  = element_text(size = 12 * text_size)
  )

print(p_other_mod)
ggsave("OtherModBarChart.png", p_other_mod, width = 6, height = 4.5, dpi = 300)


# Enzymatic bar chart
uninteresting_mods <- c("Carbamidomethyl", "Carbamyl", "Carboxymethylation",
                        "Formylation", "Sodium", "Potassium", "Magnesium",
                        "Calcium", "Deamidation", "Deamidated", "Fe[III", "Iron",
                        "Ammonia", "Water", "Water Loss", "Carboxylation", "Oxidation",
                        "Carboxymethyllysine", "Nitrosylation")

bio_mod_df <- table(
  GetModsTextWReplacement(de_splice$PeptideSequence[de_splice$Category == "Enzymatic Mod"]) %>%
    unlist() %>% str_split(., ", ") %>% unlist() %>%
    .[!(. %in% uninteresting_mods)]
) %>% as.data.frame()
bio_mod_df <- bio_mod_df[order(bio_mod_df$Freq, decreasing = TRUE), ]
bio_mod_df$mod <- factor(bio_mod_df$Var1, levels = bio_mod_df$Var1)

p_bio_mod <- ggplot(bio_mod_df, aes(x = reorder(mod, Freq), y = Freq)) +
  geom_bar(stat = "identity", color = "black", fill = "#d73027", alpha = 0.9, linewidth = 0.5) +
  coord_flip() +
  labs(title = "Enzymatic Modifications in\nDifferentially Modified Peptides",
       x = "", y = "Instances of Modification Observed") +
  theme_minimal() +
  theme(
    plot.title   = element_text(size = 14 * text_size, hjust = 0,
                                margin = margin(b = 10, l = 16, unit = "pt")),
    axis.title.x = element_text(size = 12 * text_size),
    axis.text.x  = element_text(size = 11 * text_size),
    axis.text.y  = element_text(size = 12 * text_size)
  )

print(p_bio_mod)
ggsave("BioModBarChart.png", p_bio_mod, width = 6, height = 4.5, dpi = 300)


# ── Protein-level bar charts (top 10 proteins with bio-mod DA peptides) ─────────────────
gene_table_bio_mod <- de_splice$Gene[de_splice$Category == "Enzymatic Mod"] %>% 
  table() |> as.data.frame() %>% .[order(.$Freq, decreasing = T), ]
top10_genes <- gene_table_bio_mod$.[1:10]

bio_mods_of_interest <- c("Citrullination", "Acetylation", "Phosphorylation",
                          "Methylation", "Dimethylation", "Trimethylation",
                          "Nitrosylation", "Carboxylation", "Glutarylation",
                          "Succinylation", "Hydroxyproline", "Crotonylation",
                          "Butyrylation", "Hydroxybutyrylation")

bio_mod_order <- c("Acetylation", "Citrullination", "Phosphorylation",
                   "Methylation", "Trimethylation", "Hydroxyproline",
                   "Butyrylation", "Glutarylation")

# Enzymatic protein chart ----
bio_expanded <- de_splice %>%
  filter(Category == "Enzymatic Mod", Gene %in% top10_genes) %>%
  dplyr::select(Gene, AllMod) %>%
  mutate(Mod = strsplit(AllMod, ", ")) %>%
  unnest(Mod) %>%
  filter(Mod %in% bio_mods_of_interest) %>%
  mutate(Mod = factor(Mod, levels = rev(bio_mod_order)))

bio_gene_counts <- bio_expanded %>%
  group_by(Gene, Mod) %>%
  summarise(Count = n(), .groups = "drop")

gene_order_bio <- bio_gene_counts %>%
  group_by(Gene) %>%
  summarise(Total = sum(Count), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Gene)

bio_gene_counts$Gene <- factor(bio_gene_counts$Gene, levels = rev(gene_order_bio))

bio_protein_colors <- c(
  "Citrullination"      = "darkblue",
  "Acetylation"         = "forestgreen",
  "Methylation"         = "darkorchid2",
  "Dimethylation"       = "darkorchid3",
  "Trimethylation"      = "darkorchid4",
  "Phosphorylation"     = "cyan3",
  # "Nitrosylation"       = "darkorange",
  "Hydroxyproline"      = "deeppink3",
  # "Carboxylation"       = "saddlebrown",
  "Glutarylation"       = "sienna",
  # "Succinylation"       = "saddlebrown",
  # "Crotonylation"       = "saddlebrown",
  "Butyrylation"        = "sienna2"
)

p_bio_protein <- ggplot(arrange(bio_gene_counts, Mod), aes(x = Gene, y = Count, fill = Mod)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "gray20", linewidth = 0.3) +
  scale_fill_manual(values = bio_protein_colors, limits = bio_mod_order) +
  coord_flip() +
  labs(x = "", y = "Number of Mods on Differentially Modified Peptides", fill = "",
       title = "Proteins with Differentially Modified\nPeptides: Enzymatic Modifications") +
  theme_minimal() +
  scale_y_continuous(breaks = scales::breaks_width(2)) +
  theme(
    plot.title   = element_text(size = 14 * text_size),
    axis.title.x = element_text(size = 12 * text_size),
    axis.text.x  = element_text(size = 11 * text_size),
    axis.text.y  = element_text(size = 12 * text_size),
    legend.position = "right",
    legend.text  = element_text(size = 11 * text_size),
    legend.title = element_blank()
  ) +
  guides(fill = guide_legend(reverse = F))

print(p_bio_protein)
ggsave("BioModProteinChart.png", p_bio_protein, width = 6.5, height = 4.5, dpi = 300)


# Carboxymethylation protein chart ----
gene_table_cml <- de_splice$Gene[de_splice$Category == "Carboxymethylation"] %>% 
  table() |> as.data.frame() %>% .[order(.$Freq, decreasing = T), ]
top10_genes_cml <- gene_table_cml$.[1:10]

carbmeth_expanded <- de_splice %>%
  filter(Category == "Carboxymethylation", Gene %in% top10_genes_cml) %>%
  mutate(
    has_carbmeth_K = grepl("Carboxymethyllysine", AllMod),
    has_carbmeth_X = grepl("Carboxymethylation",  AllMod) & !grepl("Carboxymethyllysine", AllMod)
  ) %>%
  filter(has_carbmeth_K | has_carbmeth_X) %>%
  pivot_longer(cols = c(has_carbmeth_K, has_carbmeth_X),
               names_to = "ModType", values_to = "HasMod") %>%
  filter(HasMod) %>%
  mutate(ModType = case_when(
    ModType == "has_carbmeth_K" ~ "Carboxymethyllysine",
    ModType == "has_carbmeth_X" ~ "N-Terminal Carboxymethylation"
  ))

carbmeth_gene_counts <- carbmeth_expanded %>%
  group_by(Gene, ModType) %>%
  summarise(Count = n(), .groups = "drop")

gene_order_carbm <- carbmeth_gene_counts %>%
  group_by(Gene) %>%
  summarise(Total = sum(Count), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Gene)

carbmeth_gene_counts$Gene <- factor(carbmeth_gene_counts$Gene, levels = rev(gene_order_carbm))
carbmeth_gene_counts$ModType <- factor(carbmeth_gene_counts$ModType, 
                                       levels = rev(c("Carboxymethyllysine", "N-Terminal Carboxymethylation")))

p_carbmeth_protein <- ggplot(carbmeth_gene_counts, aes(x = Gene, y = Count, fill = ModType)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "gray20", linewidth = 0.3) +
  scale_fill_manual(values = c("Carboxymethyllysine"           = "#f46d43",
                               "N-Terminal Carboxymethylation" = "#fdae61")) +
  coord_flip() +
  labs(x = "", y = "Number of Mods on Differentially Modified Peptides", fill = "",
       title = "Proteins with Differentially Modified\nPeptides: Carboxymethylation") +
  theme_minimal() +
  scale_y_continuous(breaks = scales::breaks_width(4)) +
  theme(
    plot.title   = element_text(size = 14 * text_size),
    axis.title.x = element_text(size = 12 * text_size),
    axis.text.x  = element_text(size = 11 * text_size),
    axis.text.y  = element_text(size = 12 * text_size),
    legend.position = "bottom",
    legend.justification = c(0, 0),
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, -55),
    legend.text  = element_text(size = 11 * text_size)
  ) +
  guides(fill = guide_legend(reverse = T))

print(p_carbmeth_protein)
ggsave("CarboxymethylProteinChart.png", p_carbmeth_protein, width = 5, height = 4.5, dpi = 300)


# Non-Enzymatic protein chart ----
other_mod_expanded_protein <- de_splice %>%
  filter(Category == "Non-Enzymatic Mod", Gene %in% top12_genes) %>%
  mutate(mod_list = GetModsWReplacement(PeptideSequence)) %>%
  unnest(mod_list) %>%
  rename(Mod = mod_list) %>%
  filter(!grepl("Carboxymeth|Carbamidomethy", Mod)) %>%
  mutate(Mod = dplyr::recode(Mod,
                             "Carbamidomethyl" = "Carbamidomethylation",
                             "Carbamyl"        = "Carbamylation",
                             "Ammonia"         = "Ammonia Loss",
                             "Water"           = "Water Loss",
                             "Ubiquitination"  = "Ubiquitination*"
  ))

other_gene_counts <- other_mod_expanded_protein %>%
  group_by(Gene, Mod) %>%
  summarise(Count = n(), .groups = "drop")

gene_order_other <- other_gene_counts %>%
  group_by(Gene) %>%
  summarise(Total = sum(Count), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Gene)

other_gene_counts$Gene <- factor(other_gene_counts$Gene, levels = rev(gene_order_other))

unique_other_mods_protein <- unique(other_gene_counts$Mod)
other_protein_colors <- setNames(scales::hue_pal()(length(unique_other_mods_protein)),
                                 unique_other_mods_protein)

p_other_protein <- ggplot(other_gene_counts, aes(x = Gene, y = Count, fill = Mod)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "gray20", linewidth = 0.3) +
  scale_fill_manual(values = other_protein_colors) +
  coord_flip() +
  labs(x = "", y = "Number of DA Peptides with Non-Enzymatic Modifications", fill = "Modification",
       title = "Top 12 Proteins with Differentially Abundant\n'Non-Enzymatic Mod' Peptides") +
  theme_minimal() +
  theme(
    plot.title   = element_text(size = 14 * text_size),
    axis.title.x = element_text(size = 12 * text_size),
    axis.text.x  = element_text(size = 11 * text_size),
    axis.text.y  = element_text(size = 12 * text_size),
    legend.text  = element_text(size = 10 * text_size),
    legend.title = element_text(size = 11 * text_size)
  )

print(p_other_protein)
ggsave("OtherModProteinChart.png", p_other_protein, width = 6.5, height = 4.5, dpi = 300)

# ── Protein-Level DA: Diverse PTMs ────────────────────────────────────────----
# Use y.protein built from the diverse-PTM search above (already filtered to >= 2 peptides)
y.protein.filt <- y.protein[y.protein$genes$NPeptides >= 2, ]

fit_pro <- dpcDE(y.protein.filt, design = design, block = targets$Donor,
                 plot = FALSE, sample.weights = TRUE)
contrast.matrix.pro <- makeContrasts(ALS_vs_CTRL = "ALS - CTRL", levels = design)
fit_pro2  <- contrasts.fit(fit_pro, contrast.matrix.pro)
eb_fit_pro <- eBayes(fit_pro2)

table_protein_diverse <- topTable(eb_fit_pro, coef = "ALS_vs_CTRL", number = Inf)
logFC_cutoff <- 1.5
pval_cutoff  <- 0.05
table_protein_diverse$DEStatus <- ifelse(
  table_protein_diverse$adj.P.Val < pval_cutoff & table_protein_diverse$logFC >= logFC_cutoff,  "UP",
  ifelse(table_protein_diverse$adj.P.Val < pval_cutoff & table_protein_diverse$logFC <= -logFC_cutoff, "DOWN", "NotDE")
)

# ── Protein-Level DA: Limited PTMs ────────────────────────────────────────----
# Load the limited-PTM search (unmodified + Carbamidomethylation only)
path_limited <- r"(Tdp43_Stratified/QuantData/LimitedPtms/QuantifiedPeptides.tsv)"

df_lim <- read.csv(path_limited, sep = '\t', row.names = NULL)
df_lim <- df_lim[grepl("Homo sapiens", df_lim$Organism), ]
df_lim <- CleanAmbiguousPeptides(df_lim)

intensity_cols_lim <- grep("^Intensity_", colnames(df_lim), value = TRUE)
expr_lim  <- GetExpressionMatrix(df_lim, intensity_cols_lim)
genes_lim <- GetGeneTable(df_lim, expr_lim)

col_keep_lim <- colSums(!is.na(expr_lim)) > 2500
expr_lim <- expr_lim[, col_keep_lim]
targets_lim <- GetMetaDataTdpStratified(expr_lim)
expr_lim <- expr_lim[, colnames(expr_lim) %in% targets_lim$Sample]
keep_lim  <- rowSums(!is.na(expr_lim)) >= 3
expr_lim  <- expr_lim[keep_lim, ]
genes_lim <- genes_lim[keep_lim, ]
stopifnot(all(targets_lim$Sample == colnames(expr_lim)))

y.peptide.lim <- structure(
  list(E = expr_lim, 
       genes = genes_lim, 
       targets = targets_lim), class = "EList")
dpcfit_lim  <- dpc(y.peptide.lim)

y.protein.lim <- dpcQuant(y.peptide.lim, "ProteinGroup", dpc = dpcfit_lim, verbose = FALSE)
y.protein.lim <- AddProteinGeneNames(y.protein.lim, y.peptide.lim)
y.protein.lim <- y.protein.lim[y.protein.lim$genes$NPeptides >= 2, ]

design_lim <- model.matrix(formula(~0 + Group + PMI + Sex), data = targets_lim)
colnames(design_lim) <- gsub("Group", "", colnames(design_lim))
colnames(design_lim)[ncol(design_lim)] <- "Sex"

fit_lim  <- dpcDE(y.protein.lim, design = design_lim, block = targets_lim$Donor,
                  plot = FALSE, sample.weights = TRUE)
fit_lim2  <- contrasts.fit(fit_lim, makeContrasts(ALS_vs_CTRL = "ALS - CTRL", levels = design_lim))
eb_fit_lim <- eBayes(fit_lim2)

table_protein_limited <- topTable(eb_fit_lim, coef = "ALS_vs_CTRL", number = Inf)
table_protein_limited$DEStatus <- ifelse(
  table_protein_limited$adj.P.Val < pval_cutoff & table_protein_limited$logFC >= logFC_cutoff,  "UP",
  ifelse(table_protein_limited$adj.P.Val < pval_cutoff & table_protein_limited$logFC <= -logFC_cutoff, "DOWN", "NotDE")
)



# ── Publication Plots: Protein-Level ──────────────────────────────────────----
text_size <- 1.1

genes_of_interest_pro <- c("PADI2", "PRKCA", "FYN", "MBP", "GFAP", "CALM1",
                            "MAP2K4", "CAMK1D", "SIRT2")

# Volcano: diverse PTMs protein DA ----
highlight_pro <- table_protein_diverse %>% filter(Gene %in% genes_of_interest_pro)

p_pro_volcano <- ggplot(table_protein_diverse, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = DEStatus), alpha = 0.6, size = 2, shape = 19) +
  scale_color_manual(values = c("NotDE" = "grey", "UP" = "red", "DOWN" = "red")) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_point(data = highlight_pro, color = "darkred", size = 4, shape = 19) +
  geom_text_repel(data = highlight_pro, aes(label = Gene),
                  color = "black", size = 5, fontface = "bold") +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value",
       title = "Differentially Abundant Proteins") +
  theme_minimal() +
  theme(
    plot.title   = element_text(size = 14 * text_size),
    axis.title   = element_text(size = 12 * text_size),
    axis.text    = element_text(size = 11 * text_size),
    legend.position = "none"
  )

print(p_pro_volcano)
ggsave("ProteinVolcano_DiversePTMs.png", p_pro_volcano, width = 6, height = 5, dpi = 300)


# Comparison: diverse vs limited PTMs ----
diverse_de_genes <- table_protein_diverse$Gene[table_protein_diverse$DEStatus != "NotDE"]
limited_de_genes <- table_protein_limited$Gene[table_protein_limited$DEStatus != "NotDE"]

merged_pro <- merge(
  table_protein_diverse,
  table_protein_limited[, c("Gene", "logFC", "adj.P.Val")],
  by = "Gene", suffixes = c("", "_NoMod")
)
merged_pro$LogFC_Diff <- merged_pro$logFC - merged_pro$logFC_NoMod
merged_pro$ModDiff <- case_when(
  merged_pro$Gene %in% setdiff(limited_de_genes, diverse_de_genes) ~ "NoModDA",
  merged_pro$Gene %in% setdiff(diverse_de_genes, limited_de_genes) ~ "ModDA",
  TRUE ~ "Shared"
)
merged_pro <- merged_pro[order(merged_pro$ModDiff, decreasing = TRUE), ]

genes_of_interest_comp <- c("NEFM", "DRG1", "CLASP1", "EIF3J", "MARCKSL1", "EEF1A1", "USP9X")
highlight_comp <- merged_pro %>% filter(Gene %in% genes_of_interest_comp)

modiff_colors <- c("Shared" = "grey80", "ModDA" = "forestgreen", "NoModDA" = "darkred")

# Volcano colored by ModDiff ----
p_moddiff_volcano <- ggplot(merged_pro, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = ModDiff), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = modiff_colors) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = highlight_comp, aes(label = Gene),
                  color = "black", size = 5, fontface = "bold", box.padding = 0.5) +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value",
       title = "Changes in Differential Abundance when\nModifications are Considered") +
  theme_minimal() +
  theme(
    plot.title   = element_text(size = 14 * text_size),
    axis.title   = element_text(size = 12 * text_size),
    axis.text    = element_text(size = 11 * text_size),
    legend.position = "none"
  )

print(p_moddiff_volcano)
ggsave("ProteinVolcano_ModDiff.png", p_moddiff_volcano, width = 6, height = 5, dpi = 300)


# logFC scatter: diverse vs limited ----
p_logfc_scatter <- ggplot(merged_pro, aes(x = logFC, y = logFC_NoMod)) +
  geom_point(aes(color = ModDiff), alpha = 0.7, size = 2.5) +
  scale_color_manual(values = modiff_colors) +
  geom_hline(yintercept = 0, color = "grey20", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "grey20", linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_text_repel(data = highlight_comp, aes(label = Gene),
                  color = "black", size = 5, fontface = "bold",
                  box.padding = 1.5, point.padding = 0.5,
                  min.segment.length = 0, segment.color = "black") +
  labs(x = "Log2 Fold Change (Diverse PTMs)",
       y = "Log2 Fold Change (Limited PTMs)",
       title = "Comparison of Protein Log2FC With\nand Without Modified Peptides") +
  theme_minimal() +
  theme(
    plot.title   = element_text(size = 14 * text_size),
    axis.title   = element_text(size = 12 * text_size),
    axis.text    = element_text(size = 11 * text_size),
    legend.position = "none"
  )

print(p_logfc_scatter)
ggsave("ProteinLogFC_DiverseVsLimited.png", p_logfc_scatter, width = 6, height = 5, dpi = 300)


# ── DE Overlap: Diverse vs Limited PTMs ────────────────────────────────────----
n_diverse_only <- length(setdiff(diverse_de_genes, limited_de_genes))
n_limited_only <- length(setdiff(limited_de_genes, diverse_de_genes))
n_both         <- length(intersect(diverse_de_genes, limited_de_genes))

cat("DE protein counts:\n")
cat(sprintf("  Diverse PTMs only : %d\n", n_diverse_only))
cat(sprintf("  Limited PTMs only : %d\n", n_limited_only))
cat(sprintf("  Shared (both)     : %d\n", n_both))
cat(sprintf("  Total diverse     : %d\n", length(diverse_de_genes)))
cat(sprintf("  Total limited     : %d\n", length(limited_de_genes)))

# Area-proportional Euler diagram (circles scale with set size)
library(eulerr)
fit_euler <- euler(list("Diverse PTMs" = diverse_de_genes,
                        "Limited PTMs" = limited_de_genes))

p_venn <- plot(fit_euler,
               quantities = list(cex = 1.3),
               labels     = list(cex = 1.1),
               fills      = list(fill = c("forestgreen", "darkred"), alpha = 0.35),
               edges      = list(col = "grey20", lwd = 1.5),
               main       = "DE Proteins: Diverse vs. Limited PTMs")

print(p_venn)
png("ProteinDA_VennDiagram.png", width = 5, height = 4, units = "in", res = 300)
print(p_venn)
dev.off()


# ── Speculative Analysis of Specific Proteins ──────────────────────────────----
# Analysis of peptides specific to GFAP ----------------------------------------

# Quick function to set the mean "expression" of control samples to 0 for each peptide, 
# so that the heatmap will show relative changes from control
NormalizeToControl <- function(eList){
  control_samples <- eList$targets$Group == "CTRL"
  control_means <- rowMeans(eList$E[, control_samples], na.rm = TRUE)
  normalized_E <- sweep(eList$E, 1, control_means, FUN = "-")
  eList$E <- normalized_E
  return(eList)
}

CheckColumnName <- function(column_name, sample_names){
  str_split(column_name, "_")[[1]][[3]] %in% sample_names
}

RemoveOutliers <- function(eList, outlier_samples) {
  keep_samples <- ! sapply(eList$targets$Sample, function(s) CheckColumnName(s, outlier_samples))
  eList$E <- eList$E[, keep_samples, drop = FALSE]
  eList$targets <- eList$targets[keep_samples, , drop = FALSE]
  return(eList)
}

# Subset to GFAP peptides and create a new EList for just those peptides
GFAP_peptides <- de_splice$PeptideSequence[de_splice$Gene == "GFAP"]
GFAP_Data <- GetNormalizedEList(GFAP_peptides, z, y.protein, normalize_to_protein = T)

# We're going to remove four of the control samples that have high exp across the board
GFAP_Data <- RemoveOutliers(GFAP_Data, outlier_samples = c("4A3", "4A6", "4A4", "4A9"))

GFAP_Data <- NormalizeToControl(GFAP_Data)

# Start here to start fresh
GFAP_Ordered <- GFAP_Data[, order(GFAP_Data$targets$Donor)]

# Check for crazy outlier columns
colMeans(GFAP_Ordered$E) %>% .[order(.)]

# Rename the columns to something meaningful
donor_table = table(GFAP_Ordered$targets$Donor)
new_names = sapply(names(donor_table), function(d) {
  seq_num = seq_len(donor_table[d])
  paste0(d, "-", seq_num)
}) %>% unlist()
colnames(GFAP_Ordered) <- new_names

# Add peptide start and end positions to the metadata
GFAP_Ordered$genes <- AddPositionColumns(GFAP_Ordered$genes)

result <- ExportToJavaTreeView(
  elist       = GFAP_Ordered,
  out_prefix  = r"(C:\Users\Alex\Source\Repos\AlsMotorNeuronAnalysis\Tdp43_Stratified\PeptideHeatmapCols)",
  name_col    = "PeptideSequence",
  yorf_col    = "ProteinGroup",
  cluster_cols = TRUE
)

result$hc_arrays$labels[result$hc_arrays$order] # the order in which the samples were clustered
# Provides a gut-check to verify that control samples are clustering together and not being split across the dendrogram
PlotClusteredHeatmap(result, main = "Differentially Modified Peptides")


result <- ExportToJavaTreeView(
  elist       = GFAP_Ordered,
  out_prefix  = r"(C:\Users\Alex\Source\Repos\AlsMotorNeuronAnalysis\Tdp43_Stratified\PeptideHeatmapNoCols)",
  name_col    = "PeptideSequence",
  yorf_col    = "ProteinGroup",
  cluster_cols = FALSE#,
  #hclust_method = "average"
)

PlotClusteredHeatmap(result, main = "Differentially Modified Peptides")


# GFAP peptide check

# Find the unmodified DGEVIKESK row in y.peptide
idx <- grep("^DGEVIKESK$", y.peptide$genes$PeptideSequence)

# Check which samples observed it (non-NA)
observed <- y.peptide$E[idx, ] > 13.78033

# Break down by condition
table(y.peptide$targets$Group[observed])

table(y.peptide$targets$TDP_Lvl[observed])

# Also show total samples per condition for context
table(y.peptide$targets$Group)


# Doublecheck metadata # Metadata double check
supp_3 <- "4C11	4C9	4C8	4C7	4C6	4C5	4C4	4C3	4B11	4B10	4B9	4B8	4B7	4B5	4B4	4B3	4A11	4A10	4A9	4A8	4A7	4A6	4A5	4A4	4A3	1A8	3C5	3B10	3B9	3A10	3A8	3A6	3A5	3A4	2B5	1C9	1C4	1B10	1B6	1A9	1A7	1A5	3C10	3B11	3B7	3B6	3A11	2C10	2C9	2C8	2C7	2B11	2B10	2B9	2B8	2B3	2A11	2A9	2A8	2A7	2A5	2A4	2A3	1C10	1C5	1B4	1A11	1A10	2C3	1A6	3C9	3C8	3C6	3B8	3B5	3B4	3A7	2C6	2C5	2B7	2B6	1B8	1A3	3C7	3C4	3B3	3A9	3A3	2C4	2A10	2A6	1C8	1C7	1C6	1C3	1B11	1B9	1B7	1B5	1B3"
supp_names <- str_split(supp_3, "\t")[[1]]

int_cols <- colnames(df)[grep("Intensity_", colnames(df))]
file_names <- str_split(int_cols, "_") %>% sapply(., function(x) x[3])
# Check if all file names from the metadata are present in the data columns

intersect(supp_names, file_names)
setdiff(supp_names, file_names)
setdiff(file_names, supp_names)

# Log2 logit scale: For supplemental -------------------------------------------

# log2-logit function: log2(f / (1 - f))
log2_logit <- function(f) log2(f / (1 - f))

# Inverse: back-transform a coefficient to a proportion
inv_log2_logit <- function(x) 2^x / (1 + 2^x)

# Plot over a range of proportions
f_seq <- seq(0.001, 0.999, by = 0.001)

ggplot(data.frame(f = f_seq, y = log2_logit(f_seq)), aes(x = f, y = y)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
  labs(
    x = "Proportion of total protein signal (f)",
    y = "log2-logit(f)",
    title = "log2-logit transform"
  ) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(-10, 10, by = 1)) +
  theme_bw()


  
  # Show how a coefficient maps to a proportion change
  # given a baseline proportion in the reference condition
baseline_f <- c(0.01, 0.05, 0.10, 0.25, 0.50)
coef_seq <- seq(-4, 4, by = 0.1)

df <- expand.grid(baseline = baseline_f, coef = coef_seq) |>
  dplyr::mutate(
    baseline_logit = log2_logit(baseline),
    new_logit      = baseline_logit + coef,
    new_f          = inv_log2_logit(new_logit)
  )

ggplot(df, aes(x = coef, y = new_f, color = factor(baseline))) +
  geom_line() +
  geom_hline(yintercept = baseline_f, linetype = "dotted", color = "grey70") +
  labs(
    x = "diffSplice coefficient",
    y = "New proportion of protein signal",
    color = "Baseline\nproportion",
    title = "Proportion shift for a given diffSplice coefficient"
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent) +
  theme_bw()







