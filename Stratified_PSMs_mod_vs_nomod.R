# PSM Counts: Modified vs Unmodified ----
library(tidyverse)

# Load the data
modified_psm_path <- r"(D:\Kelly_ALS_motor_nueron_dataset\MM1p1p4_GPTMD_Search_Carboxymethyl_Carbamido_2\Task2-SearchTask\AllPSMs.psmtsv)"
notmodified_psm_path <- r"(D:\Kelly_ALS_motor_nueron_dataset\MM1p1p4_GPTMDSearch_PubMods_HumanFasta\Task2-SearchTask\AllPSMs.psmtsv)"

output_path <- r"(C:\Users\Alex\Source\Repos\AlsMotorNeuronAnalysis\PSM_Plots)"
setwd(output_path)

# Read both PSM files
modified_df <- read.csv(modified_psm_path, sep = '\t', row.names = NULL)
notmodified_df <- read.csv(notmodified_psm_path, sep = '\t', row.names = NULL)

# Filter by pep-q values
filter_mm_results <- function(mm_results, q_threshold = 1, pep_q_threshold = 0.01) {
  filtered_results <- mm_results %>%
    filter(PEP_QValue <= pep_q_threshold) %>%
    filter(QValue <= q_threshold) %>%
    filter(!grepl("C|D", Decoy.Contaminant.Target))
  return(filtered_results)
}

modified_df <- filter_mm_results(modified_df)
notmodified_df <- filter_mm_results(notmodified_df)

# Count PSMs per file for each dataset
modified_counts <- modified_df %>%
  group_by(`File.Name`) %>%
  summarise(count = n()) %>%
  mutate(Type = "Search 2: Diverse PTMs")

notmodified_counts <- notmodified_df %>%
  group_by(`File.Name`) %>%
  summarise(count = n()) %>%
  mutate(Type = "Search 1: Limited PTMs")

# Combine the counts
combined_counts <- bind_rows(modified_counts, notmodified_counts)

# Extract just the file name (without path) for cleaner labels
combined_counts <- combined_counts %>%
  mutate(File = basename(`File.Name`))

combined_counts$File <- combined_counts$File %>%
  str_split(., "_") %>% sapply(., function(x) x[2])

mod_counts <- combined_counts[combined_counts$Type == "Search 2: Diverse PTMs", ]
mod_counts <- mod_counts[order(mod_counts$count), ]
mod_counts$File <- factor(mod_counts$File, levels = mod_counts$File)

combined_counts$File <- factor(combined_counts$File, levels = levels(mod_counts$File))
combined_counts <- combined_counts[order(combined_counts$Type, decreasing = T), ]

# Create overlaid bar chart
ggplot(combined_counts, aes(x = File, y = count, fill = Type)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1.8, "lines")) +
  labs(
    title = "PSMs per Cell: Diverse PTMs vs. Limited PTMs",
    x = "Sample Number",
    y = "",
    fill = "PSM Type"
  ) +
  scale_fill_manual(values = c("Search 2: Diverse PTMs" = "#78CFDE", "Search 1: Limited PTMs" = "#44648E"),
                    guide = guide_legend(reverse = TRUE))
ggsave("PSM_Counts_per_Cell.png", width = 7, height = 5, dpi = 300)

# Total PSM counts: Modified vs Unmodified ----
total_psm_counts <- combined_counts %>%
  group_by(Type) %>%
  summarise(Total = sum(count))

ggplot(total_psm_counts, aes(x = 1, y = Total, fill = Type)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.7, width = 0.6) +
  theme_minimal() +
  labs(
    title = "Total\nPSMs",
    x = "",
    y = "Number of PSMs",
    fill = "PSM Type"
  ) +
  scale_fill_manual(values = c("Modified" = "#78CFDE", "Unmodified" = "#44648E")) +
  coord_cartesian(clip = "off") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.position = "none",
        plot.title.position = "plot") +
  geom_text(aes(label = scales::comma(Total)), vjust = -0.5, size = 4,
            position = position_identity())
ggsave("Total_PSM_Counts.png", width = 1.5, height = 5, dpi = 300)

# Repeat the process but for Peptides ----
# Load the data
modified_peptide_path <- r"(D:\Kelly_ALS_motor_nueron_dataset\MM1p1p4_GPTMD_Search_Carboxymethyl_Carbamido_2\Task2-SearchTask\AllPeptides.psmtsv)"
notmodified_peptide_path <- r"(D:\Kelly_ALS_motor_nueron_dataset\MM1p1p4_GPTMDSearch_PubMods_HumanFasta\Task2-SearchTask\AllPeptides.psmtsv)"

# Read both Peptides files
modified_pep <- read.csv(modified_peptide_path, sep = '\t', row.names = NULL)
notmodified_pep <- read.csv(notmodified_peptide_path, sep = '\t', row.names = NULL)

modified_pep <- filter_mm_results(modified_pep)
notmodified_pep <- filter_mm_results(notmodified_pep)


# Count unique peptides per file: filter PSMs to sequences validated in the
# peptides file, then count distinct full sequences per file
modified_counts <- modified_df %>%
  filter(`Full.Sequence` %in% modified_pep$`Full.Sequence`) %>%
  group_by(`File.Name`) %>%
  distinct(`Full.Sequence`, .keep_all = TRUE) %>%
  summarise(count = n()) %>%
  mutate(Type = "Search 2: Diverse PTMs")

notmodified_counts <- notmodified_df %>%
  filter(`Full.Sequence` %in% notmodified_pep$`Full.Sequence`) %>%
  group_by(`File.Name`) %>%
  distinct(`Full.Sequence`, .keep_all = TRUE) %>%
  summarise(count = n()) %>%
  mutate(Type = "Search 1: Limited PTMs")

# Combine the counts
combined_counts_pep <- bind_rows(modified_counts, notmodified_counts)

# Extract just the file name (without path) for cleaner labels
combined_counts_pep <- combined_counts_pep %>%
  mutate(File = basename(`File.Name`))

combined_counts_pep$File <- combined_counts_pep$File %>%
  str_split(., "_") %>% sapply(., function(x) x[2])

mod_counts <- combined_counts_pep[combined_counts_pep$Type == "Search 2: Diverse PTMs", ]
mod_counts <- mod_counts[order(mod_counts$count), ]
mod_counts$File <- factor(mod_counts$File, levels = mod_counts$File)

combined_counts_pep$File <- factor(combined_counts_pep$File, levels = levels(mod_counts$File))
combined_counts_pep <- combined_counts_pep[order(combined_counts_pep$Type, decreasing = T), ]

# Create overlaid bar chart
ggplot(combined_counts_pep, aes(x = File, y = count, fill = Type)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1.8, "lines")) +
  labs(
    title = "Peptides per Cell: Diverse PTMs vs. Limited PTMs",
    x = "Sample Number",
    y = "",
    fill = "Peptide Type"
  ) +
  scale_fill_manual(values = c("Search 2: Diverse PTMs" = "#78CFDE", "Search 1: Limited PTMs" = "#44648E"),
                    guide = guide_legend(reverse = TRUE))
ggsave("Peptide_Counts_per_Cell.png", width = 7, height = 5, dpi = 300)

# Total Peptide counts: Modified vs Unmodified ----
total_pep_counts <- combined_counts_pep %>%
  group_by(Type) %>%
  summarise(Total = sum(count))

ggplot(total_pep_counts, aes(x = 1, y = Total, fill = Type)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.7, width = 0.6) +
  theme_minimal() +
  labs(
    title = "Total\nPeptides",
    x = "",
    y = "Total Number of Peptides",
    fill = "Peptide Type"
  ) +
  scale_fill_manual(values = c("Modified" = "#78CFDE", "Unmodified" = "#44648E")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.position = "none",
        plot.title.position = "plot") +
  geom_text(aes(label = scales::comma(Total)), vjust = -0.5, size = 4,
            position = position_identity())
ggsave("Total_Peptide_Counts.png", width = 1.5, height = 5, dpi = 300)


# Alyklation analysis ----

# Count the number of off-target alkylation psms 
alkylated_pep <- modified_pep[grepl("Carbamidomethyl on [^C]", modified_pep$Full.Sequence), ]
alkylated_pep$CleanMods <- gsub("Carbamidomethyl on [A-Z]", "", alkylated_pep$Mods)

# Count rows where alkylated_pep$CleanMods is not empty or whitespace
alkylated_plus <- alkylated_pep[!grepl("^\\s*$", alkylated_pep$CleanMods), ]

# Remove all psms without a corresponding entry in the respective peptides file
modified_df <- modified_df[modified_df$`Full.Sequence` %in% modified_pep$`Full.Sequence`, ]
notmodified_df <- notmodified_df[notmodified_df$`Full.Sequence` %in% notmodified_pep$`Full.Sequence`, ]


# Modification Type Analysis ----
source(r"(C:\Users\Alex\Source\Repos\AlsMotorNeuronAnalysis\Supplemental\CustomScripts.R)")

# Filter for high-confidence PSMs (PEP_QValue <= 0.01)
modified_filtered <- modified_df %>%
  filter(PEP_QValue <= 0.01)

# Extract and count modifications
# Use Full.Sequence column which contains modification annotations like [Common:Acetylation on K]
mod_list <- GetModsTextWReplacement(modified_filtered$Full.Sequence, all_mods = T)

mod_counts <- mod_list %>%
  unlist() %>% 
  str_split(., ", ") %>%
  unlist() %>%
  table() %>% # Flatten the list and count occurrences
  as.data.frame() %>% # Count each modification type
  arrange(desc(Freq))

# Name column
names(mod_counts)[1] <- "all_mods"

# Remove un-mod row
mod_counts <- mod_counts[mod_counts$all_mods != "", ]

# Classify each modification as Biological or Other
mod_counts$Category <- ifelse(mod_counts$all_mods %in% biological_mods,
                               "Biological Mod", "Other Mod")

# Rename for readability
mod_counts$all_mods_renamed <- as.character(mod_counts$all_mods)
mod_counts$all_mods_renamed[mod_counts$all_mods_renamed == "Carbamidomethyl"] <- "Carbamido-\nmethylation"
mod_counts$all_mods_renamed[mod_counts$all_mods_renamed == "Carboxymethyllysine"] <- "Carboxy-\nmethyllysine"
mod_counts$all_mods_renamed[mod_counts$all_mods_renamed == "Carboxymethylation"] <- "Carboxy-\nmethylation"
mod_counts$all_mods_renamed[mod_counts$all_mods_renamed == "Carbamyl"] <- "Carbamylation"
mod_counts$all_mods_renamed[mod_counts$all_mods_renamed == "Ammonia"] <- "Ammonia\nLoss"
mod_counts$all_mods_renamed[mod_counts$all_mods_renamed == "Water"] <- "Water\nLoss"
mod_counts$all_mods_renamed[mod_counts$all_mods_renamed == "Iron"] <- "Iron\nAdduct"
mod_counts$all_mods_renamed[mod_counts$all_mods_renamed == "Sodium"] <- "Sodium\nAdduct"
mod_counts$all_mods_renamed[mod_counts$all_mods_renamed == "Calcium"] <- "Calcium\nAdduct"

# Category colors from PTMOccupancy
category_colors <- c("Biological Mod" = "#d73027", "Other Mod" = "#f46d43")

# Biological Modifications plot ----
bio_mods_detected <- mod_counts %>% filter(Category == "Biological Mod")

bio_mods_detected$all_mods_renamed <- factor(bio_mods_detected$all_mods_renamed,
  levels = bio_mods_detected$all_mods_renamed[order(bio_mods_detected$Freq)])

# Take 12 most common
bio_mods_detected <- bio_mods_detected[1:12,]

ggplot(bio_mods_detected, aes(x = all_mods_renamed, y = Freq)) +
  geom_bar(stat = "identity", fill = "#d73027", alpha = 0.9, color = "gray20", linewidth = 0.8) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Biological Modifications Detected",
    x = "",
    y = "Number of PTMs\nObserved Across PSMs"
  ) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 14, margin = margin(l = 5,  unit = "pt")), # Margin to make bio mod plot match "other" mod plot
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    plot.title.position = "plot"
  ) +
  scale_y_log10()
ggsave("Biological_Modifications_Detected.png", width = 5, height = 5, dpi = 300)

# Other Modifications plot ----
other_mods_detected <- mod_counts %>% filter(Category == "Other Mod")

other_mods_detected$all_mods_renamed <- factor(other_mods_detected$all_mods_renamed,
  levels = other_mods_detected$all_mods_renamed[order(other_mods_detected$Freq)])

# Take 10 most common
other_mods_detected <- other_mods_detected[1:10, ]

ggplot(other_mods_detected, aes(x = all_mods_renamed, y = Freq)) +
  geom_bar(stat = "identity", fill = "#f46d43", alpha = 0.9, color = "gray20", linewidth = 0.8) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Other Modifications Detected",
    x = "",
    y = "Number of PTMs\nObserved Across PSMs"
  ) +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    plot.title.position = "plot",
  ) +
  scale_y_log10()
ggsave("Other_Modifications_Detected.png", width = 5, height = 5, dpi = 300)



# Carboxymethylation analysis

# Load the data
modified_psm_path <- r"(D:\Kelly_ALS_motor_nueron_dataset\MM1p1p4_GPTMD_Search_Carboxymethyl_Carbamido_2\Task2-SearchTask\AllPSMs.psmtsv)"

# Read both PSM files
modified_df <- read.csv(modified_psm_path, sep = '\t', row.names = NULL)

# Filter by pep-q values
filter_mm_results <- function(mm_results, q_threshold = 1, pep_q_threshold = 0.01) {
  filtered_results <- mm_results %>%
    filter(PEP_QValue <= pep_q_threshold) %>%
    filter(QValue <= q_threshold) %>%
    filter(!grepl("C|D", Decoy.Contaminant.Target))
  return(filtered_results)
}

modified_df <- filter_mm_results(modified_df)

term_carboxymethyl <- modified_df[grepl("Carboxymethylation on X", modified_df$Full.Sequence), ]
