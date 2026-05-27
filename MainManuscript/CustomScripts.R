library(dplyr)
library(stringr)

# Functions for FlashLFQ + Limpa scripts

# Data Loading and Cleaning ----

# Clipboard <- function(data)
# {
#   write.table(data, "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
# }

Clipboard <- function(data)
{
  clip <- pipe("pbcopy", "w")                                                 
  write.table(data, clip, sep="\t", row.names=FALSE, col.names=FALSE)         
  close(clip)
}

# Function to remove non-human attributions from ambiguous peptides
# i.e., if a peptide is ambiguous between a human and contaminant protein, the 
# contaminant protein annotation is removed
RemoveAmbiguity <- function(row){
  protein_groups <- row["Protein.Groups"] %>% strsplit(., ";") %>% unlist()
  species <- row["Organism"] %>% strsplit(., ";") %>% unlist()
  human_entries <- sapply(species, function(x) grep(x, "Homo sapiens")) %>% unlist()
  if(length(human_entries) != 1) # if there are multiple human entries, then the peptide is ambiguous between multiple human proteins and should be discarded
    return(row)
  human_position <- match("Homo sapiens", species)
  # remove the non-human entries to disambiguate the peptide
  row["Protein.Groups"] <- protein_groups[human_position]
  row["Organism"] <- "Homo sapiens"
  return(row)
}

# Takes in the QuantifiedPeptides data frame, 
# And handles ambiguity, removing human/human ambigous peptides to prevent
# Use in downstream quantification
CleanAmbiguousPeptides <- function(df){
  # Here, we process peptides that are ambiguous between a human and contaminant protein,
  # and assign it to be soley human
  df_processed <- apply(df, MARGIN = 1, FUN = RemoveAmbiguity)
  df <- df_processed %>% t() %>% as.data.frame()
  df <- df[!grepl(";", df$Protein.Groups), ] # remove all peptides that are ambiguous between two or more human proteins
  return(df)
}

# Takes in a QuantifiedPeptides data frame and a character vector
# that identifies the columns where intensity is reported. Returnes
# a log2 transformed expression matrix
GetExpressionMatrix <- function(df, intensity_cols){
  expr <- apply(df[, intensity_cols], 2, as.numeric)
  rownames(expr) <- df$Sequence
  expr[expr == 0] <- NA
  expr <- log2(expr)
  return(expr)
}

# Strips away the "primary:" prefix from gene labels and returns the first gene listed
CleanGeneNames <- function(gene_names){
  gene_names <- gene_names %>% str_split(., ":") %>% sapply(., function(x) x[2])
  gene_names <- gene_names %>% str_split(., ";") %>% sapply(., function(x) x[1])
  gene_names <- gene_names %>% str_split(., ",") %>% sapply(., function(x) x[1])
  return(gene_names)
}

# Returns a table that contains metadata about the protein and gene corresponding 
# to each peptide
GetGeneTable <- function(df, expr_matrix = expr){
  genes <- df[, c("Sequence", "Base.Sequence", "Protein.Groups", "Gene.Names")]
  colnames(genes) <- c("PeptideSequence", "BaseSequence", "ProteinGroup", "Gene")
  rownames(genes) <- df$Sequence
  genes$Gene <- CleanGeneNames(genes$Gene)
  genes <- genes[rownames(expr_matrix), , drop = FALSE]
}

# Adds ProteinGroup (from E rownames) and Gene to y.protein$genes using
# y.peptide$genes as the accession -> gene name lookup table
AddProteinGeneNames <- function(y.protein, y.peptide) {
  accessions <- rownames(y.protein$E)
  if (!"ProteinGroup" %in% colnames(y.protein$genes)) {
    y.protein$genes$ProteinGroup <- accessions
  }
  lookup <- y.peptide$genes[!duplicated(y.peptide$genes$ProteinGroup),
                             c("ProteinGroup", "Gene")]
  idx <- match(y.protein$genes$ProteinGroup, lookup$ProteinGroup)
  y.protein$genes$Gene <- lookup$Gene[idx]
  return(y.protein)
}

GetMetaDataPilotStudy <- function(expr, keep_b6 = F){
  # Assign conditions - Load in metadata (copied and pasted from the supplemental #2)
  short_name_str <- "1A3	1A4	1A5	2A6	2A7	2A8	1B3	1B4	1B5	2B6	2B7	2B8	1A6	1A7	1A8	2A3	2A4	2A5	1B6	1B7	1B8	2B3	2B4	2B5	1C6	1C7	1C8	2C3	2C4	2C5	1C3	1C4	1C5	2C6	2C7	2C8"
  short_names <- str_split(short_name_str, "\t") %>% unlist()
  
  condition_str <- "ALS	ALS	ALS	ALS	ALS	ALS	ALS	ALS	ALS	ALS	ALS	ALS	ALS	ALS	ALS	ALS	ALS	ALS	CTRL	CTRL	CTRL	CTRL	CTRL	CTRL	CTRL	CTRL	CTRL	CTRL	CTRL	CTRL	CTRL	CTRL	CTRL	CTRL	CTRL	CTRL"
  conditions <- str_split(condition_str, "\t") %>% unlist()
  
  donor_id_str <- "ALS #1	ALS #1	ALS #1	ALS #1	ALS #1	ALS #1	ALS #2	ALS #2	ALS #2	ALS #2	ALS #2	ALS #2	ALS #3	ALS #3	ALS #3	ALS #3	ALS #3	ALS #3	CTL #1	CTL #1	CTL #1	CTL #1	CTL #1	CTL #1	CTL #2	CTL #2	CTL #2	CTL #2	CTL #2	CTL #2	CTL #3	CTL #3	CTL #3	CTL #3	CTL #3	CTL #3"
  donor_ids <- str_split(donor_id_str, "\t") %>% unlist()

  name_condition_table <- data.frame(
    ShortName = short_names,
    Group = conditions,
    Donor = donor_ids
  )
  
  # Get shortened names of columns to match up to the short names from the supplemental data
  column_names <- colnames(expr) %>% strsplit("_") %>% sapply(function(x) x[5]) %>% unlist()
  
  # remove metadata that doesn't correspond to any files
  name_condition_table <- name_condition_table[name_condition_table$ShortName %in% column_names,]
  
  # order the new metadata table
  name_condition_table <- name_condition_table[match(column_names, name_condition_table$ShortName), ]
  
  # Metadata (Group + Donor) 
  targets <- data.frame(
    Sample = colnames(expr),
    Group = factor(name_condition_table$Group),
    Donor = factor(name_condition_table$Donor)
  )
  
  # Removes pooled samples
  targets <- targets[!is.na(targets$Group), ]
  
  if(!keep_b6){
    #Removes one definite outlier
    targets <- targets[targets$Sample != "Intensity_Biogen_ALS_Pilot_1B6_020420", ]
  }
  

  pmi_info <- list(
    "ALS #1" = 9.3,
    "ALS #2" = 6.3,
    "ALS #3" = 12,
    "CTL #1" = 15, 
    "CTL #2" = 14,
    "CTL #3" = 12.91
  )
  
  # Attempt to control for PMI
  age_info <- list(
    "ALS #1" = 80,
    "ALS #2" = 59,
    "ALS #3" = 66,
    "CTL #1" = 56, 
    "CTL #2" = 58,
    "CTL #3" = 63
  )
  
  # Attempt to control for PMI
  sex_info <- list(
    "ALS #1" = "M",
    "ALS #2" = "M",
    "ALS #3" = "F",
    "CTL #1" = "M", 
    "CTL #2" = "M",
    "CTL #3" = "M"
  )
  
  targets$PMI <- unlist(pmi_info[targets$Donor])
  targets$Age <- unlist(age_info[targets$Donor])
  targets$Sex <- unlist(sex_info[targets$Donor]) %>% as.factor()
  
  return(targets)
}

GetMetaDataTdpStratified <- function(expr){
  # Assign conditions - Load in metadata (copied and pasted from the supplemental #3)
  short_name_str <- "4C11	4C9	4C8	4C7	4C6	4C5	4C4	4C3	4B11	4B10	4B9	4B8	4B7	4B5	4B4	4B3	4A11	4A10	4A9	4A8	4A7	4A6	4A5	4A4	4A3	1A8	3C5	3B10	3B9	3A10	3A8	3A6	3A5	3A4	2B5	1C9	1C4	1B10	1B6	1A9	1A7	1A5	3C10	3B11	3B7	3B6	3A11	2C10	2C9	2C8	2C7	2B11	2B10	2B9	2B8	2B3	2A11	2A9	2A8	2A7	2A5	2A4	2A3	1C10	1C5	1B4	1A11	1A10	2C3	1A6	3C9	3C8	3C6	3B8	3B5	3B4	3A7	2C6	2C5	2B7	2B6	1B8	1A3	3C7	3C4	3B3	3A9	3A3	2C4	2A10	2A6	1C8	1C7	1C6	1C3	1B11	1B9	1B7	1B5	1B3"
  short_names <- str_split(short_name_str, "\t") %>% unlist()
  
  condition_str <- "CTL	CTL	CTL	CTL	CTL	CTL	CTL	CTL	CTL	CTL	CTL	CTL	CTL	CTL	CTL	CTL	CTL	CTL	CTL	CTL	CTL	CTL	CTL	CTL	CTL	NON	NON	NON	NON	NON	NON	NON	NON	NON	NON	NON	NON	NON	NON	NON	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MLD	MOD	MOD	MOD	MOD	MOD	MOD	MOD	MOD	MOD	MOD	MOD	MOD	MOD	SEV	SEV	SEV	SEV	SEV	SEV	SEV	SEV	SEV	SEV	SEV	SEV	SEV	SEV	SEV	SEV	SEV	SEV"
  conditions <- str_split(condition_str, "\t") %>% unlist()
  
  donor_id_str <- "CTL #4	CTL #4	CTL #4	CTL #4	CTL #4	CTL #4	CTL #4	CTL #4	CTL #4	CTL #4	CTL #4	CTL #4	CTL #4	CTL #5	CTL #5	CTL #5	CTL #5	CTL #5	CTL #5	CTL #5	CTL #5	CTL #5	CTL #5	CTL #5	CTL #5	ALS #3	ALS #1	ALS #1	ALS #1	ALS #1	ALS #1	ALS #1	ALS #1	ALS #1	ALS #4	ALS #3	ALS #3	ALS #3	ALS #3	ALS #3	ALS #3	ALS #3	ALS #1	ALS #1	ALS #1	ALS #1	ALS #1	ALS #4	ALS #4	ALS #4	ALS #4	ALS #4	ALS #4	ALS #4	ALS #4	ALS #4	ALS #4	ALS #4	ALS #4	ALS #4	ALS #4	ALS #4	ALS #4	ALS #3	ALS #3	ALS #3	ALS #3	ALS #3	ALS #4	ALS #3	ALS #1	ALS #1	ALS #1	ALS #1	ALS #1	ALS #1	ALS #1	ALS #4	ALS #4	ALS #4	ALS #4	ALS #3	ALS #3	ALS #1	ALS #1	ALS #1	ALS #1	ALS #1	ALS #4	ALS #4	ALS #4	ALS #3	ALS #3	ALS #3	ALS #3	ALS #3	ALS #3	ALS #3	ALS #3	ALS #3"
  donor_ids <- str_split(donor_id_str, "\t") %>% unlist()
  
  name_condition_table <- data.frame(
    ShortName = short_names,
    Condition = conditions,
    Disease = c(rep("CTRL", 25), rep("ALS", length(conditions) - 25)),
    Donor = donor_ids
  )
  
  # Get shortened names of columns to match up to the short names from the supplemental data
  column_names <- colnames(expr) %>% strsplit("_") %>% sapply(function(x) x[3]) %>% unlist()
  
  # remove metadata that doesn't correspond to any files
  name_condition_table <- name_condition_table[name_condition_table$ShortName %in% column_names,]
  
  # order the new metadata table
  name_condition_table <- name_condition_table[match(column_names, name_condition_table$ShortName), ]
  name_condition_table <- name_condition_table[!is.na(name_condition_table$ShortName), ]
  
  # Metadata (Group + Donor) 
  targets <- data.frame(
    Sample = colnames(expr)[column_names %in% name_condition_table$ShortName],
    Group = factor(name_condition_table$Disease), 
    TDP_Lvl = factor(name_condition_table$Condition, levels = c("CTL", "NON", "MLD", "MOD", "SEV")),
    Donor = factor(name_condition_table$Donor)
  )
  
  pmi_info <- list(
    "ALS #1" = 9.3,
    "ALS #3" = 12,
    "ALS #4" = 28.87, 
    "CTL #4" = 21,
    "CTL #5" = 13
  )
  
  # Attempt to control for PMI
  age_info <- list(
    "ALS #1" = 80,
    "ALS #3" = 66,
    "ALS #4" = 78, 
    "CTL #4" = 46,
    "CTL #5" = 50
  )
  
  # Attempt to control for PMI
  sex_info <- list(
    "ALS #1" = "M",
    "ALS #3" = "F",
    "ALS #4" = "M", 
    "CTL #4" = "F",
    "CTL #5" = "M"
  )
  
  targets$PMI <- unlist(pmi_info[targets$Donor])
  targets$Age <- unlist(age_info[targets$Donor])
  targets$Sex <- unlist(sex_info[targets$Donor]) %>% as.factor()
  
  return(targets)
}

# Custom Normalization for PTM Occupancy ----

# Finds control peptides based on those peptides that are most consistently correlated with the 
# Abundance of other peptides from that protein
GetControlPeptides <- function(y.peptide, pg, mod_seq, num_peps_to_avg = 3)
{
  # Select all peptides from the given protein group
  pep <- y.peptide[y.peptide$genes$ProteinGroup == pg, ]
  if(nrow(pep) < num_peps_to_avg) return(NULL)
  
  # Remove the modified peptide
  pep <- pep[!pep$genes$PeptideSequence == mod_seq, ]
  
  # Find the average for all peptide expression
  avg_exp <- colMeans(pep$E)
  
  # find the degree to which different peptides expression correlates with average expression
  corrs <- apply(pep$E, MARGIN = 1, FUN = function(x) {cor(x, avg_exp, method = "pearson")})
  corrs <- corrs[order(corrs, decreasing = T)]
  if(sum(corrs > 0.8) > num_peps_to_avg){
    top_rows_by_corr <- names(corrs)[corrs > 0.8]
  } else {
    top_rows_by_corr <- names(corrs)[1:num_peps_to_avg]
  }

  # Take the top N most correlated peptides and return them
  pep <- pep[row.names(pep$E) %in% top_rows_by_corr, ]
  
  return(pep)
}

# Create a new eList object where modified peptide intensities are normalized ----
NormalizeExpressionToPeptides <- function(y.peptide, mod_seq, base_seq, pg)
{
  #Find the modified entry and it's expression levels
  modified_row <- which(y.peptide$genes$PeptideSequence == mod_seq)
  modified_expr <- y.peptide$E[modified_row,] 
  
  # Find control peptides with stable expression and high observability
  if(is.null(pg)) {
    modified_expr <- 0
    return(modified_expr)
  }
  control_peps <- GetControlPeptides(y.peptide, pg, mod_seq)
  if(is.null(control_peps)) {
    modified_expr <- 0
    return(modified_expr)
  }
  
  unmodified_expr <- colMeans(control_peps$E)
  
  stopifnot(all(y.peptide$targets$Sample == control_peps$targets$Sample))
  
  # Normalize the expression values
  # Normalization procedure: find the ratio of intensities between the 
  # modified peptidoform and the control peptides (modified_expr - unmodified_expr), 
  # then multiply by a constant scaling factor (median(unmodified_expr))
  normed_expr <- median(unmodified_expr) + (modified_expr - unmodified_expr)
  return(normed_expr)
}
# Example use case
# # Normalize expression/calculate PTM Occupancy, then run eBayes
# all.mod.normed.pep <- GetNormalizedEList(all.mod.peptides, y.peptide)

NormalizeExpressionToProteins <- function(peptide_e, protein_e, mod_seq, pg)
{
  #Find the modified entry and it's expression levels
  modified_row <- which(peptide_e$genes$PeptideSequence == mod_seq)
  modified_expr <- peptide_e$E[modified_row, ] 
  
  # Find protein corresponding to modified peptide
  if(is.null(pg)) {
    modified_expr <- 0
    return(modified_expr)
  }
  protein_row <- which(protein_e$genes$ProteinGroup == pg)
  if(length(protein_row) == 0){
    modified_expr <- 0
    return(modified_expr)
  }
  protein_expr <- protein_e$E[protein_row, ] %>% as.vector()
  
  stopifnot(all(peptide_e$targets$Sample == y.protein$targets$Sample))
  
  # Normalize the expression values
  # Normalization procedure: find the ratio of intensities between the 
  # modified peptidoform and the parent protein (modified_expr - protein_expr ), 
  # then multiply by a constant scaling factor (median(protein_expr ))
  normed_expr <- median(protein_expr) + (modified_expr - protein_expr)
  return(normed_expr)
}

GetNormalizedEList <- function(peptides_of_interest, y.peptide, y.protein = NULL, normalize_to_protein = F){
  
  mod.peptide <- y.peptide[y.peptide$genes$PeptideSequence %in% peptides_of_interest, ]
  if(!normalize_to_protein){
    normed_expr = apply(mod.peptide$genes, MARGIN = 1, FUN = function(x) 
      NormalizeExpressionToPeptides(y.peptide, x[["PeptideSequence"]], x[["BaseSequence"]], x[["ProteinGroup"]]))
  }else if(!is.null(y.protein))
  {
    filtered_protein <- y.protein[y.protein$genes$NPeptides >= 4, ]
    normed_expr = apply(mod.peptide$genes, MARGIN = 1, FUN = function(x) 
      NormalizeExpressionToProteins(y.peptide, filtered_protein, x[["PeptideSequence"]], x[["ProteinGroup"]]))
  }
  
  # If there are any null values, the apply functions above will return a list
  # if there are no null values, its a matrix
  # this if/ifelse clause handles both data types
  if(is.list(normed_expr)){
    filtered_normed_expr <- normed_expr[!is.null(normed_expr)]
    normed_elist = do.call(rbind, filtered_normed_expr)
  } else if (is.numeric(normed_expr)){
    normed_elist <- t(normed_expr)
  }
  
  normed.peptide <- mod.peptide
  stopifnot((normed.peptide$E %>% dim()) == (normed_elist %>% dim()))
  normed.peptide$E <- normed_elist 
  # Filter out rows that could not be normalized
  normed.peptide <- normed.peptide[!apply(normed_elist, MARGIN = 1, FUN = function(x) all(x==0) ), ]
  return(normed.peptide)
}

GetCtrlNormalizedOccupancy <- function(normed.peptide){
  # Hierarchical Clustering ----
  ctrl_mean_occupancies <- normed.peptide$E[,normed.peptide$targets$Group == "CTRL"] %>% rowMeans()
  normed_occupancies <- normed.peptide$E - ctrl_mean_occupancies
  
  row.names(normed_occupancies) <- normed.peptide$genes$PeptideSequence
  normed.peptide$targets$TDP_Lvl <- factor(normed.peptide$targets$TDP_Lvl, levels = c("CTL", "NON", "MLD", "MOD", "SEV"))
  colnames(normed_occupancies) <- normed.peptide$targets$TDP_Lvl
  normed_occupancies <- normed_occupancies[,order(normed.peptide$targets$TDP_Lvl)]
  
  return(normed_occupancies)
}


# Functions for parsing modifications ----

biological_mods <- c("Phosphorylation", "Acetylation", "Citrullination", 
                     "Hydroxyproline", "Methylation", "Diphthamide",
                     "Trimethylation", "Succinylation", "Dimethylation",
                     "Glutarylation", "Palmitoylation",
                     "Succinylation", "Malonylation",  "Crotonylation",
                     "Butyrylation", "Hydroxybutyrylation") #, "Carboxymethyllysine")

biological_mod_regex <- paste0(biological_mods, collapse = "|")


uninteresting_mods <- c("Carbamidomethyl", "Carbamyl", "Formylation",
                        "Sodium", "Potassium", "Magnesium", "Calcium",
                        "Deamidation", "Deamidated", "Fe[III", "Iron",
                        "Ammonia", "Water", "Water Loss", "Carboxylation",
                        "Nitrosylation")

# Takes in a vector of peptide full sequences and returns a named list
# where every element is a character vector of the mods contained in the given
# peptide
GetMods <- function(peptides)
{
  mods <- str_match_all(peptides, "\\[.*?:([^\\]]+)\\]") %>% lapply(., function(x) x[,2])
  names(mods) <- peptides
  return( mods )
}

# Takes in a list of peptide full sequences and returns a named list
GetModsTextWReplacement <- function(peptides, all_mods = T)
{
  mods_to_remove <- uninteresting_mods
  if(all_mods) {
    mods_to_remove <- c()
  }
  mods <- GetMods(peptides)
  replacedMods <- lapply(mods, function(x) lapply(x, .replaceMods) %>%
                           unlist()  %>%
                           .[!. %in% mods_to_remove] %>%
                           paste0(., collapse = ", "))
  return(replacedMods)
}

# This function replaces modification substrings with more general terms
# It removes all instances of carbamidomethyl on C and keeps carbamidomethyl elsewhere
.replaceMods <- function(x)
{
  # Define replacements as a named vector:
  #   names = substrings to look for
  #   values = replacement words
  replacements <- c(
    "Phospho" = "Phosphorylation",
    "Citrulli" ="Citrullination",
    "acetyl"="Acetylation",
    "Acetylation" = "Acetylation",
    "xidation on P"="Hydroxyproline",
    "Hydroxylation on P" = "Hydroxyproline",
    "4-hydroxypro" = "Hydroxyproline",
    "Hydroxylation" = "Oxidation",
    "xidation" = "Oxidation",
    "hydroxy" = "Oxidation",
    "GG" = "Ubiquitination",
    "Omega-N" = "Methylation",
    "Tele-methly" = "Methylation",
    "-methyl" = "Methylation",
    "trimethyl" = "Trimethylation",
    "dimethyl" = "Dimethylation",
    "Symmetric" = "Dimethylation",
    "Asymmetric" = "Dimethylation",
    "Dimethylated" = "Dimethylation",
    "-nitro" = "Nitrosylation",
    "-succin" = "Succinylation",
    "-malonyl" = "Malonylation",
    "-glutamyl" = "Glutamylation",
    "butyryll" = "Butyrylation",    
    "palmitoyl" = "Palmitoylation",
    "myristoyl" = "Myristoylation",
    "glutary" = "Glutarylation",
    "crotonyl" = "Crotonylation",
    "lipoyll" = "Lipoylation",
    "Pyrrolidone" = "Pyroglutamate",
    "Cysteine" = "Cysteine Sulfinic Acid",
    "N-Pyruvate" = "N-pyruvate 2-iminyl-valine",
    "Fe\\[" = "Iron",
    "Cu\\[" = "Copper",
    "Deamidat" = "Deamidation",
    "Water" = "Water Loss",
    "Carboxymethylation on K" = "Carboxymethyllysine"
  )
  for (pattern in names(replacements)) {
    if (grepl(pattern, x)) {
      x <- replacements[pattern]
      break # stops at first matching pattern
    }
  }
  # Remove carbamidomethyl on cysteine
  x <- x[x != "Carbamidomethyl on C"]
  # Trim off everything after the first whitespace
  x <- gsub("\\s.*", "", x)
  return(x)
}

# Classify peptides as either unmodified, biologically modified, or "other"
# Fucntion depends on the biological_mods vector defined above
ClassifyPeptideMods <- function(replacedMods)
{
  classifications <- sapply(replacedMods, function(x) {
    if(x == "") {
      return("Unmodified")
    } else if(grepl(biological_mod_regex, x)) {
      return("Enzymatic Mod")
    } else {
      return("Non-Enzymatic Mod")
    }
  })
  return(classifications)
}

ClassifyPeptideModsSpecific <- function(replacedMods)
{
  classifications <- sapply(replacedMods, function(x) {
    if(x == "") {
      return("Unmodified")
    } else if(grepl(biological_mod_regex, x)) {
      return("Enzymatic Mod")
    } else if(AllCarbamidomethylation(x)) {
      return("Carbamidomethylation")
    } else if(AllCarboxymethyl(x)) {
      return("Carboxymethylation") 
    } else {
      return("Non-Enzymatic Mod")
    }
  })
  return(classifications)
}

AllCarbamidomethylation <- function(modstring)
{
  mods <- str_split(modstring, ", ") %>% unlist()
  if (all(mods == "Carbamidomethyl")){
    return (T)
  } else {
    return (F)
  }
}

AllCarboxymethyl<- function(modstring)
{
  mods <- str_split(modstring, ", ") %>% unlist()
  if (all(mods == "Carboxymethyllysine" | mods == "Carboxymethylation" | mods == "Carbamidomethyl")){
    return (T)
  } else {
    return (F)
  }
}

# Returns TRUE if a single raw mod string (from GetMods) is present in a
# limited-PTM search: Carbamidomethyl on Cys (fixed), Oxidation on Met
# (variable), or Deamidation on Asn/Gln (variable).
# Must operate on raw strings BEFORE .replaceMods, which collapses residue info
# (e.g. "Hydroxylation on N" and "Oxidation on M" both become "Oxidation").
.isLimitedSearchRawMod <- function(m) {
  grepl("xidation on M$", m) ||
    (grepl("[Dd]eamidat", m) && grepl(" on [NQ]$", m)) ||
    m == "Carbamidomethyl on C"
}

# Returns a logical vector the same length as raw_mods_list (named list from
# GetMods()): TRUE if every modification in that peptide is limited-search
# compatible. Unmodified peptides (empty mod vector) return FALSE.
IsLimitedSearchPeptide <- function(raw_mods_list) {
  sapply(raw_mods_list, function(mods) {
    length(mods) > 0 && all(vapply(mods, .isLimitedSearchRawMod, logical(1)))
  })
}

ClassifyPeptideModsSpecificPrivilegeCarboxymethyl <- function(replacedMods)
{
  classifications <- sapply(replacedMods, function(x) {
    if(x == "") {
      return("Unmodified")
    } else if(grepl("Carboxymethyl", x)) {
      return("Carboxymethylation") 
    } else if(grepl(biological_mod_regex, x)) {
      return("Enzymatic Mod")
    } else if(AllCarbamidomethylation(x)) {
      return("Carbamidomethylation")
    }  else {
      return("Non-Enzymatic Mod")
    }
  })
  return(classifications)
}

ClassifyPeptideModsSpecificIgnoreCarboxymethyl <- function(replacedMods)
{
  classifications <- sapply(replacedMods, function(x) {
    if(x == "") {
      return("Unmodified")
    } else if(grepl(biological_mod_regex, x)) {
      return("Enzymatic Mod")
    } else if(AllCarbamidomethylation(x)) {
      return("Carbamidomethylation")
    } else {
      return("Non-Enzymatic Mod")
    }
  })
  return(classifications)
}

# Takes in a vector of peptide full sequences and returns a named list
# of bools that is True if the peptide sequence contains only Carbamamidomethylation
# or if the sequence contains no modifications at all
# (A fixed modification where cysteine is alkylated, a common step in sample prep)
# And false if any other mod is present
OnlyCysteineAlkylation <- function(peptides)
{
  mods <- GetMods(peptides)
  all_fixed = sapply(mods, function(x) all(x == "Carbamidomethyl") )
  names(all_fixed) <- names(mods)
  return(all_fixed)
}


# Takes in a vector of peptide full sequences and returns a named list
# of bools that is True if the peptide sequence contains only Carbamamidomethylation or Formylation
# or if the sequence contains no modifications at all
# (A fixed modification where cysteine is alkylated, a common step in sample prep)
# And false if any other mod is present
OnlyCysteineAlkylationOrFormylation <- function(peptides)
{
  mods <- GetMods(peptides)
  all_fixed = sapply(mods, function(x) all(x %in% c("Carbamidomethyl", "Formylation") ))
  names(all_fixed) <- names(mods)
  return(all_fixed)
}

# Takes in a vector of peptide full sequences and returns a named list
# of bools that is True if the peptide sequence contains only Carbamamidomethylation or Formylation
# or if the sequence contains no modifications at all
# (A fixed modification where cysteine is alkylated, a common step in sample prep)
# And false if any other mod is present
OnlyCysteineAlkylationOrFormylationOrVariableOxidationOfMethionine <- function(peptides)
{
  mods <- GetMods(peptides)
  all_fixed = sapply(mods, function(x) all(x %in% c("Carbamidomethyl", "Formylation", "Oxidation") ))
  names(all_fixed) <- names(mods)
  return(all_fixed)
}

# Returns TRUE if peptide has ONLY non-enzymatic modifications (should be excluded from bio analysis)
# Oxidation on P (hydroxyproline) is biological; oxidation on other residues is artifact
# Peptides with co-occurring non-enzymatic mods (e.g., Oxidation on M + Formylation) will be excluded
OnlyNonEnzymaticMods <- function(peptides) {
  non_enzymatic <- c("Carbamidomethyl", "Formylation", "Other")
  non_enzymatic <- c(non_enzymatic, uninteresting_mods)

  sapply(peptides, function(pep) {
    # Get raw modification strings (includes residue info like "Oxidation on M")
    mods_raw <- str_match_all(pep, "\\[.*?:([^\\]]+)\\]")[[1]][,2]

    if (length(mods_raw) == 0) return(TRUE)  # No mods = exclude

    # Check each mod
    all_non_enzymatic <- all(sapply(mods_raw, function(mod) {
      mod_name <- gsub(" on .*", "", mod)  # Extract just the mod type

      # Oxidation/Hydroxylation on P is biological (hydroxyproline)
      if (grepl("xidation|Hydroxylation", mod) && grepl("on P", mod)) {
        return(FALSE)  # This is biological, not non-enzymatic
      }
      # Other oxidations are non-enzymatic artifacts
      if (grepl("xidation|Hydroxylation", mod)) {
        return(TRUE)  # Non-enzymatic
      }
      # Check against known non-enzymatic list
      return(mod_name %in% non_enzymatic)
    }))

    return(all_non_enzymatic)
  })
}

# Finding sequence positions ----

# The SequencePositionTable.tsv file was generated by parsing the full sequences of all peptides in the dataset
# The source file can be found here: r"(D:\Kelly_ALS_motor_nueron_dataset\MM1p1p4_GPTMD_Search_Carboxymethyl_Carbamido_2\Task2-SearchTask\AllPeptides.psmtsv)"

GetSequencePosition <- function(base_sequence){
  if(!exists("sequencePositionTable")){
    path <- "Data/SequencePositionTable.tsv"
    sequencePositionTable <- read.csv(path, sep = '\t', row.names = NULL)
  }
  return(sequencePositionTable[sequencePositionTable$Base.Sequence == base_sequence, c("Start", "End")])
}

ParsePosition <- function(position_in_brackets){
  matches <- regmatches(position_in_brackets, regexec("\\[(\\d+)\\s+to\\s+(\\d+)\\]", position_in_brackets))[[1]]
  
  first_num  <- as.numeric(matches[2])
  second_num <- as.numeric(matches[3])
  return(c(Start = first_num, End = second_num))
}

GetPositionVector <- function(base_sequence){
  GetSequencePosition(base_sequence) %>% ParsePosition() %>% return()
}

AddPositionColumns <- function(df){
  if (!exists("sequencePositionTable")) {
    sequencePositionTable <<- read.csv("Data/SequencePositionTable.tsv",
                                       sep = '\t', row.names = NULL)
  }
  idx <- match(df$BaseSequence, sequencePositionTable$Base.Sequence)
  df$Start <- sequencePositionTable$Start[idx]
  df$End   <- sequencePositionTable$End[idx]
  return(df)
}

# Find Opposing peptides ----

# Expects two rows and two columns labeled "Start" and "End"
IsOverlappingTwoTable <- function(df1, df2){
  if(df1$Gene[1] != df2$Gene[1]) return(F)
  start1 <- df1[1, "Start"]
  start2 <- df2[1, "Start"]
  end1 <- df1[1, "End"]
  end2 <- df2[1, "End"]
  
  if(length(end2) == 0 | length(end1) == 0) return(F)
  if(length(start2) == 0 | length(start1) == 0) return(F)
  
  if(start1 >= start2 & start1 <= end2) return(T)
  if(start2 >= start1 & start2 <= end1) return(T)
  return(F)
}

# Takes in a table with the following columns: Gene, Start, End, logFC, PeptideSequence
CheckForOpposingPeptidesTwoTable <- function(queryTable, dbTable){
  opposing <- rep(F, nrow(queryTable))
  names(opposing) <- queryTable$PeptideSequence
  
  for(gene in unique(queryTable$Gene))
  {
    singleGeneQueryDf <- queryTable[grepl(gene, queryTable$Gene), ]
    singleGeneDb <- dbTable[grepl(gene, dbTable$Gene), ]
    
    for(i in 1:nrow(singleGeneQueryDf)){
      for(j in 1:nrow(singleGeneDb)){
        if( IsOverlappingTwoTable(singleGeneQueryDf[i,], singleGeneDb[j,] ) &
            sign(singleGeneQueryDf[i, "logFC"]) !=  sign(singleGeneDb[j, "logFC"]) ){
          opposing[singleGeneQueryDf$PeptideSequence[i]] <- T
        }
      }
    }
  }
  
  return(opposing)
}

# TODO: Function that extend check for opposing to return full opposing table

# Expects two rows and two columns labeled "Start" and "End"
IsOverlapping <- function(df){
  if(df$Gene[1] != df$Gene[2]) return(F)
  start1 <- df[1, "Start"]
  start2 <- df[2, "Start"]
  end1 <- df[1, "End"]
  end2 <- df[2, "End"]
  
  if(start1 >= start2 & start1 <= end2) return(T)
  if(start2 >= start1 & start2 <= end1) return(T)
  return(F)
}

# Takes in a table with the following columns: Gene, Start, End, logFC, PeptideSequence
CheckForOpposingPeptides <- function(table){
  opposing <- rep(F, nrow(table))
  names(opposing) <- table$PeptideSequence
  
  for(gene in unique(table$Gene))
  {
    singleGeneDf <- table[table$Gene == gene, ]
    
    for(i in 1:nrow(singleGeneDf)){
      for(j in 1:nrow(singleGeneDf)){
        if( IsOverlapping(singleGeneDf[c(i,j), ]) & sign(singleGeneDf[i, "logFC"]) !=  sign(singleGeneDf[j, "logFC"]) ){
          opposing[singleGeneDf$PeptideSequence[i]] <- T
          opposing[singleGeneDf$PeptideSequence[j]] <- T
        }
      }
    }
  }
  
  return(opposing)
}
