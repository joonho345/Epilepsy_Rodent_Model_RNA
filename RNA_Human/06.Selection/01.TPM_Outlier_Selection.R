#############################
# Outlier Selection Script
# z_scores_mtle < -2.0 (low expression in MTLE)
# z_scores_nl > 2.0 (high expression in NL)
#############################

library(dplyr)

###############################################################
##################### Use TPM Count matrix ####################
###############################################################

#### Read TPM matrix ####
tpm_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/03.Normalization/adjusted_merged_matrix_1_TPM.txt"
Target_matrix <- read.table(tpm_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

#### Read coldata ####
coldata_df <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/02.Quantification/filtered_coldata.txt",
                         sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE, 
                         fill = TRUE, quote = "", comment.char = "")

#### Apply filtering logic (matching 00.Subgrouping_H.R lines 59-69) ####
# Filter to MTLEHS, MTLE, NL and create MTLEALL grouping
Target_coldata <- coldata_df[coldata_df$Diagnosis %in% c("MTLEHS", "MTLE", "NL"), ]
Target_coldata$Diagnosis_group <- ifelse(Target_coldata$Diagnosis %in% c("MTLEHS", "MTLE"), "MTLEALL", "NL")
Target_coldata$Diagnosis <- NULL
colnames(Target_coldata)[which(names(Target_coldata) == "Diagnosis_group")] <- "Diagnosis"
Target_coldata <- Target_coldata %>% 
  filter(Brain_Location == "HIPPO") %>% 
  mutate(across(everything(), ~ replace(.x, .x == "", "N/A"))) %>%
  mutate(Diagnosis_Sub = ifelse(Diagnosis_Sub == "TLE", "MTLE", Diagnosis_Sub))

# Ensure samples match
common_samples <- intersect(rownames(Target_coldata), colnames(Target_matrix))
Target_matrix <- Target_matrix[, common_samples, drop = FALSE]
Target_coldata <- Target_coldata[common_samples, , drop = FALSE]

# Add sample_name column for merging
Target_coldata$sample_name <- rownames(Target_coldata)

specified_groups <- c("MTLEHS", "MTLE", "NL")
Target_groups <- 'HIPPO'

# Select genes
genes <- c("GFAP", "NOTCH1", "AQP4")

# Initialize empty lists to store outliers
mtle_outliers <- list()
nl_outliers <- list()

# Loop through each gene
for (gene in genes) {
  # Extract expression data for the current gene
  gene_expression <- Target_matrix[gene, ]
  gene_expression <- t(gene_expression)
  gene_expression <- as.data.frame(gene_expression)
  gene_expression$sample_name <- rownames(gene_expression)
  colnames(gene_expression)[1] <- "target_gene_expression"
  
  # Merge with metadata
  gene_expression_merged <- merge(gene_expression, Target_coldata, by = "sample_name")
  
  # Filter for specified groups
  gene_expression_filtered <- gene_expression_merged[gene_expression_merged$Diagnosis %in% specified_groups, ]
  gene_expression_filtered$Diagnosis <- factor(gene_expression_filtered$Diagnosis, levels = specified_groups)
  
  # Split data into MTLEHS + MTLE and NL groups
  mtle_data <- gene_expression_filtered[gene_expression_filtered$Diagnosis %in% c("MTLEHS", "MTLE"), ]
  nl_data <- gene_expression_filtered[gene_expression_filtered$Diagnosis == "NL", ]
  
  # Calculate Z-scores for MTLEHS + MTLE
  if (nrow(mtle_data) > 0 && sd(mtle_data$target_gene_expression) > 0) {
    z_scores_mtle <- (mtle_data$target_gene_expression - mean(mtle_data$target_gene_expression)) / sd(mtle_data$target_gene_expression)
    outliers_mtle <- mtle_data$sample_name[which(z_scores_mtle < -2.0)]
  } else {
    outliers_mtle <- character(0)
  }
  
  # Calculate Z-scores for NL
  if (nrow(nl_data) > 0 && sd(nl_data$target_gene_expression) > 0) {
    z_scores_nl <- (nl_data$target_gene_expression - mean(nl_data$target_gene_expression)) / sd(nl_data$target_gene_expression)
    outliers_nl <- nl_data$sample_name[which(z_scores_nl > 2.0)]
  } else {
    outliers_nl <- character(0)
  }
  
  # Store outliers in lists
  mtle_outliers[[gene]] <- outliers_mtle
  nl_outliers[[gene]] <- outliers_nl
}

# Create output directory
output_dir <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/06.Selection/01.TPM_Outliers"

# Common MTLE outliers (intersection of all three genes)
if (all(genes %in% names(mtle_outliers))) {
  common_mtle_outliers <- Reduce(intersect, list(
    mtle_outliers[[genes[1]]],
    mtle_outliers[[genes[2]]],
    mtle_outliers[[genes[3]]]
  ))
} else {
  common_mtle_outliers <- character(0)
}

# Common NL outliers (intersection of all three genes)
if (all(genes %in% names(nl_outliers))) {
  common_nl_outliers <- Reduce(intersect, list(
    nl_outliers[[genes[1]]],
    nl_outliers[[genes[2]]],
    nl_outliers[[genes[3]]]
  ))
} else {
  common_nl_outliers <- character(0)
}

# Save MTLE outliers for each gene
for (gene in genes) {
  if (gene %in% names(mtle_outliers)) {
    file_name <- paste0(output_dir, "/01.MTLE_Outliers_", gene, ".txt")
    write.table(mtle_outliers[[gene]], file = file_name, 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

# Save NL outliers for each gene
for (gene in genes) {
  if (gene %in% names(nl_outliers)) {
    file_name <- paste0(output_dir, "/01.NL_Outliers_", gene, ".txt")
    write.table(nl_outliers[[gene]], file = file_name, 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

# Save common outliers
if (length(common_mtle_outliers) > 0) {
  file_name <- paste0(output_dir, "/01.MTLE_Outliers_Common_All_Genes.txt")
  write.table(common_mtle_outliers, file = file_name, 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

if (length(common_nl_outliers) > 0) {
  file_name <- paste0(output_dir, "/01.NL_Outliers_Common_All_Genes.txt")
  write.table(common_nl_outliers, file = file_name, 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Union of all MTLE outliers across all genes
if (all(genes %in% names(mtle_outliers))) {
  all_mtle_outliers <- unique(c(
    mtle_outliers[[genes[1]]],
    mtle_outliers[[genes[2]]],
    mtle_outliers[[genes[3]]]
  ))
} else {
  all_mtle_outliers <- character(0)
}

# Union of all NL outliers across all genes
if (all(genes %in% names(nl_outliers))) {
  all_nl_outliers <- unique(c(
    nl_outliers[[genes[1]]],
    nl_outliers[[genes[2]]],
    nl_outliers[[genes[3]]]
  ))
} else {
  all_nl_outliers <- character(0)
}

# Combined list of all outliers (MTLE + NL)
ALL_outliers <- unique(c(all_mtle_outliers, all_nl_outliers))

if (length(all_mtle_outliers) > 0) {
  file_name <- paste0(output_dir, "/01.MTLE_Outliers_All_Genes.txt")
  write.table(all_mtle_outliers, file = file_name, 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

if (length(all_nl_outliers) > 0) {
  file_name <- paste0(output_dir, "/01.NL_Outliers_All_Genes.txt")
  write.table(all_nl_outliers, file = file_name, 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

if (length(ALL_outliers) > 0) {
  file_name <- paste0(output_dir, "/01.ALL_Outliers_All_Genes.txt")
  write.table(ALL_outliers, file = file_name, 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
