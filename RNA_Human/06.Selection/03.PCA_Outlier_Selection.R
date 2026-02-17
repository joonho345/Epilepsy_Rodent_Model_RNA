#############################
# PCA Outlier Selection Script
# MTLE outliers: PC2 < 0
# NL outliers: PC2 > 0
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

# Ensure sample order matches
if (!all(colnames(Target_matrix) == rownames(Target_coldata))) {
  Target_coldata <- Target_coldata[colnames(Target_matrix), , drop = FALSE]
}

# Add sample_name column for merging
Target_coldata$sample_name <- rownames(Target_coldata)

###############################################################
##################### PCA Analysis ####################
###############################################################

# Filter to only H_gene_set_W genes (matching the plotting script)
genes_in_matrix <- rownames(Target_matrix)
genes_to_keep <- intersect(genes_in_matrix, H_gene_set_W)
tpm_matrix_filtered <- Target_matrix[genes_to_keep, , drop = FALSE]

# Filter samples
filtered_samples <- intersect(rownames(Target_coldata), colnames(tpm_matrix_filtered))
tpm_matrix_filtered <- tpm_matrix_filtered[, filtered_samples, drop = FALSE]
coldata_subset <- Target_coldata[filtered_samples, , drop = FALSE]

# Ensure sample order matches
if (!all(colnames(tpm_matrix_filtered) == rownames(coldata_subset))) {
  coldata_subset <- coldata_subset[colnames(tpm_matrix_filtered), , drop = FALSE]
}

#### Perform PCA ####
pca_data <- as.matrix(tpm_matrix_filtered)
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)

summary_df <- summary(pca_result)$importance
variance_explained_pct <- round(summary_df["Proportion of Variance", ] * 100, 2)

# Create PCA data frame using $rotation (sample coordinates when genes are rows)
pca_df <- as.data.frame(pca_result$rotation)
pca_df$Sample <- rownames(pca_df)

# Add metadata
pca_df$Diagnosis <- coldata_subset[pca_df$Sample, "Diagnosis"]
if ("Diagnosis_Sub" %in% colnames(coldata_subset)) {
  pca_df$Diagnosis_Sub <- coldata_subset[pca_df$Sample, "Diagnosis_Sub"]
}

###############################################################
##################### Identify Outliers ####################
###############################################################

# Split data into MTLE (MTLEHS + MTLE) and NL groups
mtle_data <- pca_df[pca_df$Diagnosis == "MTLEALL", ]
nl_data <- pca_df[pca_df$Diagnosis == "NL", ]

# MTLE outliers: PC2 < 0
mtle_outliers <- mtle_data$Sample[mtle_data$PC2 < 0]

# NL outliers: PC2 > 0
nl_outliers <- nl_data$Sample[nl_data$PC2 > 0]

# Combined list of all outliers
ALL_outliers <- unique(c(mtle_outliers, nl_outliers))

###############################################################
##################### Save Results ####################
###############################################################

# Create output directory
output_dir <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/06.Selection/03.PCA_Outliers"

# Save MTLE outliers
if (length(mtle_outliers) > 0) {
  file_name <- paste0(output_dir, "/03.MTLE_Outliers_PCA.txt")
  write.table(mtle_outliers, file = file_name, 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Save NL outliers
if (length(nl_outliers) > 0) {
  file_name <- paste0(output_dir, "/03.NL_Outliers_PCA.txt")
  write.table(nl_outliers, file = file_name, 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Save all outliers
if (length(ALL_outliers) > 0) {
  file_name <- paste0(output_dir, "/03.ALL_Outliers_PCA.txt")
  write.table(ALL_outliers, file = file_name, 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
