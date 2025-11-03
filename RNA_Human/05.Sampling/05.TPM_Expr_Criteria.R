#############################
# z_scores_mtle < -1.5
# z_scores_nl > 1.5
#############################

###############################################################
##################### Use TPM Count matrix ####################
###############################################################

# choose target groups (hippo)
Target_matrix <- adjusted_df_hippo
Target_coldata <- coldata_df_hippo
specified_groups <- c("MTLEHS", "MTLE", "NL")
Target_groups <- 'HIPPO'

# Select genes
genes <- c("GFAP", "NOTCH1")

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
  gene_expression_merged <- merge(gene_expression, Target_coldata, by.x = "sample_name", by.y = "Run")
  
  # Filter for specified groups
  gene_expression_filtered <- gene_expression_merged[gene_expression_merged$Diagnosis %in% specified_groups, ]
  gene_expression_filtered$Diagnosis <- factor(gene_expression_filtered$Diagnosis, levels = specified_groups)
  
  # Split data into MTLEHS + MTLE and NL groups
  mtle_data <- gene_expression_filtered[gene_expression_filtered$Diagnosis %in% c("MTLEHS", "MTLE"), ]
  nl_data <- gene_expression_filtered[gene_expression_filtered$Diagnosis == "NL", ]
  
  # Calculate Z-scores for MTLEHS + MTLE
  z_scores_mtle <- (mtle_data$target_gene_expression - mean(mtle_data$target_gene_expression)) / sd(mtle_data$target_gene_expression)
  outliers_mtle <- mtle_data$sample_name[which(z_scores_mtle < -1.5)]
  
  # Calculate Z-scores for NL
  z_scores_nl <- (nl_data$target_gene_expression - mean(nl_data$target_gene_expression)) / sd(nl_data$target_gene_expression)
  outliers_nl <- nl_data$sample_name[which(z_scores_nl > 1.5)]
  
  # Store outliers in lists
  mtle_outliers[[gene]] <- outliers_mtle
  nl_outliers[[gene]] <- outliers_nl
}

# Identify shared and union outliers
shared_mtle_outliers <- intersect(mtle_outliers[[genes[1]]], mtle_outliers[[genes[2]]])
shared_nl_outliers <- intersect(nl_outliers[[genes[1]]], nl_outliers[[genes[2]]])

union_mtle_outliers <- union(mtle_outliers[[genes[1]]], mtle_outliers[[genes[2]]])
union_nl_outliers <- union(nl_outliers[[genes[1]]], nl_outliers[[genes[2]]])

# Output the results
cat("Outliers for GFAP and NOTCH1 in MTLEHS + MTLE:\n")
print(mtle_outliers)

cat("\nOutliers for GFAP and NOTCH1 in NL:\n")
print(nl_outliers)

cat("\nShared outliers for GFAP and NOTCH1 in MTLEHS + MTLE:\n")
print(shared_mtle_outliers)

cat("\nShared outliers for GFAP and NOTCH1 in NL:\n")
print(shared_nl_outliers)

cat("\nUnion of outliers for GFAP and NOTCH1 in MTLEHS + MTLE:\n")
print(union_mtle_outliers)

cat("\nUnion of outliers for GFAP and NOTCH1 in NL:\n")
print(union_nl_outliers)

ALL_outliers <- c(union_mtle_outliers, union_nl_outliers)
ALL_outliers
