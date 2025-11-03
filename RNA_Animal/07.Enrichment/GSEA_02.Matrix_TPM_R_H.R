library(readr)
library(tidyverse)
library(dplyr)

#### import TPM Count matrix with gene id ####
adjusted_TPM_id_df_R <- read.table("/home/joonho345/3_RNA/RNA_Animal/05.Enrichment/GSEA_01.Matrix_GTF/adjusted_merged_matrix_id_TPM_R.txt", 
                                   sep = "\t", header = TRUE, row.names = 1)
#### import coldata ####
coldata_df_R <- read.table("/home/joonho345/3_RNA/RNA_Animal/02.Quantification/filtered_coldata_R.txt",
                           sep = "\t", header = TRUE, row.names = 1, fill = TRUE)


## G: PPS & CTL_PPS
coldata_df_R_G <- coldata_df_R %>% filter(cohort == "G")
filtered_samples_R_G <- rownames(coldata_df_R_G)
adjusted_df_R_G <- adjusted_TPM_id_df_R[, colnames(adjusted_TPM_id_df_R) %in% filtered_samples_R_G]

### Choose one matrix to analyze ###
coldata_df_R_PPS_O_AC <- coldata_df_R_G %>% filter(TYPE_1 == "PPS_O_AC" | TYPE_1 == "PPS_CTL")
filtered_samples_R_PPS_O_AC <- rownames(coldata_df_R_PPS_O_AC)
adjusted_df_R_PPS_O_AC <- adjusted_df_R_G[, colnames(adjusted_df_R_G) %in% filtered_samples_R_PPS_O_AC]
adjusted_df_GSEA <- adjusted_df_R_PPS_O_AC
coldata_df_GSEA <- coldata_df_R_PPS_O_AC
Target_comparison <- 'R_PPS_O_AC'
contrast_vector <- c("Treatment_Lateral","PPS","CTL_PPS")

coldata_df_R_PPS_O_IM <- coldata_df_R_G %>% filter(TYPE_1 == "PPS_O_IM" | TYPE_1 == "PPS_CTL")
filtered_samples_R_PPS_O_IM <- rownames(coldata_df_R_PPS_O_IM)
adjusted_df_R_PPS_O_IM <- adjusted_df_R_G[, colnames(adjusted_df_R_G) %in% filtered_samples_R_PPS_O_IM]
adjusted_df_GSEA <- adjusted_df_R_PPS_O_IM
coldata_df_GSEA <- coldata_df_R_PPS_O_IM
Target_comparison <- 'R_PPS_O_IM'
contrast_vector <- c("Treatment_Lateral","PPS","CTL_PPS")

coldata_df_R_PPS_O_CR <- coldata_df_R_G %>% filter(TYPE_1 == "PPS_O_CR" | TYPE_1 == "PPS_CTL")
filtered_samples_R_PPS_O_CR <- rownames(coldata_df_R_PPS_O_CR)
adjusted_df_R_PPS_O_CR <- adjusted_df_R_G[, colnames(adjusted_df_R_G) %in% filtered_samples_R_PPS_O_CR]
adjusted_df_GSEA <- adjusted_df_R_PPS_O_CR
coldata_df_GSEA <- coldata_df_R_PPS_O_CR
Target_comparison <- 'R_PPS_O_CR'
contrast_vector <- c("Treatment_Lateral","PPS","CTL_PPS")

coldata_df_R_PPS_O_DOFS <- coldata_df_R_G %>% filter(TYPE_1 == "PPS_O_DOFS" | TYPE_1 == "PPS_CTL")
filtered_samples_R_PPS_O_DOFS <- rownames(coldata_df_R_PPS_O_DOFS)
adjusted_df_R_PPS_O_DOFS <- adjusted_df_R_G[, colnames(adjusted_df_R_G) %in% filtered_samples_R_PPS_O_DOFS]
adjusted_df_GSEA <- adjusted_df_R_PPS_O_DOFS
coldata_df_GSEA <- coldata_df_R_PPS_O_DOFS
Target_comparison <- 'R_PPS_O_DOFS'
contrast_vector <- c("Treatment_Lateral","PPS","CTL_PPS")


# Check if coldata and count matrix have the same order with respect to samples
if (!all(rownames(coldata_df_GSEA) == colnames(adjusted_df_GSEA))) {
  stop("Rownames of coldata do not match colnames of count data")
}


###################################################################
#### GSEA Matrix ####
gsea_expr_df <- adjusted_df_GSEA
gsea_coldata_df <- coldata_df_GSEA

#### GSEA Matrix - Expression matrix ####
num_genes <- nrow(gsea_expr_df)
num_datasets <- ncol(gsea_expr_df)
vector_colnames <- colnames(gsea_expr_df)

gsea_expr_df$NAME <- rownames(gsea_expr_df)   # Insert the 'NAME' column with gene IDs
gsea_expr_df$description <- "NA"         # Insert the 'description' column with 'NA'
gsea_expr_df <- gsea_expr_df[, c("NAME", "description", vector_colnames)]

# Create the header lines
header_lines <- c(
  "#1.2",
  paste(num_genes, num_datasets, sep = "\t"),
  paste(colnames(gsea_expr_df), collapse = "\t")
)

# Write the GSEA formatted file
gsea_output_path <- paste0("/home/joonho345/3_RNA/RNA_Animal/05.Enrichment/GSEA_02.Matrix_GSEA/gsea_formatted_matrix_",Target_comparison,".txt")
writeLines(header_lines, con = gsea_output_path)
write.table(gsea_expr_df, file = gsea_output_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)


#### GSEA Matrix - Phenotype matrix ####
diagnosis_vector <- gsea_coldata_df$Treatment_Lateral
num_samples <- length(diagnosis_vector)
num_groups <- length(unique(diagnosis_vector))
first_line <- paste(num_samples, num_groups, 1, sep = " ")

group_names <- unique(diagnosis_vector)
second_line <- paste("#", paste(group_names, collapse = " "), sep = " ")
third_line <- paste(diagnosis_vector, collapse = "\t")

phenotype_matrix <- c(first_line, second_line, third_line)

# Write the phenotype matrix to a file
output_file_path <- paste0("/home/joonho345/3_RNA/RNA_Animal/05.Enrichment/GSEA_02.Matrix_GSEA/phenotype_matrix_",Target_comparison,".txt")
writeLines(phenotype_matrix, con = output_file_path)

