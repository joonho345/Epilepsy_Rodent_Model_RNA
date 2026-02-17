library(readr)
library(tidyverse)
library(dplyr)

#### import TPM Count matrix with gene id ####
adjusted_TPM_id_df_M <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/GSEA_01.Matrix_GTF/adjusted_merged_matrix_id_TPM_M.txt", 
                                 sep = "\t", header = TRUE, row.names = 1)
#### import coldata ####
coldata_df_M <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/02.Quantification/filtered_coldata_M.txt",
                           sep = "\t", header = TRUE, row.names = 1, fill = TRUE)


# M_PILO_IP_A_HA
coldata_df_M_PILO_IP_A_HA <- coldata_df_M %>% filter(TYPE_2 == "M_PILO_IP_A_HA" | TYPE_2 == "CTL")
filtered_samples_M_PILO_IP_A_HA <- rownames(coldata_df_M_PILO_IP_A_HA)
adjusted_df_M_PILO_IP_A_HA <- adjusted_TPM_id_df_M[, colnames(adjusted_TPM_id_df_M) %in% filtered_samples_M_PILO_IP_A_HA]
adjusted_df_GSEA <- adjusted_df_M_PILO_IP_A_HA
coldata_df_GSEA <- coldata_df_M_PILO_IP_A_HA
Target_comparison <- 'M_PILO_IP_A_HA'
contrast_vector <- c("Treatment_Lateral","PILO","CTL_SAL_IP")

# M_PILO_IP_A_AC
coldata_df_M_PILO_IP_A_AC <- coldata_df_M %>% filter(TYPE_2 == "M_PILO_IP_A_AC")
filtered_samples_M_PILO_IP_A_AC <- rownames(coldata_df_M_PILO_IP_A_AC)
adjusted_df_M_PILO_IP_A_AC <- adjusted_TPM_id_df_M[, colnames(adjusted_TPM_id_df_M) %in% filtered_samples_M_PILO_IP_A_AC]
adjusted_df_GSEA <- adjusted_df_M_PILO_IP_A_AC
coldata_df_GSEA <- coldata_df_M_PILO_IP_A_AC
Target_comparison <- 'M_PILO_IP_A_AC'
contrast_vector <- c("Treatment_Lateral","PILO","CTL_SAL_IP")

# M_PILO_IP_A_IM
coldata_df_M_PILO_IP_A_IM <- coldata_df_M %>% filter(TYPE_2 == "M_PILO_IP_A_IM" | TYPE_2 == "CTL")
filtered_samples_M_PILO_IP_A_IM <- rownames(coldata_df_M_PILO_IP_A_IM)
adjusted_df_M_PILO_IP_A_IM <- adjusted_TPM_id_df_M[, colnames(adjusted_TPM_id_df_M) %in% filtered_samples_M_PILO_IP_A_IM]
adjusted_df_GSEA <- adjusted_df_M_PILO_IP_A_IM
coldata_df_GSEA <- coldata_df_M_PILO_IP_A_IM
Target_comparison <- 'M_PILO_IP_A_IM'
contrast_vector <- c("Treatment_Lateral","PILO","CTL_SAL_IP")

# M_PILO_IP_A_CR
coldata_df_M_PILO_IP_A_CR <- coldata_df_M %>% filter(TYPE_2 == "M_PILO_IP_A_CR" | TYPE_2 == "CTL")
filtered_samples_M_PILO_IP_A_CR <- rownames(coldata_df_M_PILO_IP_A_CR)
adjusted_df_M_PILO_IP_A_CR <- adjusted_TPM_id_df_M[, colnames(adjusted_TPM_id_df_M) %in% filtered_samples_M_PILO_IP_A_CR]
adjusted_df_GSEA <- adjusted_df_M_PILO_IP_A_CR
coldata_df_GSEA <- coldata_df_M_PILO_IP_A_CR
Target_comparison <- 'M_PILO_IP_A_CR'
contrast_vector <- c("Treatment_Lateral","PILO","PILO_SAL_IP")

# M_KAI_IH_IPSI_A_HA
coldata_df_M_KAI_IH_IPSI_A_HA <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IH_IPSI_A_HA")
filtered_samples_M_KAI_IH_IPSI_A_HA <- rownames(coldata_df_M_KAI_IH_IPSI_A_HA)
adjusted_df_M_KAI_IH_IPSI_A_HA <- adjusted_TPM_id_df_M[, colnames(adjusted_TPM_id_df_M) %in% filtered_samples_M_KAI_IH_IPSI_A_HA]
adjusted_df_GSEA <- adjusted_df_M_KAI_IH_IPSI_A_HA
coldata_df_GSEA <- coldata_df_M_KAI_IH_IPSI_A_HA
Target_comparison <- 'M_KAI_IH_IPSI_A_HA'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_IPSI","CTL_SAL_HIPPO_IPSI")

# M_KAI_IH_IPSI_A_AC
coldata_df_M_KAI_IH_IPSI_A_AC <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IH_IPSI_A_AC")
filtered_samples_M_KAI_IH_IPSI_A_AC <- rownames(coldata_df_M_KAI_IH_IPSI_A_AC)
adjusted_df_M_KAI_IH_IPSI_A_AC <- adjusted_TPM_id_df_M[, colnames(adjusted_TPM_id_df_M) %in% filtered_samples_M_KAI_IH_IPSI_A_AC]
adjusted_df_GSEA <- adjusted_df_M_KAI_IH_IPSI_A_AC
coldata_df_GSEA <- coldata_df_M_KAI_IH_IPSI_A_AC
Target_comparison <- 'M_KAI_IH_IPSI_A_AC'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_IPSI","CTL_SAL_HIPPO_IPSI")

# M_KAI_IH_IPSI_A_IM
coldata_df_M_KAI_IH_IPSI_A_IM <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IH_IPSI_A_IM")
filtered_samples_M_KAI_IH_IPSI_A_IM <- rownames(coldata_df_M_KAI_IH_IPSI_A_IM)
adjusted_df_M_KAI_IH_IPSI_A_IM <- adjusted_TPM_id_df_M[, colnames(adjusted_TPM_id_df_M) %in% filtered_samples_M_KAI_IH_IPSI_A_IM]
adjusted_df_GSEA <- adjusted_df_M_KAI_IH_IPSI_A_IM
coldata_df_GSEA <- coldata_df_M_KAI_IH_IPSI_A_IM
Target_comparison <- 'M_KAI_IH_IPSI_A_IM'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_IPSI","CTL_SAL_HIPPO_IPSI")

# M_KAI_IH_IPSI_A_CR
coldata_df_M_KAI_IH_IPSI_A_CR <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IH_IPSI_A_CR")
filtered_samples_M_KAI_IH_IPSI_A_CR <- rownames(coldata_df_M_KAI_IH_IPSI_A_CR)
adjusted_df_M_KAI_IH_IPSI_A_CR <- adjusted_TPM_id_df_M[, colnames(adjusted_TPM_id_df_M) %in% filtered_samples_M_KAI_IH_IPSI_A_CR]
adjusted_df_GSEA <- adjusted_df_M_KAI_IH_IPSI_A_CR
coldata_df_GSEA <- coldata_df_M_KAI_IH_IPSI_A_CR
Target_comparison <- 'M_KAI_IH_IPSI_A_CR'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_IPSI","CTL_SAL_HIPPO_IPSI")

# M_KAI_IA_IPSI_A_AC
coldata_df_M_KAI_IA_IPSI_A_AC <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IA_IPSI_A_AC")
filtered_samples_M_KAI_IA_IPSI_A_AC <- rownames(coldata_df_M_KAI_IA_IPSI_A_AC)
adjusted_df_M_KAI_IA_IPSI_A_AC <- adjusted_TPM_id_df_M[, colnames(adjusted_TPM_id_df_M) %in% filtered_samples_M_KAI_IA_IPSI_A_AC]
adjusted_df_GSEA <- adjusted_df_M_KAI_IA_IPSI_A_AC
coldata_df_GSEA <- coldata_df_M_KAI_IA_IPSI_A_AC
Target_comparison <- 'M_KAI_IA_IPSI_A_AC'
contrast_vector <- c("Treatment_Lateral","KAI_AMG_IPSI","CTL_SAL_AMG_IPSI")

# M_KAI_IA_IPSI_A_IM
coldata_df_M_KAI_IA_IPSI_A_IM <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IA_IPSI_A_IM")
filtered_samples_M_KAI_IA_IPSI_A_IM <- rownames(coldata_df_M_KAI_IA_IPSI_A_IM)
adjusted_df_M_KAI_IA_IPSI_A_IM <- adjusted_TPM_id_df_M[, colnames(adjusted_TPM_id_df_M) %in% filtered_samples_M_KAI_IA_IPSI_A_IM]
adjusted_df_GSEA <- adjusted_df_M_KAI_IA_IPSI_A_IM
coldata_df_GSEA <- coldata_df_M_KAI_IA_IPSI_A_IM
Target_comparison <- 'M_KAI_IA_IPSI_A_IM'
contrast_vector <- c("Treatment_Lateral","KAI_AMG_IPSI","CTL_SAL_AMG_IPSI")

# M_KAI_IA_IPSI_A_CR
coldata_df_M_KAI_IA_IPSI_A_CR <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IA_IPSI_A_CR")
filtered_samples_M_KAI_IA_IPSI_A_CR <- rownames(coldata_df_M_KAI_IA_IPSI_A_CR)
adjusted_df_M_KAI_IA_IPSI_A_CR <- adjusted_TPM_id_df_M[, colnames(adjusted_TPM_id_df_M) %in% filtered_samples_M_KAI_IA_IPSI_A_CR]
adjusted_df_GSEA <- adjusted_df_M_KAI_IA_IPSI_A_CR
coldata_df_GSEA <- coldata_df_M_KAI_IA_IPSI_A_CR
Target_comparison <- 'M_KAI_IA_IPSI_A_CR'
contrast_vector <- c("Treatment_Lateral","KAI_AMG_IPSI","CTL_SAL_AMG_IPSI")

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
gsea_output_path <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/GSEA_02.Matrix_GSEA/gsea_formatted_matrix_",Target_comparison,".txt")
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
output_file_path <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/GSEA_02.Matrix_GSEA/phenotype_matrix_",Target_comparison,".txt")
writeLines(phenotype_matrix, con = output_file_path)

