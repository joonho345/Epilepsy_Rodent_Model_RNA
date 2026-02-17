library(readr)
library(tidyverse)
library(dplyr)

#### import TPM Count matrix with gene id ####
adjusted_TPM_id_df <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/08.Enrichment/GSEA_01.Matrix_GTF/adjusted_merged_matrix_1_id_TPM.txt", 
                                 sep = "\t", header = TRUE, row.names = 1)
#### import coldata ####
coldata_df <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/02.Quantification/filtered_coldata.txt",
                         sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE, 
                         fill = TRUE, quote = "", comment.char = "")


#### Exclude outliers ####
OUTLIER_samples <- c("SRR10867969", "SRR10868146", "SRR16522693SRR16522703", "SRR16522694SRR16522704", "SRR28007092", "SRR3628384", "SRR10868129", "SRR15406573", "SRR10868083", "SRR16522692SRR16522702")
coldata_df <- coldata_df[!rownames(coldata_df) %in% OUTLIER_samples, ]
adjusted_TPM_id_df <- adjusted_TPM_id_df[, !colnames(adjusted_TPM_id_df) %in% OUTLIER_samples]


#### Export only HIPPO ####
coldata_df_hippo <- coldata_df %>%
  filter(Diagnosis == "MTLEHS" | Diagnosis == "MTLE" | Diagnosis == "NL") %>%
  filter(Brain_Location == "HIPPO") %>%
  mutate(across(everything(), ~ replace(.x, .x == "", "N/A")))
filtered_samples_hippo <- rownames(coldata_df_hippo)
adjusted_TPM_id_df_hippo <- adjusted_TPM_id_df[, colnames(adjusted_TPM_id_df) %in% filtered_samples_hippo]


#### Grouping Diagnoses for GSEA ####
coldata_df_hippo_GSEA <- coldata_df_hippo %>%
  mutate(Diagnosis_GSEA = case_when(
    Diagnosis %in% c("MTLEHS", "MTLE") ~ "MTLEALL",
    TRUE ~ Diagnosis  # Keep other diagnoses unchanged
  ))

# Extract the filtered sample names for GSEA
filtered_samples_GSEA <- rownames(coldata_df_hippo_GSEA)
adjusted_df_hippo_GSEA <- adjusted_TPM_id_df_hippo[, colnames(adjusted_TPM_id_df_hippo) %in% filtered_samples_GSEA]
# Ensure that `coldata_df_hippo_GSEA` and `adjusted_df_GSEA` have matched sample names
adjusted_df_hippo_GSEA <- adjusted_df_hippo_GSEA[, filtered_samples_GSEA]


###################################################################
#### GSEA Matrix ####
gsea_expr_df <- adjusted_df_hippo_GSEA
gsea_coldata_df <- coldata_df_hippo_GSEA
group_name <- "GSEA"

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
gsea_output_path <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/08.Enrichment/GSEA_02.Matrix_GSEA/gsea_formatted_matrix_",group_name,"_hippo.txt")
writeLines(header_lines, con = gsea_output_path)
write.table(gsea_expr_df, file = gsea_output_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)


#### GSEA Matrix - Phenotype matrix ####
diagnosis_vector <- gsea_coldata_df$Diagnosis_GSEA
num_samples <- length(diagnosis_vector)
num_groups <- length(unique(diagnosis_vector))
first_line <- paste(num_samples, num_groups, 1, sep = " ")

group_names <- unique(diagnosis_vector)
second_line <- paste("#", paste(group_names, collapse = " "), sep = " ")
third_line <- paste(diagnosis_vector, collapse = "\t")

phenotype_matrix <- c(first_line, second_line, third_line)

# Write the phenotype matrix to a file
output_file_path <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/08.Enrichment/GSEA_02.Matrix_GSEA/phenotype_matrix_",group_name,"_hippo.txt")
writeLines(phenotype_matrix, con = output_file_path)

