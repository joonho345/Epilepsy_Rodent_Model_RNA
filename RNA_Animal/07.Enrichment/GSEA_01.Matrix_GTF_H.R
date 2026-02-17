library(dplyr)
library(tidyr)

##################### MOUSE ##########################
#### Create GTF Matrix ####
gtf_path <- "/data/resource/reference/mouse/Mus_musculus.GRCm39.112.gtf"
gtf_table <- read.table(gtf_path, header = FALSE, skip = 5, sep = "\t", quote = "", comment.char = "")

gtf_table <- gtf_table %>% dplyr::select(-c(V1, V4, V5, V6, V7, V8))
gtf_table <- gtf_table[gtf_table$V3 == "gene", ]

# Create gene_id column
extract_gene_id <- function(v9_string) {
  # Use a regular expression to find and extract the gene_id value
  matches <- regmatches(v9_string, regexpr('gene_id "\\S+"', v9_string))
  gene_id <- ifelse(length(matches) > 0, sub('gene_id "(\\S+)"', '\\1', matches), NA)
  return(gene_id)
}
gene_ids <- sapply(gtf_table$V9, extract_gene_id)
gtf_table$gene_id <- gene_ids

# Create gene_name column
extract_gene_name <- function(v9_string) {
  # Use a regular expression to find and extract the gene_id value
  matches <- regmatches(v9_string, regexpr('gene_name "\\S+"', v9_string))
  gene_name <- ifelse(length(matches) > 0, sub('gene_name "(\\S+)"', '\\1', matches), NA)
  return(gene_name)
}
gene_names <- sapply(gtf_table$V9, extract_gene_name)
gtf_table$gene_name <- gene_names

# Remove rows where 'gene_name' is NA
gtf_table <- gtf_table %>% filter(!is.na(gene_name))


#### Expression matrix ####
adjusted_TPM_df <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/03.Normalization/adjusted_merged_matrix_TPM_M.txt", 
                              sep = "\t", header = TRUE, row.names = 1)
adjusted_TPM_id_df <- adjusted_TPM_df

# Change the gene name into gene id
adjusted_TPM_id_df$gene_name <- rownames(adjusted_TPM_df)
adjusted_TPM_id_df <- merge(adjusted_TPM_id_df, gtf_table[, c("gene_id", "gene_name")], by = "gene_name", all.x = TRUE)
adjusted_TPM_id_df <- adjusted_TPM_id_df %>% dplyr::select(-gene_name)
adjusted_TPM_id_df <- aggregate(. ~ gene_id, data = adjusted_TPM_id_df, FUN = sum)
rownames(adjusted_TPM_id_df) <- adjusted_TPM_id_df$gene_id
adjusted_TPM_id_df <- adjusted_TPM_id_df %>% dplyr::select(-gene_id)

# SAVE expression matrix with gene id
adjusted_TPM_id_matrix <- as.matrix(adjusted_TPM_id_df)
write.table(adjusted_TPM_id_matrix, file = "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/GSEA_01.Matrix_GTF/adjusted_merged_matrix_id_TPM_M.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


##################### RAT ##########################
#### Create GTF Matrix ####
gtf_path <- "/data/resource/reference/mouse/Rattus_norvegicus.mRatBN7.2.112.gtf"
gtf_table <- read.table(gtf_path, header = FALSE, skip = 5, sep = "\t", quote = "", comment.char = "")

gtf_table <- gtf_table %>% dplyr::select(-c(V1, V4, V5, V6, V7, V8))
gtf_table <- gtf_table[gtf_table$V3 == "gene", ]

# Create gene_id column
extract_gene_id <- function(v9_string) {
  # Use a regular expression to find and extract the gene_id value
  matches <- regmatches(v9_string, regexpr('gene_id "\\S+"', v9_string))
  gene_id <- ifelse(length(matches) > 0, sub('gene_id "(\\S+)"', '\\1', matches), NA)
  return(gene_id)
}
gene_ids <- sapply(gtf_table$V9, extract_gene_id)
gtf_table$gene_id <- gene_ids

# Create gene_name column
extract_gene_name <- function(v9_string) {
  # Use a regular expression to find and extract the gene_id value
  matches <- regmatches(v9_string, regexpr('gene_name "\\S+"', v9_string))
  gene_name <- ifelse(length(matches) > 0, sub('gene_name "(\\S+)"', '\\1', matches), NA)
  return(gene_name)
}
gene_names <- sapply(gtf_table$V9, extract_gene_name)
gtf_table$gene_name <- gene_names

# Remove rows where 'gene_name' is NA
gtf_table <- gtf_table %>% filter(!is.na(gene_name))


#### Expression matrix ####
adjusted_TPM_df <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/03.Normalization/adjusted_merged_matrix_TPM_R.txt", 
                              sep = "\t", header = TRUE, row.names = 1)
adjusted_TPM_id_df <- adjusted_TPM_df

# Change the gene name into gene id
adjusted_TPM_id_df$gene_name <- rownames(adjusted_TPM_df)
adjusted_TPM_id_df <- merge(adjusted_TPM_id_df, gtf_table[, c("gene_id", "gene_name")], by = "gene_name", all.x = TRUE)
adjusted_TPM_id_df <- adjusted_TPM_id_df %>% dplyr::select(-gene_name)
adjusted_TPM_id_df <- aggregate(. ~ gene_id, data = adjusted_TPM_id_df, FUN = sum)
rownames(adjusted_TPM_id_df) <- adjusted_TPM_id_df$gene_id
adjusted_TPM_id_df <- adjusted_TPM_id_df %>% dplyr::select(-gene_id)

# SAVE expression matrix with gene id
adjusted_TPM_id_matrix <- as.matrix(adjusted_TPM_id_df)
write.table(adjusted_TPM_id_matrix, file = "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/GSEA_01.Matrix_GTF/adjusted_merged_matrix_id_TPM_R.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

