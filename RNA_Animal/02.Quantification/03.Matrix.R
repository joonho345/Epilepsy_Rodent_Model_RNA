library(dplyr)


#### coldata ####
coldata_all <- read.csv("/home/joonho345/3_RNA/RNA_Animal/Raw_Data/Animal_sheet_FTP.csv", header = TRUE) # total 349
coldata_M <- coldata_all %>% filter(Species == "Musmusculus") # M 277
coldata_R <- coldata_all %>% filter(Species == "Rattusnorvegicus") # R 72
rownames(coldata_M) <- coldata_M$Run
rownames(coldata_R) <- coldata_R$Run
coldata_df_M <- as.data.frame(coldata_M, stringsAsFactors = FALSE)
coldata_df_R <- as.data.frame(coldata_R, stringsAsFactors = FALSE)

coldata_matrix_M <- as.matrix(coldata_df_M)
write.table(coldata_matrix_M, file = "/home/joonho345/3_RNA/RNA_Animal/02.Quantification/filtered_coldata_M.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
coldata_matrix_R <- as.matrix(coldata_df_R)
write.table(coldata_matrix_R, file = "/home/joonho345/3_RNA/RNA_Animal/02.Quantification/filtered_coldata_R.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

#### DESeq_coldata_1 ####
DESeq_coldata_all <- read.csv("/home/joonho345/3_RNA/RNA_Animal/Raw_Data/Animal_sheet_FTP_1.csv", header = TRUE) # total 349
DESeq_coldata_M <- DESeq_coldata_all %>% filter(Species == "Musmusculus")
DESeq_coldata_R <- DESeq_coldata_all %>% filter(Species == "Rattusnorvegicus")
rownames(DESeq_coldata_M) <- DESeq_coldata_M$Group
rownames(DESeq_coldata_R) <- DESeq_coldata_R$Group
DESeq_coldata_df_M <- as.data.frame(DESeq_coldata_M, stringsAsFactors = FALSE)
DESeq_coldata_df_R <- as.data.frame(DESeq_coldata_R, stringsAsFactors = FALSE)

DESeq_coldata_matrix_M <- as.matrix(DESeq_coldata_df_M)
write.table(DESeq_coldata_matrix_M, file = "/home/joonho345/3_RNA/RNA_Animal/03.DESeq/DESeq_coldata_1_M.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
DESeq_coldata_matrix_R <- as.matrix(DESeq_coldata_df_R)
write.table(DESeq_coldata_matrix_R, file = "/home/joonho345/3_RNA/RNA_Animal/03.DESeq/DESeq_coldata_1_R.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

#### DESeq_coldata_2 ####
DESeq_coldata_all <- read.csv("/home/joonho345/3_RNA/RNA_Animal/Raw_Data/Animal_sheet_FTP_2.csv", header = TRUE) # total 349
DESeq_coldata_M <- DESeq_coldata_all %>% filter(Species == "M")
DESeq_coldata_R <- DESeq_coldata_all %>% filter(Species == "R")
rownames(DESeq_coldata_M) <- DESeq_coldata_M$Group
rownames(DESeq_coldata_R) <- DESeq_coldata_R$Group
DESeq_coldata_df_M <- as.data.frame(DESeq_coldata_M, stringsAsFactors = FALSE)
DESeq_coldata_df_R <- as.data.frame(DESeq_coldata_R, stringsAsFactors = FALSE)

DESeq_coldata_matrix_M <- as.matrix(DESeq_coldata_df_M)
write.table(DESeq_coldata_matrix_M, file = "/home/joonho345/3_RNA/RNA_Animal/03.DESeq_2/DESeq_coldata_2_M.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
DESeq_coldata_matrix_R <- as.matrix(DESeq_coldata_df_R)
write.table(DESeq_coldata_matrix_R, file = "/home/joonho345/3_RNA/RNA_Animal/03.DESeq_2/DESeq_coldata_2_R.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


######################################################################
#### matrix -> Total 349 Samples####
# Define the path to your files and list all .htseq.count.txt files
path <- "/data/project/HS/RNA_Animal/02.Quantification/01.HTseq/"
files <- list.files(path, pattern = "*.htseq.count.txt", full.names = TRUE)
files_M <- files[basename(files) %in% paste0(rownames(coldata_df_M), ".htseq.count.txt")]
files_R <- files[basename(files) %in% paste0(rownames(coldata_df_R), ".htseq.count.txt")]
# Read each file and store it in a list of data frames
read_and_format <- function(file) {
  df <- read.table(file, sep = "\t", header = FALSE, col.names = c("Gene", gsub(".htseq.count.txt", "", basename(file))))
  return(df)
}

# Read and merge data frames for Musmusculus
dfs_M <- lapply(files_M, read_and_format)
merged_df_M <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), dfs_M)
merged_df_M <- merged_df_M[-(1:5), ]  # Remove the first 5 rows (assumed to be non-gene data)
rownames(merged_df_M) <- merged_df_M$Gene
merged_df_M <- merged_df_M[, -1]

# Read and merge data frames for Rattusnorvegicus
dfs_R <- lapply(files_R, read_and_format)
merged_df_R <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), dfs_R)
merged_df_R <- merged_df_R[-(1:5), ]  # Remove the first 5 rows (assumed to be non-gene data)
rownames(merged_df_R) <- merged_df_R$Gene
merged_df_R <- merged_df_R[, -1]


#### Ensure the order of sample IDs ####
# Musmusculus: Reorder the columns of the expression matrix to match the rownames of coldata
merged_df_M <- merged_df_M[, rownames(coldata_df_M)]
# Rattusnorvegicus: Reorder the columns of the expression matrix to match the rownames of coldata
merged_df_R <- merged_df_R[, rownames(coldata_df_R)]


#### Gene_id: Save the matrices and coldata ####
merged_matrix_M <- as.matrix(merged_df_M)
write.table(merged_matrix_M, file = "/home/joonho345/3_RNA/RNA_Animal/02.Quantification/merged_matrix_id_M.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
merged_matrix_R <- as.matrix(merged_df_R)
write.table(merged_matrix_R, file = "/home/joonho345/3_RNA/RNA_Animal/02.Quantification/merged_matrix_id_R.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


#### Change the gene id into gene name ####
# Mouse
merged_df_M$gene_ids_M <- rownames(merged_df_M)
merged_df_M <- merge(merged_df_M, gtf_table_M[, c("gene_ids_M", "gene_names_M")], by = "gene_ids_M", all.x = TRUE)
merged_df_M <- merged_df_M %>% dplyr::select(-gene_ids_M)
merged_df_M <- aggregate(. ~ gene_names_M, data = merged_df_M, FUN = sum)
rownames(merged_df_M) <- merged_df_M$gene_names_M
merged_df_M <- merged_df_M %>% dplyr::select(-gene_names_M)
# Rat
merged_df_R$gene_ids_R <- rownames(merged_df_R)
merged_df_R <- merge(merged_df_R, gtf_table_R[, c("gene_ids_R", "gene_names_R")], by = "gene_ids_R", all.x = TRUE)
merged_df_R <- merged_df_R %>% dplyr::select(-gene_ids_R)
merged_df_R <- aggregate(. ~ gene_names_R, data = merged_df_R, FUN = sum)
rownames(merged_df_R) <- merged_df_R$gene_names_R
merged_df_R <- merged_df_R %>% dplyr::select(-gene_names_R)


#### Gene_name: Save the matrices and coldata ####
# Musmusculus
merged_matrix_M <- as.matrix(merged_df_M)
write.table(merged_matrix_M, file = "/home/joonho345/3_RNA/RNA_Animal/02.Quantification/merged_matrix_M.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
# Rattusnorvegicus
merged_matrix_R <- as.matrix(merged_df_R)
write.table(merged_matrix_R, file = "/home/joonho345/3_RNA/RNA_Animal/02.Quantification/merged_matrix_R.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

