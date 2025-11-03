# Load necessary libraries
library(Seurat)
library(Matrix)

# Set working directory
setwd('/home/joonho345/3_RNA/scRNA_Human/')
GSE <- '2021Trans'


#### Get All Genes for hippo_trans ####
# Get the sparse matrix (counts)
expression_matrix <- hippo_trans@assays$originalexp@counts
expression_df <- as.data.frame(expression_matrix, check.names = FALSE)
features <- rownames(hippo_trans)
rownames(expression_df) <- features
expression_df$Gene <- rownames(expression_df)
expression_df <- expression_df[, c(ncol(expression_df), 1:(ncol(expression_df) - 1))]
head(expression_df[, 1:50])
length(rownames(expression_df))
# cluster
cluster_annotations <- as.character(hippo_trans@meta.data$cellType)
header <- c("Gene", cluster_annotations)
output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.2021Trans_Seurat_Deconv/2021Trans_scRNA_matrix_no_colnames.txt"
write.table(expression_df, file = output_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
formatted_output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.2021Trans_Seurat_Deconv/2021Trans_scRNA_matrix_all_cluster_label.txt"
writeLines(paste(header, collapse = "\t"), con = formatted_output_path)
file.append(formatted_output_path, output_path)
# cluster_1
cluster_annotations <- as.character(hippo_trans@meta.data$cluster_1)
header <- c("Gene", cluster_annotations)
formatted_output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.2021Trans_Seurat_Deconv/2021Trans_scRNA_matrix_all_cluster_1.txt"
writeLines(paste(header, collapse = "\t"), con = formatted_output_path)
file.append(formatted_output_path, output_path)


#### Extract 2000 Top Variable Genes for hippo_trans ####
hippo_trans_1 <- FindVariableFeatures(hippo_trans, nfeatures = 2000)
top_variable_genes_1 <- VariableFeatures(hippo_trans_1)
expression_matrix <- hippo_trans_1@assays$originalexp@counts
expression_df_top_1 <- as.data.frame(expression_matrix, check.names = FALSE)
features <- rownames(hippo_trans_1)
rownames(expression_df_top_1) <- features
expression_df_top_1$Gene <- rownames(expression_df_top_1)
expression_df_top_1 <- expression_df_top_1[, c(ncol(expression_df_top_1), 1:(ncol(expression_df_top_1) - 1))]
expression_df_top_1 <- expression_df_top_1[rownames(expression_df_top_1) %in% top_variable_genes_1, ]
head(expression_df_top_1[, 1:50])
length(rownames(expression_df_top_1))
# cluster
cluster_annotations <- as.character(hippo_trans_1@meta.data$cellType)
header <- c("Gene", cluster_annotations)
output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.2021Trans_Seurat_Deconv/2021Trans_scRNA_matrix_no_colnames.txt"
write.table(expression_df_top_1, file = output_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
formatted_output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.2021Trans_Seurat_Deconv/2021Trans_scRNA_matrix_top2000_cluster_label.txt"
writeLines(paste(header, collapse = "\t"), con = formatted_output_path)
file.append(formatted_output_path, output_path)
# cluster_1
cluster_annotations <- as.character(hippo_trans_1@meta.data$cluster_1)
header <- c("Gene", cluster_annotations)
formatted_output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.2021Trans_Seurat_Deconv/2021Trans_scRNA_matrix_top2000_cluster_1.txt"
writeLines(paste(header, collapse = "\t"), con = formatted_output_path)
file.append(formatted_output_path, output_path)


#### Extract 1000 Top Variable Genes for hippo_trans ####
hippo_trans_2 <- FindVariableFeatures(hippo_trans, nfeatures = 1000)
top_variable_genes_2 <- VariableFeatures(hippo_trans_2)
expression_matrix <- hippo_trans_2@assays$originalexp@counts
expression_df_top_2 <- as.data.frame(expression_matrix, check.names = FALSE)
features <- rownames(hippo_trans_2)
rownames(expression_df_top_2) <- features
expression_df_top_2$Gene <- rownames(expression_df_top_2)
expression_df_top_2 <- expression_df_top_2[, c(ncol(expression_df_top_2), 1:(ncol(expression_df_top_2) - 1))]
expression_df_top_2 <- expression_df_top_2[rownames(expression_df_top_2) %in% top_variable_genes_2, ]
head(expression_df_top_2[, 1:50])
length(rownames(expression_df_top_2))
# cluster
cluster_annotations <- as.character(hippo_trans_2@meta.data$cluster)
header <- c("Gene", cluster_annotations)
output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.2021Trans_Seurat_Deconv/2021Trans_scRNA_matrix_no_colnames.txt"
write.table(expression_df_top_2, file = output_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
formatted_output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.2021Trans_Seurat_Deconv/2021Trans_scRNA_matrix_top1000_cluster_label.txt"
writeLines(paste(header, collapse = "\t"), con = formatted_output_path)
file.append(formatted_output_path, output_path)
# cluster_1
cluster_annotations <- as.character(hippo_trans_2@meta.data$cluster_1)
header <- c("Gene", cluster_annotations)
formatted_output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.2021Trans_Seurat_Deconv/2021Trans_scRNA_matrix_top1000_cluster_1.txt"
writeLines(paste(header, collapse = "\t"), con = formatted_output_path)
file.append(formatted_output_path, output_path)

