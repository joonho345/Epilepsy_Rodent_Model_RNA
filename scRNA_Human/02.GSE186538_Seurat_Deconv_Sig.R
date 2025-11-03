# Load necessary libraries
library(Seurat)
library(Matrix)

# Set working directory
setwd('/home/joonho345/3_RNA/scRNA_Human/')
GSE <- 'GSE186538'


#### Get All Genes for human_hippo ####
# Get the sparse matrix (counts)
expression_matrix <- human_hippo@assays$RNA@layers$data
expression_df <- as.data.frame(expression_matrix, check.names = FALSE)
features <- rownames(human_hippo)
rownames(expression_df) <- features
expression_df$Gene <- rownames(expression_df)
expression_df <- expression_df[, c(ncol(expression_df), 1:(ncol(expression_df) - 1))]
head(expression_df[, 1:50])
length(rownames(expression_df))
# cluster
cluster_annotations <- as.character(human_hippo@meta.data$cluster)
header <- c("Gene", cluster_annotations)
output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/GSE186538_scRNA_matrix_no_colnames.txt"
write.table(expression_df, file = output_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
formatted_output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/GSE186538_scRNA_matrix_all_cluster_label.txt"
writeLines(paste(header, collapse = "\t"), con = formatted_output_path)
file.append(formatted_output_path, output_path)
# cluster_1
cluster_annotations <- as.character(human_hippo@meta.data$cluster_1)
header <- c("Gene", cluster_annotations)
formatted_output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/GSE186538_scRNA_matrix_all_cluster_1.txt"
writeLines(paste(header, collapse = "\t"), con = formatted_output_path)
file.append(formatted_output_path, output_path)
# cluster_2
cluster_annotations <- as.character(human_hippo@meta.data$cluster_2)
header <- c("Gene", cluster_annotations)
formatted_output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/GSE186538_scRNA_matrix_all_cluster_2.txt"
writeLines(paste(header, collapse = "\t"), con = formatted_output_path)
file.append(formatted_output_path, output_path)
file.remove(output_path)


#### Extract 2000 Top Variable Genes for human_hippo ####
human_hippo_1 <- FindVariableFeatures(human_hippo, nfeatures = 2000)
top_variable_genes_1 <- VariableFeatures(human_hippo_1)
expression_matrix <- human_hippo_1@assays$RNA@layers$data
expression_df_top_1 <- as.data.frame(expression_matrix, check.names = FALSE)
features <- rownames(human_hippo_1)
rownames(expression_df_top_1) <- features
expression_df_top_1$Gene <- rownames(expression_df_top_1)
expression_df_top_1 <- expression_df_top_1[, c(ncol(expression_df_top_1), 1:(ncol(expression_df_top_1) - 1))]
expression_df_top_1 <- expression_df_top_1[rownames(expression_df_top_1) %in% top_variable_genes_1, ]
head(expression_df_top_1[, 1:50])
length(rownames(expression_df_top_1))
# cluster
cluster_annotations <- as.character(human_hippo_1@meta.data$cluster)
header <- c("Gene", cluster_annotations)
output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/GSE186538_scRNA_matrix_no_colnames.txt"
write.table(expression_df_top_1, file = output_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
formatted_output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/GSE186538_scRNA_matrix_top2000_cluster_label.txt"
writeLines(paste(header, collapse = "\t"), con = formatted_output_path)
file.append(formatted_output_path, output_path)
# cluster_1
cluster_annotations <- as.character(human_hippo_1@meta.data$cluster_1)
header <- c("Gene", cluster_annotations)
formatted_output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/GSE186538_scRNA_matrix_top2000_cluster_1.txt"
writeLines(paste(header, collapse = "\t"), con = formatted_output_path)
file.append(formatted_output_path, output_path)
# cluster_2
cluster_annotations <- as.character(human_hippo_1@meta.data$cluster_2)
header <- c("Gene", cluster_annotations)
formatted_output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/GSE186538_scRNA_matrix_top2000_cluster_2.txt"
writeLines(paste(header, collapse = "\t"), con = formatted_output_path)
file.append(formatted_output_path, output_path)
file.remove(output_path)


#### Extract 1000 Top Variable Genes for human_hippo ####
human_hippo_2 <- FindVariableFeatures(human_hippo, nfeatures = 1000)
top_variable_genes_2 <- VariableFeatures(human_hippo_2)
expression_matrix <- human_hippo_2@assays$RNA@layers$data
expression_df_top_2 <- as.data.frame(expression_matrix, check.names = FALSE)
features <- rownames(human_hippo_2)
rownames(expression_df_top_2) <- features
expression_df_top_2$Gene <- rownames(expression_df_top_2)
expression_df_top_2 <- expression_df_top_2[, c(ncol(expression_df_top_2), 1:(ncol(expression_df_top_2) - 1))]
expression_df_top_2 <- expression_df_top_2[rownames(expression_df_top_2) %in% top_variable_genes_2, ]
head(expression_df_top_2[, 1:50])
length(rownames(expression_df_top_2))
# cluster
cluster_annotations <- as.character(human_hippo_2@meta.data$cluster)
header <- c("Gene", cluster_annotations)
output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/GSE186538_scRNA_matrix_no_colnames.txt"
write.table(expression_df_top_2, file = output_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
formatted_output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/GSE186538_scRNA_matrix_top1000_cluster_label.txt"
writeLines(paste(header, collapse = "\t"), con = formatted_output_path)
file.append(formatted_output_path, output_path)
# cluster_1
cluster_annotations <- as.character(human_hippo_2@meta.data$cluster_1)
header <- c("Gene", cluster_annotations)
formatted_output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/GSE186538_scRNA_matrix_top1000_cluster_1.txt"
writeLines(paste(header, collapse = "\t"), con = formatted_output_path)
file.append(formatted_output_path, output_path)
# cluster_2
cluster_annotations <- as.character(human_hippo_2@meta.data$cluster_2)
header <- c("Gene", cluster_annotations)
formatted_output_path <- "/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/GSE186538_scRNA_matrix_top1000_cluster_2.txt"
writeLines(paste(header, collapse = "\t"), con = formatted_output_path)
file.append(formatted_output_path, output_path)
file.remove(output_path)
