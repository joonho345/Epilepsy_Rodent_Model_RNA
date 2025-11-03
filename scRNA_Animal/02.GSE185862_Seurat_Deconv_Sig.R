# Load necessary libraries
library(Seurat)
library(Matrix)
library(Cairo)

# Set working directory
setwd('/home/joonho345/3_RNA/scRNA_Animal/')
GSE <- 'GSE185862'

#### Function to Write Sparse Matrix Once ####
write_sparse_matrix_once <- function(expression_matrix, genes, output_path) {
  file_conn <- file(output_path, open = "wt")
  for (i in seq_len(nrow(expression_matrix))) {
    gene_expression <- as.numeric(expression_matrix[i, ])
    gene_line <- paste(c(genes[i], gene_expression), collapse = "\t")
    writeLines(gene_line, con = file_conn)
  }
  close(file_conn)
}
#### Function to Append Cluster Annotations ####
write_cluster_annotations <- function(cluster_annotations, expression_matrix_path, output_cluster_path) {
  cluster_header <- paste(c("Gene", cluster_annotations), collapse = "\t")
  expression_data <- readLines(expression_matrix_path)
  expression_data <- c(cluster_header, expression_data)
  writeLines(expression_data, con = output_cluster_path)
}
#### Remove the intermediate expression matrix file ####
remove_intermediate_file <- function(file_path) {
  if (file.exists(file_path)) {
    file.remove(file_path)
  }
}

#### Get All Genes for hip_subset ####
# Get the sparse matrix (counts)
expression_matrix <- hip_subset@assays$RNA@counts
genes <- rownames(expression_matrix)
# Cluster label annotations
cluster_annotations <- as.character(hip_subset@meta.data$cluster_label)
cluster_annotations_all_1 <- as.character(hip_subset@meta.data$cluster_1)
cluster_annotations_all_2 <- as.character(hip_subset@meta.data$cluster_2)
# Write expression matrix once
output_matrix_path <- "/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/GSE185862_scRNA_matrix_expression_all.txt"
write_sparse_matrix_once(expression_matrix, genes, output_matrix_path)
# Create separate cluster files by appending annotations to the same expression matrix
write_cluster_annotations(cluster_annotations, output_matrix_path, "/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/GSE185862_scRNA_matrix_all_cluster_label.txt")
write_cluster_annotations(cluster_annotations_all_1, output_matrix_path, "/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/GSE185862_scRNA_matrix_all_cluster_1.txt")
write_cluster_annotations(cluster_annotations_all_2, output_matrix_path, "/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/GSE185862_scRNA_matrix_all_cluster_2.txt")
# Remove intermediate expression matrix file
remove_intermediate_file(output_matrix_path)


#### Extract 2000 Top Variable Genes for hip_subset_1 ####
hip_subset_1 <- FindVariableFeatures(hip_subset, nfeatures = 2000)
top_variable_genes_1 <- VariableFeatures(hip_subset_1)
# Filter the expression matrix to only the top 2000 variable genes
expression_matrix_1 <- expression_matrix[top_variable_genes_1, ]
# Cluster label annotations
cluster_annotations_1 <- as.character(hip_subset_1@meta.data$cluster_label)
cluster_annotations_1_1 <- as.character(hip_subset_1@meta.data$cluster_1)
cluster_annotations_1_2 <- as.character(hip_subset_1@meta.data$cluster_2)
# Write expression matrix once for top 2000 genes
output_matrix_path_1 <- "/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/GSE185862_scRNA_matrix_expression_top2000.txt"
write_sparse_matrix_once(expression_matrix_1, top_variable_genes_1, output_matrix_path_1)
# Create separate cluster files for top 2000 variable genes
write_cluster_annotations(cluster_annotations_1, output_matrix_path_1, "/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/GSE185862_scRNA_matrix_top2000_cluster_label.txt")
write_cluster_annotations(cluster_annotations_1_1, output_matrix_path_1, "/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/GSE185862_scRNA_matrix_top2000_cluster_1.txt")
write_cluster_annotations(cluster_annotations_1_2, output_matrix_path_1, "/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/GSE185862_scRNA_matrix_top2000_cluster_2.txt")
# Remove intermediate expression matrix file
remove_intermediate_file(output_matrix_path_1)


#### Extract 1000 Top Variable Genes for hip_subset_2 ####
hip_subset_2 <- FindVariableFeatures(hip_subset, nfeatures = 1000)
top_variable_genes_2 <- VariableFeatures(hip_subset_2)
# Filter the expression matrix to only the top 1000 variable genes
expression_matrix_2 <- expression_matrix[top_variable_genes_2, ]
# Cluster label annotations
cluster_annotations_2 <- as.character(hip_subset_2@meta.data$cluster_label)
cluster_annotations_2_1 <- as.character(hip_subset_2@meta.data$cluster_1)
cluster_annotations_2_2 <- as.character(hip_subset_2@meta.data$cluster_2)
# Write expression matrix once for top 1000 genes
output_matrix_path_2 <- "/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/GSE185862_scRNA_matrix_expression_top1000.txt"
write_sparse_matrix_once(expression_matrix_2, top_variable_genes_2, output_matrix_path_2)
# Create separate cluster files for top 1000 variable genes
write_cluster_annotations(cluster_annotations_2, output_matrix_path_2, "/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/GSE185862_scRNA_matrix_top1000_cluster_label.txt")
write_cluster_annotations(cluster_annotations_2_1, output_matrix_path_2, "/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/GSE185862_scRNA_matrix_top1000_cluster_1.txt")
write_cluster_annotations(cluster_annotations_2_2, output_matrix_path_2, "/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/GSE185862_scRNA_matrix_top1000_cluster_2.txt")
# Remove intermediate expression matrix file
remove_intermediate_file(output_matrix_path_2)
