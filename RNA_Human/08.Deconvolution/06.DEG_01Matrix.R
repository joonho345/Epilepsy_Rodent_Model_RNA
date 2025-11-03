library(dplyr)


#### import coldata_df from 00.Matrix_H.R
coldata_df

#### import HiRES matrix ####
## cluster_1 ##
celltype <- 'ExN'
celltype <- 'InN'
celltype <- 'Astro'
celltype <- 'Micro'
celltype <- 'Oligo'
# 2019CA all cluster_1
#Dataset_H <- '2019CA'
#type <- 'all'
Dataset_H <- 'GSE186538'
type <- 'top1000'
clustertype <- 'cluster_1'
file_path <- paste0("/home/joonho345/3_RNA/RNA_Human/08.Deconvolution/03.CIBERSORT_", Dataset_H,
                    "_04HiRes/", Dataset_H, "_scRNA_matrix_", type, "_", clustertype, "/Merged_", celltype, "_1.txt")

## cluster_2 ##
celltype <- 'ExN_DG'
celltype <- 'ExN_CA'
celltype <- 'ExN_SUB'
# GSE186538 top1000 cluster_2
Dataset_H <- 'GSE186538'
type <- 'top1000'
clustertype <- 'cluster_2'
file_path <- paste0("/home/joonho345/3_RNA/RNA_Human/08.Deconvolution/03.CIBERSORT_", Dataset_H,
                      "_04HiRes/", Dataset_H, "_scRNA_matrix_", type, "_", clustertype, "/Merged_", celltype, "_1.txt")

# print version
version <- paste0(celltype, Dataset_H, type, clustertype)
print(version)
Expr_matrix <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


#### choose groups to compare #### 
filtered_samples_1 <- rownames(coldata_df_1)
Expr_matrix_1 <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_1]

filtered_samples_2 <- rownames(coldata_df_2)
Expr_matrix_2 <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_2]

filtered_samples_3 <- rownames(coldata_df_3)
Expr_matrix_3 <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_3]




