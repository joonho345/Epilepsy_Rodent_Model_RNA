library(dplyr)


#### import coldata_df from 00.Matrix_H.R
coldata_df

#### import HiRES matrix ####
## CellType1 (Integration version) ##
celltype <- 'ExN_H'
celltype <- 'InN_H'
celltype <- 'Astro_H'
celltype <- 'Micro_H'
celltype <- 'Oligo_H'
celltype <- 'OPC_H'
celltype <- 'Endo_H'

Dataset_H <- 'GSE160189_186538_Integration'
type <- 'DEGs'
clustertype <- 'CellType1'
file_path <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/03.CIBERSORT_", Dataset_H,
                    "_04HiRes/", Dataset_H, "_scRNA_matrix_", type, "_", clustertype, "/Merged_", celltype, ".txt")

# print version
version <- paste0(celltype, "_", Dataset_H, "_", type, "_", clustertype)
print(version)
Expr_matrix <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


#### choose groups to compare #### 
filtered_samples_1 <- rownames(coldata_df_1)
Expr_matrix_1 <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_1]
