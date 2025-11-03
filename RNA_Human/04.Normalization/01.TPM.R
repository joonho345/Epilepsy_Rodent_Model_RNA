library(GeoTcgaData)
library(dplyr)

# https://github.com/YuLab-SMU/GeoTcgaData
load("/home/joonho345/resources/GeoTcgaData/gene_cov.rda")
# v1.992

#### import matrix ####
merged_matrix <- as.matrix(read.table("/home/joonho345/3_RNA/RNA_Human/02.Quantification/merged_matrix.txt", sep = "\t", header = TRUE, row.names = 1))
assign("merged_matrix_TPM", countToTpm(merged_matrix, keyType = "SYMBOL", gene_cov = gene_cov))
write.table(merged_matrix_TPM, file = "/home/joonho345/3_RNA/RNA_Human/04.Normalization/merged_matrix_TPM.txt", sep = "\t", quote = FALSE)

adjusted_merged_matrix <- as.matrix(read.table("/home/joonho345/3_RNA/RNA_Human/02.Quantification/adjusted_merged_matrix.txt", sep = "\t", header = TRUE, row.names = 1))
assign("adjusted_merged_matrix_TPM", countToTpm(adjusted_merged_matrix, keyType = "SYMBOL", gene_cov = gene_cov))
write.table(adjusted_merged_matrix_TPM, file = "/home/joonho345/3_RNA/RNA_Human/04.Normalization/adjusted_merged_matrix_TPM.txt", sep = "\t", quote = FALSE)

adjusted_merged_matrix_1 <- as.matrix(read.table("/home/joonho345/3_RNA/RNA_Human/02.Quantification/adjusted_merged_matrix_1.txt", sep = "\t", header = TRUE, row.names = 1))
assign("adjusted_merged_matrix_1_TPM", countToTpm(adjusted_merged_matrix_1, keyType = "SYMBOL", gene_cov = gene_cov))
write.table(adjusted_merged_matrix_1_TPM, file = "/home/joonho345/3_RNA/RNA_Human/04.Normalization/adjusted_merged_matrix_1_TPM.txt", sep = "\t", quote = FALSE)

adjusted_merged_matrix_2 <- as.matrix(read.table("/home/joonho345/3_RNA/RNA_Human/02.Quantification/adjusted_merged_matrix_2.txt", sep = "\t", header = TRUE, row.names = 1))
assign("adjusted_merged_matrix_2_TPM", countToTpm(adjusted_merged_matrix_2, keyType = "SYMBOL", gene_cov = gene_cov))
write.table(adjusted_merged_matrix_2_TPM, file = "/home/joonho345/3_RNA/RNA_Human/04.Normalization/adjusted_merged_matrix_2_TPM.txt", sep = "\t", quote = FALSE)

# CountToTpm
# assign(result_Name, countToTpm(Target_Matrix, keyType = "SYMBOL", gene_cov = gene_cov))
