library(GeoTcgaData)
library(dplyr)

# https://github.com/YuLab-SMU/GeoTcgaData
load("/home/joonho345/resources/GeoTcgaData/gene_cov.rda")
# v1.992

#### import matrix ####
merged_matrix <- as.matrix(read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/02.Quantification/merged_matrix.txt", sep = "\t", header = TRUE, row.names = 1))
assign("merged_matrix_TPM", countToTpm(merged_matrix, keyType = "SYMBOL", gene_cov = gene_cov))
write.table(merged_matrix_TPM, file = paste0(output_dir, "/merged_matrix_TPM.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

adjusted_merged_matrix_1 <- as.matrix(read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/02.Quantification/adjusted_merged_matrix_1.txt", sep = "\t", header = TRUE, row.names = 1))
assign("adjusted_merged_matrix_1_TPM", countToTpm(adjusted_merged_matrix_1, keyType = "SYMBOL", gene_cov = gene_cov))
write.table(adjusted_merged_matrix_1_TPM, file = paste0(output_dir, "/adjusted_merged_matrix_1_TPM.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

