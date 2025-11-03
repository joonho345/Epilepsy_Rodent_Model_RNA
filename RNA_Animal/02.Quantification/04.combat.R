library(devtools)
library(sva)

#### matrix and coldata ####
merged_matrix_M
merged_matrix_R
coldata_df_M
coldata_df_R

#### vectors ####
## Mouse
vector_Sequencing_M <- coldata_df_M$Sequencing
vector_Treatment_Lateral_M <- coldata_df_M$Treatment_Lateral

# gene name
adjusted_merged_matrix_M <- ComBat_seq(merged_matrix_M, batch=vector_Sequencing_M)
adjusted_merged_matrix_M_1 <- ComBat_seq(merged_matrix_M, batch=vector_Sequencing_M, group=vector_Treatment_Lateral_M)

write.table(adjusted_merged_matrix_M, file = "/home/joonho345/3_RNA/RNA_Animal/02.Quantification/adjusted_merged_matrix_M.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(adjusted_merged_matrix_M_1, file = "/home/joonho345/3_RNA/RNA_Animal/02.Quantification/adjusted_merged_matrix_M_1.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

## Rat
vector_Sequencing_R <- coldata_df_R$Sequencing
vector_Treatment_Lateral_R <- coldata_df_R$Treatment_Lateral

# gene name
adjusted_merged_matrix_R <- ComBat_seq(merged_matrix_R, batch=vector_Sequencing_R)
adjusted_merged_matrix_R_1 <- ComBat_seq(merged_matrix_R, batch=vector_Sequencing_R, group=vector_Treatment_Lateral_R)

write.table(adjusted_merged_matrix_R, file = "/home/joonho345/3_RNA/RNA_Animal/02.Quantification/adjusted_merged_matrix_R.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(adjusted_merged_matrix_R_1, file = "/home/joonho345/3_RNA/RNA_Animal/02.Quantification/adjusted_merged_matrix_R_1.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

