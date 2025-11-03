library(devtools)
library(sva)

#### matrix and coldata ####
merged_matrix
merged_matrix_id
coldata_df

#### vectors ####
vector_PRJNA <- coldata_df$PRJNA
vector_Diagnosis <- coldata_df$Diagnosis
vector_Brain_Location <- coldata_df$Brain_Location

numeric_Brain_Location <- as.numeric(factor(vector_Brain_Location))
numeric_Diagnosis <- as.numeric(factor(vector_Diagnosis))
covar_mat <- cbind(numeric_Diagnosis, numeric_Brain_Location)


# gene name
adjusted_merged_matrix <- ComBat_seq(merged_matrix, batch=vector_PRJNA)
adjusted_merged_matrix_1 <- ComBat_seq(merged_matrix, batch=vector_PRJNA, group=vector_Diagnosis)
adjusted_merged_matrix_2 <- ComBat_seq(merged_matrix, batch=vector_PRJNA, group=NULL, covar_mod=covar_mat)

write.table(adjusted_merged_matrix, file = "/home/joonho345/3_RNA/RNA_Human/02.Quantification/adjusted_merged_matrix.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(adjusted_merged_matrix_1, file = "/home/joonho345/3_RNA/RNA_Human/02.Quantification/adjusted_merged_matrix_1.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(adjusted_merged_matrix_2, file = "/home/joonho345/3_RNA/RNA_Human/02.Quantification/adjusted_merged_matrix_2.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


# gene id
adjusted_merged_matrix_id <- ComBat_seq(merged_matrix_id, batch=vector_PRJNA)
adjusted_merged_matrix_id_1 <- ComBat_seq(merged_matrix_id, batch=vector_PRJNA, group=vector_Diagnosis)
adjusted_merged_matrix_id_2 <- ComBat_seq(merged_matrix_id, batch=vector_PRJNA, group=NULL, covar_mod=covar_mat)

write.table(adjusted_merged_matrix_id, file = "/home/joonho345/3_RNA/RNA_Human/02.Quantification/adjusted_merged_matrix_id.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(adjusted_merged_matrix_id_1, file = "/home/joonho345/3_RNA/RNA_Human/02.Quantification/adjusted_merged_matrix_id_1.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(adjusted_merged_matrix_id_2, file = "/home/joonho345/3_RNA/RNA_Human/02.Quantification/adjusted_merged_matrix_id_2.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
