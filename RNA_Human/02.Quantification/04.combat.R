library(devtools)
library(sva)

#### matrix and coldata ####
# Read data from text files saved by 03.Matrix.R
merged_matrix <- as.matrix(read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/02.Quantification/merged_matrix.txt", 
                                       sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
coldata_df <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/02.Quantification/filtered_coldata.txt", 
                         sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE, 
                         fill = TRUE, quote = "", comment.char = "")

#### vectors ####
# Extract vectors in the same order as matrix columns
vector_PRJNA <- coldata_df$PRJNA
vector_Diagnosis <- coldata_df$Diagnosis
vector_Brain_Location <- coldata_df$Brain_Location

numeric_Brain_Location <- as.numeric(factor(vector_Brain_Location))
numeric_Diagnosis <- as.numeric(factor(vector_Diagnosis))
covar_mat <- cbind(numeric_Diagnosis, numeric_Brain_Location)


# gene name
adjusted_merged_matrix_1 <- ComBat_seq(merged_matrix, batch=vector_PRJNA, group=vector_Diagnosis)
write.table(adjusted_merged_matrix_1, file = "/home/joonho345/1_Epilepsy_RNA/RNA_Human/02.Quantification/adjusted_merged_matrix_1.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

