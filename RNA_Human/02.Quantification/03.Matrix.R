library(dplyr)


#### coldata ####
coldata <- read.csv("/data/project/1_Epilepsy_RNA/RNA_Human/Raw_Data/Patients_sheet_FTP_0.csv", header = TRUE)
rownames(coldata) <- coldata$Run
ncol(coldata)
nrow(coldata)
coldata_df <- as.data.frame(coldata, stringsAsFactors = FALSE)

#### matrix ####
# Define the path to your files and list all .htseq.count.txt files
path <- "/data/project/1_Epilepsy_RNA/RNA_Human/02.Quantification/01.HTseq"
files <- list.files(path, pattern = "*.htseq.count.txt", full.names = TRUE)
# Read each file and store it in a list of data frames
dfs <- lapply(files, function(file) {
  df <- read.table(file, sep = "\t", header = FALSE, col.names = c("Gene", gsub(".htseq.count.txt", "", basename(file))))
  return(df)
})
# Merge all data frames by "Gene"
merged_df <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), dfs)
merged_df <- merged_df[-(1:5), ]
# Remove genes whose names start with 'ENSG'
merged_df <- merged_df[!grepl("^ENSG", merged_df$Gene), ]
rownames(merged_df) <- merged_df$Gene
merged_df <- merged_df[, -1]


#### Order of samples ####
# coldata & expression matrix -> 661 Samples 
sorted_colnames <- sort(colnames(merged_df))
sorted_rownames <- sort(rownames(coldata_df))
merged_df <- merged_df[, sorted_colnames]
coldata_df <- coldata_df[sorted_rownames, ]
filtered_samples <- rownames(coldata_df)
merged_df <- merged_df[, colnames(merged_df) %in% filtered_samples]


#### Gene_name: Save the matrices and coldata ####
merged_matrix <- as.matrix(merged_df)
write.table(merged_matrix, file = "/home/joonho345/1_Epilepsy_RNA/RNA_Human/02.Quantification/merged_matrix.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
coldata_matrix <- as.matrix(coldata_df)
write.table(coldata_matrix, file = "/home/joonho345/1_Epilepsy_RNA/RNA_Human/02.Quantification/filtered_coldata.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
