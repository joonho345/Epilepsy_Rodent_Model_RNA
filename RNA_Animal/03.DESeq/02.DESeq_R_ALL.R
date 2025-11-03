# Import 
csv_file <- "/home/joonho345/3_RNA/RNA_Animal/Raw_Data/Animal_sheet_FTP_2.csv"
file_data <- read.csv(csv_file, stringsAsFactors = FALSE)
file_data <- file_data %>% filter(Species == "R")
files_and_labels <- setNames(file_data$Filepath, file_data$Group)

# Initialize an empty list to store dataframes
list_dfs <- list()

# Loop through each file path and label
for (label in names(files_and_labels)) {
  file_path <- files_and_labels[[label]]
  df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df_filtered <- df[, "log2FoldChange", drop = FALSE]
  colnames(df_filtered) <- label
  list_dfs[[label]] <- df_filtered
}

# Matrix_DESeq_ALL
common_rows <- Reduce(intersect, lapply(list_dfs, rownames))
list_dfs <- lapply(list_dfs, function(df) df[common_rows, , drop = FALSE])
df_DESeq_ALL <- do.call(cbind, list_dfs)
Matrix_DESeq_ALL <- as.matrix(df_DESeq_ALL)

write.table(Matrix_DESeq_ALL, file = "/home/joonho345/3_RNA/RNA_Animal/03.DESeq_2/DESeq_R_ALL.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

