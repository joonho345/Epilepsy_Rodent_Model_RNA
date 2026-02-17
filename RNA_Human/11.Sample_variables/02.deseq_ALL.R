# Load necessary libraries
library(dplyr)

# Define input directory and variable list
directory_path <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/11.Sample_variables/02.deseq"
variable_list <- c("HS", "withoutHS", "IP", "NO", "IP_Sei", "IP_Other")

all_merged_dfs <- list()  # List to store matrices for final merging

for (var in variable_list) {
  # Get the list of files matching the current variable
  cell_files <- list.files(path = directory_path, pattern = paste0("DESeq_", var, ".txt$"), full.names = TRUE)
  if (length(cell_files) == 0) {
    cat("No files found for", var, "- Skipping...\n")
    next  # Skip if no files found
  }
  print(cell_files)
  list_dfs <- list()
  for (file_path in cell_files) {
    label <- tools::file_path_sans_ext(basename(file_path))  # Extract filename without extension
    df <- read.table(file_path, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)  # Read file
    if (!"log2FoldChange" %in% colnames(df)) {
      cat("Skipping", file_path, "- 'log2FoldChange' column not found.\n")
      next  # Skip file if it does not contain the expected column
    }
    df_filtered <- df[, "log2FoldChange", drop = FALSE]  # Extract 'log2FoldChange'
    colnames(df_filtered) <- label  # Rename column with file label
    list_dfs[[label]] <- df_filtered  # Store in list
  }
  
  if (length(list_dfs) == 0) {
    cat("No valid data frames for", var, "- Skipping merge step.\n")
    next
  }
  common_rows <- Reduce(intersect, lapply(list_dfs, rownames))
  list_dfs <- lapply(list_dfs, function(df) df[common_rows, , drop = FALSE])
  df_DESeq_ALL <- do.call(cbind, list_dfs)
  all_merged_dfs[[var]] <- df_DESeq_ALL  # Store merged dataset
}

if (length(all_merged_dfs) > 0) {
  common_genes <- Reduce(intersect, lapply(all_merged_dfs, rownames))
  all_merged_dfs <- lapply(all_merged_dfs, function(df) df[common_genes, , drop = FALSE])
  Final_Matrix_DESeq_ALL <- do.call(cbind, all_merged_dfs)
  colnames(Final_Matrix_DESeq_ALL) <- gsub("DESeq_", "", colnames(Final_Matrix_DESeq_ALL))
  final_output_file <- paste0(directory_path, "/DESeq_ALL.txt")
  write.table(Final_Matrix_DESeq_ALL, file = final_output_file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  cat("Final merged matrix has been created and saved as", final_output_file, "\n")
} else {
  cat("No valid matrices found for merging.\n")
}
