# Load necessary libraries
library(dplyr)

directory_path <- "/home/joonho345/3_RNA/RNA_Human/08.Deconvolution/06.DEG_02DESeq"
cohort_2019CA <- c("Astro", "ExN", "InN", "Oligo", "Micro")
cohort_GSE186538 <- c("Astro", "ExN", "InN", "Oligo", "Micro", "ExN_CA", "ExN_DG", "ExN_SUB")

cohorts <- list(
  "GSE186538" = cohort_GSE186538
)

# Loop through each cohort and cell type
for (cohort in names(cohorts)) {
  cell_types <- cohorts[[cohort]]
  
  for (cell_type in cell_types) {
    cell_files <- list.files(path = directory_path, pattern = paste0("^", cell_type, "_", cohort, ".*\\.txt$"), full.names = TRUE)
    print(cell_files)
    list_dfs <- list()

    # Loop through each file for the current cell type
    for (file_path in cell_files) {
      label <- tools::file_path_sans_ext(basename(file_path))  # Extract file name without extension
      df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)  # Read the file
      df_filtered <- df[, "log2FoldChange", drop = FALSE]  # Filter for 'log2FoldChange' column
      colnames(df_filtered) <- label  # Rename the column with the label
      list_dfs[[label]] <- df_filtered  # Add dataframe to the list
    }
    
    common_rows <- Reduce(intersect, lapply(list_dfs, rownames))
    list_dfs <- lapply(list_dfs, function(df) df[common_rows, , drop = FALSE])
    df_DESeq_ALL <- do.call(cbind, list_dfs)
    colnames(df_DESeq_ALL) <- gsub(paste0("^", cell_type, "_", cohort, "_"), "", colnames(df_DESeq_ALL))
    Matrix_DESeq_ALL <- as.matrix(df_DESeq_ALL)
    
    # Write the matrix to a file
    output_file <- paste0(directory_path, "/DESeq_", cell_type, "_ALL.txt")
    write.table(Matrix_DESeq_ALL, file = output_file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    
    # Print message indicating completion for the current cell type and cohort
    cat("Matrix for", cell_type, "in cohort", cohort, "has been created and saved as", output_file, "\n")
  }
}
