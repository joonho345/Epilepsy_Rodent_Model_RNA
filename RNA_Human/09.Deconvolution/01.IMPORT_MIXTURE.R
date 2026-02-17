#!/usr/bin/env Rscript
# Import TPM matrix for CIBERSORTx
# Copy TPM matrix from 03.Normalization, add GeneSymbol header, and save to 01.IMPORT_MIXTURE

# Source and destination directories
source_dir <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/03.Normalization"
dest_dir <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/01.IMPORT_MIXTURE"

# File names
source_file <- file.path(source_dir, "adjusted_merged_matrix_1_TPM.txt")
dest_file <- file.path(dest_dir, "CIBERSORTx_adjusted_merged_matrix_1_TPM.txt")


# Read the TPM matrix
# First, check the first line to see if it has GeneSymbol header
first_line <- readLines(source_file, n = 1)
first_field <- strsplit(first_line, "\t")[[1]][1]

if (first_field == "GeneSymbol" || first_field == "NAME") {
  # File already has GeneSymbol or NAME header
  cat("File already has gene column header:", first_field, "\n")
  data <- read.table(source_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                     check.names = FALSE)
  if (colnames(data)[1] != "GeneSymbol") {
    # Rename NAME to GeneSymbol if needed
    colnames(data)[1] <- "GeneSymbol"
  }
  output_data <- data
} else {
  # The first column contains gene names (row names) without header
  cat("Detected row names in first column without header...\n")
  # Read with row.names = 1 to properly handle the structure
  data <- read.table(source_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                     row.names = 1, check.names = FALSE)
  
  # Convert to data frame with GeneSymbol as first column
  output_data <- data.frame(
    GeneSymbol = rownames(data),
    data,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(output_data) <- NULL
}

# Write the output file with GeneSymbol as first column
write.table(output_data, file = dest_file, sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = TRUE)
