# Set the output directory
output_directory <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/02.IMPORT_GO/"

file_directory <- "/home/joonho345/resources/GeneSet/BioMart_GO/"
file_list <- list.files(file_directory, pattern = "\\.txt$", full.names = TRUE)


####### Create individual GO terms txt files####### 
for (file_path in file_list) {
  data <- read.table(file_path, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  # human gene set (unique 'Gene name' without duplicated entries)
  H_gene_set <- unique(data$Gene.name)
  H_gene_set <- H_gene_set[H_gene_set != ""]  # Remove empty entries
  
  # Create file paths for each species
  human_gene_file <- paste0(output_directory, "H_gene_set_", file_name, ".txt")
  
  # Write each gene set to a file with newline delimiter
  writeLines(H_gene_set, human_gene_file)
}


####### Create Upper categories txt files####### 
upper_categories <- c("IC", "NT", "ND", "NI", "TNF", "MF", "NG", "GG", "SP")
for (upper_category in upper_categories) {
  human_combined_genes <- character()  # Initialize as character vector
  
  # Loop through all files and collect genes from files containing the upper category in the name
  for (file_path in file_list) {
    file_name <- tools::file_path_sans_ext(basename(file_path))
    
    # Check if the file name contains the upper category
    if (grepl(upper_category, file_name)) {
      data <- read.table(file_path, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
      
      # Collect unique human genes
      H_gene_set <- unique(data$Gene.name)
      H_gene_set <- H_gene_set[H_gene_set != ""]  # Remove empty entries
      
      # Combine with the existing gene set for the upper category
      human_combined_genes <- unique(c(human_combined_genes, H_gene_set))
    }
  }
  
  # Only write if there are genes to write
  if (length(human_combined_genes) > 0) {
  upper_category_file <- paste0(output_directory, "H_gene_set_", upper_category, ".txt")
    writeLines(human_combined_genes, upper_category_file)
    cat("H_gene_set_", upper_category, ": ", length(human_combined_genes), " genes\n", sep = "")
  } else {
    cat("Warning: No genes found for category", upper_category, "\n")
  }
}
