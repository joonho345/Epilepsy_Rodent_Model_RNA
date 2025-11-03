# Set the output directory
output_directory <- "/home/joonho345/3_RNA/RNA_Human/06.Deconvolution/02.IMPORT_GO/"

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
upper_categories <- c("NT", "ERK", "NG", "MF", "SP", "TLR", "TNF", "IC", "ND", "NI", "CC", "NR")
for (upper_category in upper_categories) {
  combined_genes <- c()
  
  # Loop through all files and collect genes from files containing the upper category in the name
  for (file_path in file_list) {
    file_name <- tools::file_path_sans_ext(basename(file_path))
    
    # Check if the file name contains the upper category
    if (grepl(upper_category, file_name)) {
      data <- read.table(file_path, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
      
      # Collect unique gene names from the current file
      gene_set <- unique(data$Gene.name)
      gene_set <- gene_set[gene_set != ""]  # Remove empty entries
      
      # Combine with the existing gene set for the upper category
      combined_genes <- unique(c(combined_genes, gene_set))
    }
  }
  
  upper_category_file <- paste0(output_directory, "H_gene_set_", upper_category, ".txt")
  writeLines(combined_genes, upper_category_file)
}
