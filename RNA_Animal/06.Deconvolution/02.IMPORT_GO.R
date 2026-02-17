# Set the output directory
output_directory <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/02.IMPORT_GO/"

# Create destination directory if it doesn't exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
}


file_directory <- "/home/joonho345/resources/GeneSet/BioMart_GO/"
file_list <- list.files(file_directory, pattern = "\\.txt$", full.names = TRUE)


####### Create individual GO terms txt files#######
for (file_path in file_list) {
  data <- read.table(file_path, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
  # Get the file name without directory and extension
  file_name <- tools::file_path_sans_ext(basename(file_path))

  # mouse gene set where 'Mouse homology type' is 'ortholog_one2one'
  M_gene_set <- data %>%
    filter(Mouse.homology.type == "ortholog_one2one" & Mouse.gene.name != "") %>%
    dplyr::select(Mouse.gene.name) %>%
    distinct() %>%
    pull()
  
  # Create file paths for each species
  mouse_gene_file <- paste0(output_directory, "M_gene_set_", file_name, ".txt")
  
  # Write each gene set to a file with newline delimiter
  writeLines(M_gene_set, mouse_gene_file)
}


####### Create Upper categories txt files#######
upper_categories <- c("NT", "NG", "MF", "SP", "TNF", "IC", "ND", "NI", "GG")
for (upper_category in upper_categories) {
  combined_genes <- c()
  
  # Loop through all files and collect genes from files containing the upper category in the name
  for (file_path in file_list) {
    file_name <- tools::file_path_sans_ext(basename(file_path))
    
    # Check if the file name contains the upper category
    if (grepl(upper_category, file_name)) {
      data <- read.table(file_path, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
      
      # Mouse gene set where 'Mouse homology type' is 'ortholog_one2one'
      gene_set <- data %>%
        filter(Mouse.homology.type == "ortholog_one2one" & Mouse.gene.name != "") %>%
        dplyr::select(Mouse.gene.name) %>%
        distinct() %>%
        pull()
      
      # Combine with the existing gene set for the upper category
      combined_genes <- unique(c(combined_genes, gene_set))
    }
  }
  
  # Create file paths for each upper category
  upper_category_file <- paste0(output_directory, "M_gene_set_", upper_category, ".txt")
  writeLines(combined_genes, upper_category_file)
}
