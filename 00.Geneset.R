library(dplyr)

## ortholog ##
Target_set <- Target_genes
Target_set <- Target_genes_recurrent

ortholog_file <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_mm_17128.txt"
ortholog_file <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_rn_16211.txt"
ortholog_file <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_mm_rn_15711.txt"
orthologs <- read.table(ortholog_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gene_set_ortholog <- Target_set
# mouse
mapped_genes <- orthologs %>%
  filter(Gene.name %in% gene_set_ortholog)
hs_mm_gene_set <- mapped_genes %>%
  pull(`Mouse.gene.name`)
print(hs_mm_gene_set)
missing_genes <- setdiff(gene_set_ortholog, mapped_genes$Gene.name)
print(missing_genes)
# rat
mapped_genes <- orthologs %>%
  filter(Gene.name %in% gene_set_ortholog)
hs_rn_gene_set <- mapped_genes %>%
  pull(`Rat.gene.name`)
print(hs_rn_gene_set)
missing_genes <- setdiff(gene_set_ortholog, mapped_genes$Gene.name)
print(missing_genes)


############## BIOMART ############## 
#### for One Set ####
file_directory <- "/home/joonho345/resources/GeneSet/BioMart_GO/"
file_list <- list.files(file_directory, pattern = "\\.txt$", full.names = TRUE)
human_gene_set_names <- c()
mouse_gene_set_names <- c()
rat_gene_set_names <- c()

for (file_path in file_list) {
  data <- read.table(file_path, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
  # Get the file name without directory and extension
  file_name <- tools::file_path_sans_ext(basename(file_path))
  # Human gene set (unique 'Gene name' without duplicated entries)
  H_gene_set <- unique(data$Gene.name)
  H_gene_set <- H_gene_set[H_gene_set != ""]  # Remove empty entries
  # Mouse gene set where 'Mouse homology type' is 'ortholog_one2one'
  M_gene_set <- data %>%
    filter(Mouse.homology.type == "ortholog_one2one" & Mouse.gene.name != "") %>%
    dplyr::select(Mouse.gene.name) %>%
    distinct() %>%
    pull()
  # Rat gene set where 'Rat homology type' is 'ortholog_one2one'
  R_gene_set <- data %>%
    filter(Rat.homology.type == "ortholog_one2one" & Rat.gene.name != "") %>%
    dplyr::select(Rat.gene.name) %>%
    distinct() %>%
    pull()
  # Assign the vectors to the global environment
  assign(paste0("H_gene_set_", file_name), H_gene_set, envir = .GlobalEnv)
  assign(paste0("M_gene_set_", file_name), M_gene_set, envir = .GlobalEnv)
  assign(paste0("R_gene_set_", file_name), R_gene_set, envir = .GlobalEnv)
  # Add the names of the gene set vectors to the respective lists of names
  human_gene_set_names <- c(human_gene_set_names, paste0("H_gene_set_", file_name))
  mouse_gene_set_names <- c(mouse_gene_set_names, paste0("M_gene_set_", file_name))
  rat_gene_set_names <- c(rat_gene_set_names, paste0("R_gene_set_", file_name))
}
print(human_gene_set_names)
print(mouse_gene_set_names)
print(rat_gene_set_names)


#### for Upper Set ####
file_directory <- "/home/joonho345/resources/GeneSet/BioMart_GO/"
file_list <- list.files(file_directory, pattern = "\\.txt$", full.names = TRUE)
upper_categories <- c("IC", "NT", "ND", "NI", "TNF", "MF", "NG", "GG", "SP")
human_gene_set_names_upper <- c()
mouse_gene_set_names_upper <- c()
rat_gene_set_names_upper <- c()

for (upper_category in upper_categories) {
  human_combined_genes <- c()
  mouse_combined_genes <- c()
  rat_combined_genes <- c()
  for (file_path in file_list) {
    file_name <- tools::file_path_sans_ext(basename(file_path))
    # Check if the file name contains the upper category
    if (grepl(upper_category, file_name)) {
      data <- read.table(file_path, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
      # Collect unique human genes
      H_gene_set <- unique(data$Gene.name)
      H_gene_set <- H_gene_set[H_gene_set != ""]  # Remove empty entries
      human_combined_genes <- unique(c(human_combined_genes, H_gene_set))
      # Collect unique mouse genes
      M_gene_set <- data %>%
        filter(Mouse.homology.type == "ortholog_one2one" & Mouse.gene.name != "") %>%
        dplyr::select(Mouse.gene.name) %>%
        distinct() %>%
        pull()
      mouse_combined_genes <- unique(c(mouse_combined_genes, M_gene_set))
      # Collect unique rat genes
      R_gene_set <- data %>%
        filter(Rat.homology.type == "ortholog_one2one" & Rat.gene.name != "") %>%
        dplyr::select(Rat.gene.name) %>%
        distinct() %>%
        pull()
      rat_combined_genes <- unique(c(rat_combined_genes, R_gene_set))
    }
  }
  # Assign the combined gene sets to the global environment
  assign(paste0("H_gene_set_", upper_category), human_combined_genes, envir = .GlobalEnv)
  assign(paste0("M_gene_set_", upper_category), mouse_combined_genes, envir = .GlobalEnv)
  assign(paste0("R_gene_set_", upper_category), rat_combined_genes, envir = .GlobalEnv)
  # Add the upper category gene set names to the name lists
  human_gene_set_names_upper <- c(human_gene_set_names_upper, paste0("H_gene_set_", upper_category))
  mouse_gene_set_names_upper <- c(mouse_gene_set_names_upper, paste0("M_gene_set_", upper_category))
  rat_gene_set_names_upper <- c(rat_gene_set_names_upper, paste0("R_gene_set_", upper_category))
}
print(human_gene_set_names_upper)  # Names of all human upper category gene set vectors
print(mouse_gene_set_names_upper)  # Names of all mouse upper category gene set vectors
print(rat_gene_set_names_upper)    # Names of all rat upper category gene set vectors

## gene_set_W -> Create union of 9 gene sets (excluding ERK) for human, mouse, and rat
H_gene_set_W <- unique(unlist(mget(human_gene_set_names_upper, envir = .GlobalEnv)))
M_gene_set_W <- unique(unlist(mget(mouse_gene_set_names_upper, envir = .GlobalEnv)))
R_gene_set_W <- unique(unlist(mget(rat_gene_set_names_upper, envir = .GlobalEnv)))
assign("H_gene_set_W", H_gene_set_W, envir = .GlobalEnv)
assign("M_gene_set_W", M_gene_set_W, envir = .GlobalEnv)
assign("R_gene_set_W", R_gene_set_W, envir = .GlobalEnv)
cat("H_gene_set_W:", length(H_gene_set_W), "genes\n")
cat("M_gene_set_W:", length(M_gene_set_W), "genes\n")
cat("R_gene_set_W:", length(R_gene_set_W), "genes\n")

gene_sets <- list(
  NT = H_gene_set_NT, IC = H_gene_set_IC, NI = H_gene_set_NI, ND = H_gene_set_ND,
  TNF = H_gene_set_TNF, NG = H_gene_set_NG, GG = H_gene_set_GG, MF = H_gene_set_MF, SP = H_gene_set_SP
)
formatted_gene_sets <- sapply(names(gene_sets), function(set_name) {
  gene_list <- gene_sets[[set_name]]
  gene_count <- length(gene_list)  # Count number of genes
  paste(set_name, "(N =", gene_count, "):", paste(gene_list, collapse = ", "))
})
output_string <- paste(formatted_gene_sets, collapse = "\n")
output_file <- "/home/joonho345/resources/GeneSet/H_gene_set_ALL.txt"
writeLines(output_string, output_file)
cat(output_string, sep = "\n")


#### Subgroup ####
file_directory <- "/home/joonho345/resources/GeneSet/BioMart_GO/"
file_list <- list.files(file_directory, pattern = "\\.txt$", full.names = TRUE)
Neurogenesis <- c()
Gliogenesis <- c()
Synaptic_plasticity <- c()
Neurotransmission <- c()
Ion_channels <- c()
Neuroinflammation <- c()
Neuronal_death <- c()
Mossy_fiber_sprouting <- c()
TNF_signaling <- c()
# 1. With _number (keeping the entire GO term including the number)
extract_go_term <- function(file_path) {
  file_name <- basename(file_path)               # Extract the file name
  go_term <- tools::file_path_sans_ext(file_name) # Remove the file extension
  return(go_term)}
go_terms <- sapply(file_list, extract_go_term)
# 2. Without _number (keeping only the first two parts of the GO term)
for (go_term in go_terms) {
  if (grepl("NG", go_term)) {
    Neurogenesis <- c(Neurogenesis, go_term)
  } else if (grepl("GG", go_term)) {
    Gliogenesis <- c(Gliogenesis, go_term)
  } else if (grepl("SP", go_term)) {
    Synaptic_plasticity <- c(Synaptic_plasticity, go_term)
  } else if (grepl("NT", go_term)) {
    Neurotransmission <- c(Neurotransmission, go_term)
  } else if (grepl("IC", go_term)) {
    Ion_channels <- c(Ion_channels, go_term)
  } else if (grepl("NI", go_term)) {
    Neuroinflammation <- c(Neuroinflammation, go_term)
  } else if (grepl("ND", go_term)) {
    Neuronal_death <- c(Neuronal_death, go_term)
  } else if (grepl("MF", go_term)) {
    Mossy_fiber_sprouting <- c(Mossy_fiber_sprouting, go_term)
  } else if (grepl("TNF", go_term)) 
}
