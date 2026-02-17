# Load necessary libraries
library(dplyr)
library(pheatmap)
library(Cairo)
library(ComplexHeatmap)
library(circlize) # For colorRamp2


################################################
################ GO terms ###################
################################################
# Choose Species
Species <- 'Mouse'
Species_abbr <- 'M'
orthologs_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_mm_17128.txt"
l2fc_df <- DESeq_l2fc_df_M # import data from 00.Matrix_M.R
coldata_df <- DESeq_coldata_df_M # import data from 00.Matrix_M.R
treatment_vector <- c("KAI-ST", "KAI-IP", "PILO-IP")
treatment_color_vector <- c("KAI-ST" = "#fc8d62", "KAI-IP" = "#8da0cb", "PILO-IP" = "#66c2a5")

# Choose Species
Species <- 'Rat'
Species_abbr <- 'R'
orthologs_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_rn_16211.txt"
l2fc_df <- DESeq_l2fc_df_R # import data from 00.Matrix_M.R
coldata_df <- DESeq_coldata_df_R # import data from 00.Matrix_M.R
treatment_vector <- c("KAI-IP", "KAI-SUB", "PILO-IP", "PPS", "AMG", "TBI")
treatment_color_vector <- c("KAI-IP" = "#8da0cb", "KAI-SUB" = "#ffd92f", "PILO-IP" = "#66c2a5", 
                            "TBI" = "#e78ac3", "AMG" = "#a6d854", "PPS" = "#e5c494")



# List of human groups to analyze
human_l2fc_df <- DESeq_l2fc_df_FILTERED
Human_groups <- c('FILTERED_1_MTLEALL_NL', 'FILTERED_2_MTLEHS_NL', 'FILTERED_3_MTLE_NL')
Human_groups_name <- '_FILTERED'

# Define upper categories
upper_categories <- c("NT", "IC", "NI", "ND", "TNF", "NG", "GG", "MF", "SP")


############## Load orthologs ############## 
ortholog_one2one <- read.table(orthologs_path, sep = "\t", header = TRUE, fill = TRUE)

# Remove duplicates in human and species gene names
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one$Gene.name), ]
Species_Gene_name <- paste0(Species, '.gene.name')
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one[[Species_Gene_name]]), ]

# Get the human and species gene names
human_genes <- ortholog_one2one$Gene.name
animal_genes <- ortholog_one2one[[Species_Gene_name]]


############## Whole Gene Set (No GO Term Filtering) ############## 
# GO data directory
GO_directory <- "/home/joonho345/resources/GeneSet/BioMart_GO/"
GO_list <- list.files(GO_directory, pattern = "\\.txt$", full.names = TRUE)

# Initialize correlation and p-value matrices
correlation_matrix <- matrix(nrow = length(upper_categories), ncol = ncol(l2fc_df))
pvalue_matrix <- matrix(nrow = length(upper_categories), ncol = ncol(l2fc_df))
rownames(correlation_matrix) <- upper_categories
rownames(pvalue_matrix) <- upper_categories
colnames(correlation_matrix) <- colnames(l2fc_df)
colnames(pvalue_matrix) <- colnames(l2fc_df)

# Loop through upper categories to compute correlations
for (upper_category in upper_categories) {
  human_gene_set <- get(paste0("H_gene_set_", upper_category))
  animal_gene_set <- get(paste0(Species_abbr, "_gene_set_", upper_category))
  
  current_human_genes <- human_genes[human_genes %in% human_gene_set]
  current_animal_genes <- animal_genes[human_genes %in% current_human_genes]
  
  human_data_filtered <- human_l2fc_df[rownames(human_l2fc_df) %in% current_human_genes, ]
  animal_data_filtered <- l2fc_df[rownames(l2fc_df) %in% current_animal_genes, ]
  rownames(animal_data_filtered) <- human_genes[match(rownames(animal_data_filtered), animal_genes)]
  
  common_genes <- intersect(rownames(human_data_filtered), rownames(animal_data_filtered))
  human_data_filtered <- human_data_filtered[common_genes, ]
  animal_data_filtered <- animal_data_filtered[common_genes, ]
  human_data_filtered <- as.matrix(human_data_filtered)
  animal_data_filtered <- as.matrix(animal_data_filtered)
  
  correlation_results <- numeric(ncol(animal_data_filtered))
  pvalue_results <- numeric(ncol(animal_data_filtered))
  
  for (i in seq_len(ncol(animal_data_filtered))) {
    correlation_results[i] <- cor(human_data_filtered[, 1], animal_data_filtered[, i], use = "complete.obs")
    pvalue_results[i] <- cor.test(human_data_filtered[, 1], animal_data_filtered[, i], use = "complete.obs")$p.value
  }
  
  correlation_matrix[upper_category, ] <- correlation_results
  pvalue_matrix[upper_category, ] <- pvalue_results
}


treatments <- setNames(coldata_df$Treatment, coldata_df$Group)
phases <- setNames(coldata_df$Phase, coldata_df$Group)
ages <- setNames(coldata_df$Age, coldata_df$Group)

colnames(correlation_matrix) <- gsub("_", "-", colnames(correlation_matrix))
colnames(pvalue_matrix) <- gsub("_", "-", colnames(pvalue_matrix))

rownames(correlation_matrix) <- gsub("^NT$", "Neurotransmission", rownames(correlation_matrix))
rownames(correlation_matrix) <- gsub("^IC$", "Ion Channel", rownames(correlation_matrix))
rownames(correlation_matrix) <- gsub("^NI$", "Neuroinflammation", rownames(correlation_matrix))
rownames(correlation_matrix) <- gsub("^ND$", "Neuronal Death", rownames(correlation_matrix))
rownames(correlation_matrix) <- gsub("^TNF$", "TNF Signaling", rownames(correlation_matrix))
rownames(correlation_matrix) <- gsub("^NG$", "Neurogenesis", rownames(correlation_matrix))
rownames(correlation_matrix) <- gsub("^GG$", "Gliogenesis", rownames(correlation_matrix))
rownames(correlation_matrix) <- gsub("^MF$", "Mossy Fiber Sprouting", rownames(correlation_matrix))
rownames(correlation_matrix) <- gsub("^SP$", "Synaptic Plasticity", rownames(correlation_matrix))
rownames(pvalue_matrix) <- gsub("^NT$", "Neurotransmission", rownames(pvalue_matrix))
rownames(pvalue_matrix) <- gsub("^IC$", "Ion Channel", rownames(pvalue_matrix))
rownames(pvalue_matrix) <- gsub("^NI$", "Neuroinflammation", rownames(pvalue_matrix))
rownames(pvalue_matrix) <- gsub("^ND$", "Neuronal Death", rownames(pvalue_matrix))
rownames(pvalue_matrix) <- gsub("^TNF$", "TNF Signaling", rownames(pvalue_matrix))
rownames(pvalue_matrix) <- gsub("^NG$", "Neurogenesis", rownames(pvalue_matrix))
rownames(pvalue_matrix) <- gsub("^GG$", "Gliogenesis", rownames(pvalue_matrix))
rownames(pvalue_matrix) <- gsub("^MF$", "Mossy Fiber Sprouting", rownames(pvalue_matrix))
rownames(pvalue_matrix) <- gsub("^SP$", "Synaptic Plasticity", rownames(pvalue_matrix))

output_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/05.GO_Matrix/01.Correlation_WholeGO_upper_", Species, "_FILTERED.txt")
write.table(correlation_matrix, file = output_file, sep = "\t", quote = FALSE, col.names = TRUE)

output_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/05.GO_Matrix/02.Pvalue_WholeGO_upper_", Species, "_FILTERED.txt")
write.table(pvalue_matrix, file = output_file, sep = "\t", quote = FALSE, col.names = TRUE)



