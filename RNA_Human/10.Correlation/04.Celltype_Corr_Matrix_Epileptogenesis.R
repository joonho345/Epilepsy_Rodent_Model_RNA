# Load necessary libraries
library(dplyr)
library(pheatmap)
library(Cairo)
library(ComplexHeatmap)
library(circlize) # For colorRamp2


################################################
################ Cell types ###################
################################################
# Define cell types to be analyzed
celltypes <- c('ExN', 'InN', 'Astro', 'Micro', 'Oligo', 'OPC', 'Endo')
Human_group <- 'FILTERED_1_MTLEALL_NL'
Human_groups_name <- '_1'

# Choose Species
Species <- 'Mouse'
orthologs_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_mm_17128.txt"
coldata_df <- DESeq_coldata_df_M # import data from 00.Matrix_M.R

# Load orthologs
ortholog_one2one <- read.table(orthologs_path, sep = "\t", header = TRUE, fill = TRUE)
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one$Gene.name), ]
Species_Gene_name <- paste0(Species, '.gene.name')
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one[[Species_Gene_name]]), ]
human_genes <- ortholog_one2one$Gene.name
animal_genes <- ortholog_one2one[[Species_Gene_name]]

# Filter to only H_gene_set_W genes (epileptogenesis GO gene sets)
human_genes <- human_genes[human_genes %in% H_gene_set_W]
ortholog_one2one <- ortholog_one2one[ortholog_one2one$Gene.name %in% H_gene_set_W, ]
animal_genes <- ortholog_one2one[[Species_Gene_name]]

correlation_matrices <- list()
pvalue_matrices <- list()
for (celltype in celltypes) {
  l2fc_df <- get(paste0("DESeq_l2fc_", celltype, "_M")) # import mouse data from 00.Matrix_M.R
  human_l2fc_df <- get(paste0("DESeq_l2fc_", celltype, "_FILTERED")) # import human data
  
  correlation_matrix <- matrix(nrow = 1, ncol = ncol(l2fc_df))  
  rownames(correlation_matrix) <- celltype 
  colnames(correlation_matrix) <- colnames(l2fc_df)
  pvalue_matrix <- matrix(nrow = 1, ncol = ncol(l2fc_df))  
  rownames(pvalue_matrix) <- celltype 
  colnames(pvalue_matrix) <- colnames(l2fc_df)
  
  human_data_filtered <- human_l2fc_df[rownames(human_l2fc_df) %in% human_genes, ]
  animal_data_filtered <- l2fc_df[rownames(l2fc_df) %in% animal_genes, ]
  
  rownames(animal_data_filtered) <- human_genes[match(rownames(animal_data_filtered), animal_genes)]
  common_genes <- intersect(rownames(human_data_filtered), rownames(animal_data_filtered))
  human_data_filtered <- human_data_filtered[common_genes, ]
  animal_data_filtered <- animal_data_filtered[common_genes, ]
  
  human_l2fc <- human_data_filtered[[Human_group]]  # Single human group column
  animal_l2fc <- animal_data_filtered
  
  correlation_results <- sapply(animal_l2fc, function(x) cor(human_l2fc, x, use = "complete.obs"))
  pvalue_results <- sapply(animal_l2fc, function(x) cor.test(human_l2fc, x, use = "complete.obs")$p.value)

  correlation_matrix[1, ] <- correlation_results
  correlation_matrices[[celltype]] <- correlation_matrix
  pvalue_matrix[1, ] <- pvalue_results
  pvalue_matrices[[celltype]] <- pvalue_matrix
}

Mouse_order <- c(
  "M_GSE185862_M_KAI_IH_IPSI_A_HA", "M_GSE185862_M_KAI_IH_IPSI_A_AC", "M_GSE185862_M_KAI_IH_IPSI_A_IM", "M_GSE185862_M_KAI_IH_IPSI_A_CR",
  "M_GSE185862_M_KAI_IH_CON_A_AC", "M_GSE185862_M_KAI_IH_CON_A_IM", "M_GSE185862_M_KAI_IH_CON_A_CR", 
  "M_GSE185862_M_KAI_IA_IPSI_A_AC", "M_GSE185862_M_KAI_IA_IPSI_A_IM", "M_GSE185862_M_KAI_IA_IPSI_A_CR", "M_GSE185862_M_KAI_IP_A_HA", 
  "M_GSE185862_M_PILO_IP_A_HA", "M_GSE185862_M_PILO_IP_A_AC", "M_GSE185862_M_PILO_IP_A_IM", "M_GSE185862_M_PILO_IP_A_CR"
)

correlation_matrices <- lapply(correlation_matrices, function(mat) mat[, Mouse_order, drop = FALSE])
merged_correlation_matrix <- do.call(rbind, correlation_matrices)
pvalue_matrices <- lapply(pvalue_matrices, function(mat) mat[, Mouse_order, drop = FALSE])
merged_pvalue_matrix <- do.call(rbind, pvalue_matrices)

# Remove M_GSE185862_ prefix and convert underscores to hyphens for column names
colnames(merged_correlation_matrix) <- gsub("^M_GSE185862_", "", colnames(merged_correlation_matrix))
colnames(merged_correlation_matrix) <- gsub("_", "-", colnames(merged_correlation_matrix))
rownames(merged_correlation_matrix) <- gsub("_", "-", rownames(merged_correlation_matrix))

colnames(merged_pvalue_matrix) <- gsub("^M_GSE185862_", "", colnames(merged_pvalue_matrix))
colnames(merged_pvalue_matrix) <- gsub("_", "-", colnames(merged_pvalue_matrix))
rownames(merged_pvalue_matrix) <- gsub("_", "-", rownames(merged_pvalue_matrix))

output_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/04.CellType_Matrix/01.Correlation_Epileptogenesis_Allcelltype", 
                      Human_groups_name, ".txt")
write.table(merged_correlation_matrix, file = output_file, sep = "\t", quote = FALSE, col.names = TRUE)

output_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/04.CellType_Matrix/02.Pvalue_Epileptogenesis_Allcelltype", 
                      Human_groups_name, ".txt")
write.table(merged_pvalue_matrix, file = output_file, sep = "\t", quote = FALSE, col.names = TRUE)




