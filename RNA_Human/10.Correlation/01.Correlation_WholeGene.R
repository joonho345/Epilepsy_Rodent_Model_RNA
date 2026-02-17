# Load necessary libraries
library(dplyr)
library(pheatmap)
library(Cairo)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(grid)


# Choose Species
Species <- 'Mouse'
orthologs_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_mm_17128.txt"
l2fc_df <- DESeq_l2fc_df_M # import data from 00.Matrix_M.R
coldata_df <- DESeq_coldata_df_M # import data from 00.Matrix_M.R
treatment_vector <- c("KAI_IH", "KAI_IA", "KAI_IP", "PILO_IP")
treatment_color_vector <- c("KAI_IH" = "#fc8d62", "KAI_IA" = "#999999", "KAI_IP" = "#e5c494", "PILO_IP" = "#66c2a5",
                            "TBI" = "#e78ac3", "AMG" = "#a6d854", "KAI_SUB" = "#ffd92f", "PPS" = "#8da0cb")

Species <- 'Rat'
orthologs_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_rn_16211.txt"
l2fc_df <- DESeq_l2fc_df_R # import data from 00.Matrix_M.R
coldata_df <- DESeq_coldata_df_R # import data from 00.Matrix_M.R
treatment_vector <- c("KAI_IP", "KAI_SUB", "PILO_IP", "PPS", "AMG", "TBI")
treatment_color_vector <- c("KAI_IH" = "#fc8d62", "KAI_IA" = "#999999", "KAI_IP" = "#e5c494", "PILO_IP" = "#66c2a5",
                            "TBI" = "#e78ac3", "AMG" = "#a6d854", "KAI_SUB" = "#ffd92f", "PPS" = "#8da0cb")

human_l2fc_df <- DESeq_l2fc_df_FILTERED
Human_groups <- c('FILTERED_1_MTLEALL_NL', 'FILTERED_2_MTLEHS_NL', 'FILTERED_3_MTLE_NL')
Human_groups_name <- '_FILTERED'

############## Load orthologs ############## 
ortholog_one2one <- read.table(orthologs_path, sep = "\t", header = TRUE, fill = TRUE)

# Remove duplicates in human and species gene names
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one$Gene.name), ]
Species_Gene_name <- paste0(Species, '.gene.name')
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one[[Species_Gene_name]]), ]

# Get the human and species gene names
human_genes <- ortholog_one2one$Gene.name
animal_genes <- ortholog_one2one[[Species_Gene_name]]

############## Robust Statistics Function ##############
# Function for bootstrapping and permutation testing
get_robust_stats <- function(human_vec, animal_vec, n_perm = 1000) {
  # Remove NA pairs
  complete_idx <- complete.cases(human_vec, animal_vec)
  human_vec <- human_vec[complete_idx]
  animal_vec <- animal_vec[complete_idx]
  
  # Check if we have enough data points
  if (length(human_vec) < 3) {
    return(list(r = NA, p_perm = NA, ci_low = NA, ci_high = NA))
  }
  
  # 1. Observed Correlation
  obs_cor <- cor(human_vec, animal_vec, use = "complete.obs")
  
  # 2. Permutation p-value
  # Shuffle one vector n times and see how often random cor >= observed cor
  perm_cors <- replicate(n_perm, {
    cor(sample(human_vec), animal_vec, use = "complete.obs")
  })
  perm_p <- sum(abs(perm_cors) >= abs(obs_cor)) / n_perm
  
  # 3. Bootstrap CI (95%)
  boot_cors <- replicate(n_perm, {
    idx <- sample(1:length(human_vec), replace = TRUE)
    cor(human_vec[idx], animal_vec[idx], use = "complete.obs")
  })
  boot_ci <- quantile(boot_cors, probs = c(0.025, 0.975), na.rm = TRUE)
  
  return(list(r = obs_cor, p_perm = perm_p, ci_low = boot_ci[1], ci_high = boot_ci[2]))
}

### correlation_matrix ###
# Initialize matrices to store correlation results
correlation_matrix <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))
pvalue_matrix <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))
perm_pvalue_matrix <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))
bootstrap_ci_low_matrix <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))
bootstrap_ci_high_matrix <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))

rownames(correlation_matrix) <- Human_groups
colnames(correlation_matrix) <- colnames(l2fc_df)
rownames(pvalue_matrix) <- Human_groups
colnames(pvalue_matrix) <- colnames(l2fc_df)
rownames(perm_pvalue_matrix) <- Human_groups
colnames(perm_pvalue_matrix) <- colnames(l2fc_df)
rownames(bootstrap_ci_low_matrix) <- Human_groups
colnames(bootstrap_ci_low_matrix) <- colnames(l2fc_df)
rownames(bootstrap_ci_high_matrix) <- Human_groups
colnames(bootstrap_ci_high_matrix) <- colnames(l2fc_df)

# Number of permutations/bootstrap iterations
n_perm <- 1000

for (Human_group in Human_groups) {
  human_data_filtered <- human_l2fc_df[rownames(human_l2fc_df) %in% human_genes, ]
  animal_data_filtered <- l2fc_df[rownames(l2fc_df) %in% animal_genes, ]
  rownames(animal_data_filtered) <- human_genes[match(rownames(animal_data_filtered), animal_genes)]
  
  common_genes <- intersect(rownames(human_data_filtered), rownames(animal_data_filtered))
  human_data_filtered <- human_data_filtered[common_genes, ]
  animal_data_filtered <- animal_data_filtered[common_genes, ]
  
  human_l2fc <- human_data_filtered[[Human_group]]  # Select the relevant human group column
  animal_l2fc <- animal_data_filtered
  
  # Calculate correlation using robust statistics (bootstrapping and permutation)
  robust_results <- lapply(animal_l2fc, function(x) get_robust_stats(human_l2fc, x, n_perm = n_perm))
  
  correlation_results <- sapply(robust_results, function(x) x$r)
  pvalue_results <- sapply(animal_l2fc, function(x) cor.test(human_l2fc, x, use = "complete.obs")$p.value)
  perm_pvalue_results <- sapply(robust_results, function(x) x$p_perm)
  bootstrap_ci_low_results <- sapply(robust_results, function(x) x$ci_low)
  bootstrap_ci_high_results <- sapply(robust_results, function(x) x$ci_high)
  
  correlation_matrix[Human_group, ] <- correlation_results
  pvalue_matrix[Human_group, ] <- pvalue_results
  perm_pvalue_matrix[Human_group, ] <- perm_pvalue_results
  bootstrap_ci_low_matrix[Human_group, ] <- bootstrap_ci_low_results
  bootstrap_ci_high_matrix[Human_group, ] <- bootstrap_ci_high_results
}

correlation_matrix
pvalue_matrix
perm_pvalue_matrix
bootstrap_ci_low_matrix
bootstrap_ci_high_matrix

treatments <- setNames(coldata_df$Treatment, coldata_df$Group)
phases <- setNames(coldata_df$Phase, coldata_df$Group)
ages <- setNames(coldata_df$Age, coldata_df$Group)

# Annotation colors
treatment_colors <- treatment_color_vector
phase_colors <- c("hyperacute" = "#FFFFCC", "acute" = "#FFCC99", "DOFS" = "#E6A876", "intermediate" = "#CC9966", "chronic" = "#996633")
age_colors <- c("adult" = "#FFCCFF", "infant" = "#FFEEFF")

rownames(correlation_matrix) <- c('FILTERED_1_MTLEALL_NL', 'FILTERED_2_MTLEHS_NL', 'FILTERED_3_MTLE_NL')
rownames(pvalue_matrix) <- c('FILTERED_1_MTLEALL_NL', 'FILTERED_2_MTLEHS_NL', 'FILTERED_3_MTLE_NL')
rownames(perm_pvalue_matrix) <- c('FILTERED_1_MTLEALL_NL', 'FILTERED_2_MTLEHS_NL', 'FILTERED_3_MTLE_NL')
rownames(bootstrap_ci_low_matrix) <- c('FILTERED_1_MTLEALL_NL', 'FILTERED_2_MTLEHS_NL', 'FILTERED_3_MTLE_NL')
rownames(bootstrap_ci_high_matrix) <- c('FILTERED_1_MTLEALL_NL', 'FILTERED_2_MTLEHS_NL', 'FILTERED_3_MTLE_NL')
colnames(correlation_matrix) <- gsub("_", "-", colnames(correlation_matrix))
colnames(pvalue_matrix) <- gsub("_", "-", colnames(pvalue_matrix))
colnames(perm_pvalue_matrix) <- gsub("_", "-", colnames(perm_pvalue_matrix))
colnames(bootstrap_ci_low_matrix) <- gsub("_", "-", colnames(bootstrap_ci_low_matrix))
colnames(bootstrap_ci_high_matrix) <- gsub("_", "-", colnames(bootstrap_ci_high_matrix))
names(treatments) <- gsub("_", "-", names(treatments))
names(phases) <- gsub("_", "-", names(phases))
names(ages) <- gsub("_", "-", names(ages))
# Keep treatments with underscores to match treatment_color_vector keys

################################################################################
############## Heatmap Plotting ##############
### Bottom annotation ###
treatment_annotation <- treatments[colnames(correlation_matrix)]
treatment_annotation <- factor(treatment_annotation, levels = treatment_vector)
phase_annotation <- phases[colnames(correlation_matrix)]
phase_annotation <- factor(phase_annotation, levels = c("hyperacute", "acute", "DOFS", "intermediate", "chronic"))
age_annotation <- ages[colnames(correlation_matrix)]
age_annotation <- factor(age_annotation, levels = c("infant", "adult"))

bottom_anno <- HeatmapAnnotation(
  Treatment = treatment_annotation, 
  Phase = phase_annotation, 
  Age = age_annotation,
  col = list(
    Treatment = treatment_colors,
    Phase = phase_colors,
    Age = age_colors
  ),
  annotation_name_side = "left", # Place the annotation name on the left
  annotation_legend_param = list(
    Treatment = list(
      title = "Treatment",
      title_gp = gpar(fontsize = 12, fontfamily = "Arial"), # Arial font
      labels_gp = gpar(fontsize = 12, fontfamily = "Arial") # Arial font
    ),
    Phase = list(
      title = "Phase",
      title_gp = gpar(fontsize = 12, fontfamily = "Arial"), # Arial font
      labels_gp = gpar(fontsize = 12, fontfamily = "Arial") # Arial font
    ),
    Age = list(
      title = "Age",
      title_gp = gpar(fontsize = 12, fontfamily = "Arial"), # Arial font
      labels_gp = gpar(fontsize = 12, fontfamily = "Arial") # Arial font
    )
  ),
  gap = unit(2, "mm")
)

### color for values ###
cols <- colorRamp2(c(-max(correlation_matrix), 0, max(correlation_matrix)), c("dodgerblue", "white", "firebrick"))
cell_width <- unit(24, "mm")  # Set the cell width
cell_height <- unit(36, "mm")  # Set the cell height

### Heatmap without clustering ###
heat_plot <- Heatmap(correlation_matrix,
                     name = "Correlation", 
                     col = cols,  
                     rect_gp = gpar(col = "white", lwd = 3), 
                     cluster_rows = FALSE, 
                     cluster_columns = FALSE,
                     show_row_dend = FALSE, 
                     show_column_dend = FALSE,
                     row_names_side = "left",
                     row_names_gp = gpar(fontsize = 15,  rot = 45, fontface = "bold", fontfamily = "Arial"),
                     column_names_rot = 45,
                     column_names_gp = gpar(fontsize = 15, fontface = "bold", fontfamily = "Arial"), 
                     bottom_annotation = bottom_anno,
                     column_title = paste0("Correlation (Human vs. ", Species, " Models)"), # Add title
                     column_title_gp = gpar(fontsize = 20, fontface = "bold", fontfamily = "Arial"), # Customize font for title
                     width = cell_width * ncol(correlation_matrix),  
                     height = cell_height * nrow(correlation_matrix)) 

# Save the heatmap
file_name <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/01.Correlation/01.Correlation_WholeGene_", 
                    Species, Human_groups_name, ".png")
CairoPNG(file = file_name, width = 1600, height = 600) 
CairoPNG(file = file_name, width = 1100, height = 600) 
draw(heat_plot, heatmap_legend_side = "right", annotation_legend_side = "right", padding = unit(c(10, 10, 10, 10), "mm"))
dev.off()


output_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/01.Correlation/01.Correlation_WholeGene_", Species, Human_groups_name, ".txt")
write.table(correlation_matrix, file = output_file, sep = "\t", quote = FALSE, col.names = TRUE)

# Save additional statistics from bootstrapping and permutation testing
output_file_pvalue <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/01.Correlation/01.Correlation_WholeGene_", Species, Human_groups_name, "_pvalue.txt")
write.table(pvalue_matrix, file = output_file_pvalue, sep = "\t", quote = FALSE, col.names = TRUE)

output_file_perm_pvalue <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/01.Correlation/01.Correlation_WholeGene_", Species, Human_groups_name, "_perm_pvalue.txt")
write.table(perm_pvalue_matrix, file = output_file_perm_pvalue, sep = "\t", quote = FALSE, col.names = TRUE)

output_file_bootstrap_ci_low <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/01.Correlation/01.Correlation_WholeGene_", Species, Human_groups_name, "_bootstrap_CI_low.txt")
write.table(bootstrap_ci_low_matrix, file = output_file_bootstrap_ci_low, sep = "\t", quote = FALSE, col.names = TRUE)

output_file_bootstrap_ci_high <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/01.Correlation/01.Correlation_WholeGene_", Species, Human_groups_name, "_bootstrap_CI_high.txt")
write.table(bootstrap_ci_high_matrix, file = output_file_bootstrap_ci_high, sep = "\t", quote = FALSE, col.names = TRUE)