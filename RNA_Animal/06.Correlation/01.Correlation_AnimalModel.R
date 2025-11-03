# Load necessary libraries
library(ComplexHeatmap)
library(Cairo)


##
Species <- 'Mouse'
orthologs_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_mm_17128.txt"
l2fc_df <- DESeq_l2fc_df_M # import data from 00.Matrix_M.R
coldata_df_M <- DESeq_coldata_df_M # import data from 00.Matrix_M.R

ortholog_one2one <- read.table(orthologs_path, sep = "\t", header = TRUE, fill = TRUE)
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one$Gene.name), ]
Species_Gene_name <- paste0(Species, '.gene.name')
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one[[Species_Gene_name]]), ]
human_genes <- ortholog_one2one$Gene.name
animal_genes <- ortholog_one2one[[Species_Gene_name]]

animal_data_filtered <- l2fc_df[rownames(l2fc_df) %in% animal_genes, ]
rownames(animal_data_filtered) <- human_genes[match(rownames(animal_data_filtered), animal_genes)]
M_l2fc_df <- animal_data_filtered

##
Species <- 'Rat'
orthologs_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_rn_16211.txt"
l2fc_df <- DESeq_l2fc_df_R # import data from 00.Matrix_M.R
coldata_df_R <- DESeq_coldata_df_R # import data from 00.Matrix_M.R

ortholog_one2one <- read.table(orthologs_path, sep = "\t", header = TRUE, fill = TRUE)
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one$Gene.name), ]
Species_Gene_name <- paste0(Species, '.gene.name')
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one[[Species_Gene_name]]), ]
human_genes <- ortholog_one2one$Gene.name
animal_genes <- ortholog_one2one[[Species_Gene_name]]

animal_data_filtered <- l2fc_df[rownames(l2fc_df) %in% animal_genes, ]
rownames(animal_data_filtered) <- human_genes[match(rownames(animal_data_filtered), animal_genes)]
R_l2fc_df <- animal_data_filtered


## M+R Combined df 
shared_genes <- intersect(rownames(M_l2fc_df), rownames(R_l2fc_df))
M_l2fc_df_filtered <- M_l2fc_df[shared_genes, ]
R_l2fc_df_filtered <- R_l2fc_df[shared_genes, ]
M_l2fc_df_filtered <- M_l2fc_df_filtered[order(rownames(M_l2fc_df_filtered)), ]
R_l2fc_df_filtered <- R_l2fc_df_filtered[order(rownames(R_l2fc_df_filtered)), ]
combined_df <- cbind(M_l2fc_df_filtered, R_l2fc_df_filtered)

Group_order <- c(
  "M_KAI_ST_BO_O_HA", "M_KAI_ST_BO_O_AC", "M_KAI_ST_BO_O_IM", "M_KAI_ST_BO_O_CR",
  "M_KAI_ST_IPSI_Y_AC", "M_KAI_ST_IPSI_Y_IM", "M_KAI_ST_IPSI_Y_CR", 
  "M_KAI_ST_CON_Y_AC", "M_KAI_ST_CON_Y_IM", "M_KAI_ST_CON_Y_CR", 
  "M_KAI_IP_O_HA", "R_KAI_IP_Y_CR", "R_KAI_SUB_I_IM", 
  "M_PILO_IP_Y_HA", "M_PILO_IP_Y_AC", "M_PILO_IP_Y_IM", 
  "M_PILO_IP_O_HA", "M_PILO_IP_O_IM", "M_PILO_IP_O_CR", 
  "R_PILO_IP_Y_CR", "R_PPS_Y_AC", "R_PPS_Y_DOFS", "R_PPS_Y_IM", "R_PPS_Y_CR", 
  "R_AMG_IPSI_Y_CR", "R_TBI_IPSI_Y_CR"
)
Group_order <- Group_order[Group_order %in% colnames(combined_df)]
combined_df <- combined_df[, Group_order, drop = FALSE]

# Calculate correlation matrix between mouse groups
correlation_matrix <- cor(combined_df, use = "complete.obs")
rownames(correlation_matrix) <- gsub("_", "-", rownames(correlation_matrix))
colnames(correlation_matrix) <- gsub("_", "-", colnames(correlation_matrix))

# Define color palette for heatmap
cols <- colorRamp2(c(-max(correlation_matrix), 0, max(correlation_matrix)), c("dodgerblue", "white", "firebrick"))

#### Heatmap without clustering ####
heat_plot <- Heatmap(correlation_matrix,
                     name = "Correlation", 
                     col = cols,  # Use the color scale
                     rect_gp = gpar(col = "white", lwd = 2),  # Style the grid
                     row_title = paste0(Species, "Groups"),
                     column_title = paste0(Species, "Groups"),
                     column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                     cluster_rows = FALSE,  # Disable row clustering
                     cluster_columns = FALSE,  # Disable column clustering
                     show_row_dend = FALSE,  # Disable row dendrogram
                     show_column_dend = FALSE,  # Disable column dendrogram
                     row_names_gp = gpar(fontsize = 12),  # Set font size for row names
                     column_names_gp = gpar(fontsize = 12))  # Set font size for column names

# Save the heatmap
file_name <- paste0("/home/joonho345/3_RNA/RNA_Animal/06.Correlation/01.Correlation_Models.png")
CairoPNG(file = file_name, width = 1600, height = 700) 
draw(heat_plot, heatmap_legend_side = "right", annotation_legend_side = "right", padding = unit(c(10, 10, 10, 40), "mm"))  # Adjust padding for better visualization
dev.off()

#### Heatmap with clustering ####
heat_plot <- Heatmap(correlation_matrix,
                     name = "Correlation", 
                     col = cols,  # Use the color scale
                     rect_gp = gpar(col = "white", lwd = 2),  # Style the grid
                     row_title = paste0(Species, "Groups"),
                     column_title = paste0(Species, "Groups"),
                     column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                     cluster_rows = TRUE,  # Disable row clustering
                     cluster_columns = TRUE,  # Disable column clustering
                     show_row_dend = TRUE,  # Disable row dendrogram
                     show_column_dend = TRUE,  # Disable column dendrogram
                     row_names_gp = gpar(fontsize = 12),  # Set font size for row names
                     column_names_gp = gpar(fontsize = 12))  # Set font size for column names

# Save the heatmap
file_name <- paste0("/home/joonho345/3_RNA/RNA_Animal/06.Correlation/01.Correlation_Models_cluster.png")
CairoPNG(file = file_name, width = 1600, height = 700) 
draw(heat_plot, heatmap_legend_side = "right", annotation_legend_side = "right", padding = unit(c(10, 10, 10, 40), "mm"))  # Adjust padding for better visualization
dev.off()

