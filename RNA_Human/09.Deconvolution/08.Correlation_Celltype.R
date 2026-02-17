# Load necessary libraries
library(ComplexHeatmap)
library(Cairo)
library(circlize)  # For colorRamp2


# Choose
Dataset_H <- 'GSE160189_186538_Integration'
Dataset_M <- 'GSE185862'

row_order <- c("ExN_H", "InN_H", "Astro_H", "Micro_H", "Oligo_H", "OPC_H", "Endo_H") 
column_order <- c("ExN_M", "InN_M", "Astro_M", "Micro_M", "Oligo_M", "OPC_M", "Endo_M")


#### Load dataset #### 
human_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/03.CIBERSORT_", Dataset_H, 
                     "_01Sig/", Dataset_H, "_scRNA_matrix_DEGs_CellType1/CIBERSORTx_",
                     Dataset_H, "_scRNA_matrix_DEGs_CellType1_inferred_phenoclasses.CIBERSORTx_",
                     Dataset_H, "_scRNA_matrix_DEGs_CellType1_inferred_refsample.bm.K999.txt")
human_data <- read.table(human_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(human_data) <- human_data$NAME
human_data <- human_data[, -1]

mouse_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_", Dataset_M, 
                     "_01Sig/", Dataset_M, "_10X_scRNA_matrix_DEGs_CellType1/CIBERSORTx_",
                     Dataset_M, "_10X_scRNA_matrix_DEGs_CellType1_inferred_phenoclasses.CIBERSORTx_",
                     Dataset_M, "_10X_scRNA_matrix_DEGs_CellType1_inferred_refsample.bm.K999.txt")
mouse_data <- read.table(mouse_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(mouse_data) <- mouse_data$NAME
mouse_data <- mouse_data[, -1]

#### Load orthologs #### 
orthologs_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_mm_17128.txt"
ortholog_one2one <- read.table(orthologs_path, sep = "\t", header = TRUE, fill = TRUE)
# Remove duplicates in human and mouse gene names
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one$Gene.name), ]
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one$Mouse.gene.name), ]
# Get human and mouse gene mappings
human_genes <- ortholog_one2one$Gene.name
mouse_genes <- ortholog_one2one$Mouse.gene.name


#### Modification of genes #### 
# Match mouse genes to their human orthologs using the ortholog file
mouse_to_human_mapping <- human_genes[match(rownames(mouse_data), mouse_genes)]
valid_rows <- !is.na(mouse_to_human_mapping) & !duplicated(mouse_to_human_mapping)
mouse_data_renamed <- mouse_data[valid_rows, ]
rownames(mouse_data_renamed) <- mouse_to_human_mapping[valid_rows]
# Examine and keep only shared genes between the two datasets
shared_genes <- intersect(rownames(human_data), rownames(mouse_data_renamed))
human_data_filtered <- human_data[shared_genes, ]
human_data_filtered <- human_data_filtered[, !(colnames(human_data_filtered) %in% "ExN_EC")]
mouse_data_filtered <- mouse_data_renamed[shared_genes, ]
mouse_data_filtered <- mouse_data_filtered[, !(colnames(mouse_data_filtered) %in% "ExN_NA")]


#### Correlation calculation #### 
correlation_matrix <- cor(human_data_filtered, mouse_data_filtered, use = "complete.obs")


#### Heatmap ####
cols <- colorRamp2(c(-max(correlation_matrix), 0, max(correlation_matrix)), c("dodgerblue", "white", "firebrick"))

# cluster_1 : width = 500, height = 400
# cluster_2 : width = 600, height = 500

# with clustering
heat_plot <- Heatmap(correlation_matrix,
                     name = "Correlation", 
                     col = cols, 
                     rect_gp = gpar(col = "white", lwd = 2),
                     row_title = paste0("Human Cell Types_", Dataset_H),
                     row_title_gp = gpar(fontsize = 15, fontface = "bold"),
                     column_title = paste0("Mouse Cell Types_", Dataset_M),
                     column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                     cluster_rows = TRUE,  # Disable row clustering
                     cluster_columns = TRUE,  # Disable column clustering
                     show_row_dend = TRUE,  # Disable row dendrogram
                     show_column_dend = TRUE,  # Disable column dendrogram
                     row_names_gp = gpar(fontsize = 15),  # Set font size for row names
                     column_names_gp = gpar(fontsize = 15))  # Set font size for column names

# Save the heatmap
file_name <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/07.Correlation_Celltype/", 
                    Dataset_M, "_", Dataset_H, "_cluster.png")
CairoPNG(file = file_name, width = 700, height = 600)
draw(heat_plot, heatmap_legend_side = "right", annotation_legend_side = "right", padding = unit(c(10, 10, 10, 40), "mm"))
dev.off()

# without clustering
heat_plot <- Heatmap(correlation_matrix,
                     name = "Correlation", 
                     col = cols, 
                     rect_gp = gpar(col = "white", lwd = 2),
                     row_title = paste0("Human Cell Types_", Dataset_H),
                     row_title_gp = gpar(fontsize = 15, fontface = "bold"),
                     column_title = paste0("Mouse Cell Types_", Dataset_M),
                     column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                     cluster_rows = FALSE,  # Disable row clustering
                     cluster_columns = FALSE,  # Disable column clustering
                     show_row_dend = FALSE,  # Disable row dendrogram
                     show_column_dend = FALSE,  # Disable column dendrogram
                     row_names_gp = gpar(fontsize = 15),  # Set font size for row names
                     column_names_gp = gpar(fontsize = 15),
                     row_order = row_order,
                     column_order = column_order) 

# Save the heatmap
file_name <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/07.Correlation_Celltype/", 
                    Dataset_M, "_", Dataset_H, ".png")
CairoPNG(file = file_name, width = 700, height = 600)
draw(heat_plot, heatmap_legend_side = "right", annotation_legend_side = "right", padding = unit(c(10, 10, 10, 40), "mm"))
dev.off()

