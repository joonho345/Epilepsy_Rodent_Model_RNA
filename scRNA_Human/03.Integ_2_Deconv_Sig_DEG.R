# Load necessary libraries
library(Seurat)
library(Matrix)
library(dplyr)

# Set working directory
setwd('/home/joonho345/1_Epilepsy_RNA/scRNA_Human/')


#########################################################
#### Load Seurat Object ####
seurat_combined_o <- readRDS("/data/project/1_Epilepsy_RNA/scRNA_Human/Deconv_objects/GSE160189_186538_orig_2_RPCA_1.rds")
DefaultAssay(seurat_combined_o) <- "RNA"
integration_name <- 'GSE160189_186538_Integration'
type <- 'o'
method <- 'rpca'
seurat_combined <- seurat_combined_o


#########################################################
#### Differentially Expressed Genes (DEGs) for matrix
## Find Markers for CellType1 ##
Idents(seurat_combined) <- seurat_combined@meta.data$CellType1
markers_CellType1 <- FindAllMarkers(seurat_combined, 
                                     only.pos = TRUE, 
                                     min.pct = 0.25, 
                                     logfc.threshold = 0.25,
                                     verbose = FALSE)

# Select top markers per cell type (e.g., top 200 based on avg_log2FC and p_val_adj)
top_markers_per_celltype <- 200
markers_CellType1_top <- markers_CellType1 %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), p_val_adj) %>%
  top_n(n = min(top_markers_per_celltype, n()), wt = avg_log2FC)

## Union of markers for CellType1 ##
union_markers_CellType1 <- unique(markers_CellType1_top$gene)


#### Save marker information ####
output_base_dir <- "/home/joonho345/1_Epilepsy_RNA/scRNA_Human/02.GSE160189_186538_Integration_Seurat_Deconv"
markers_output_dir <- paste0(output_base_dir, "/05.Markers/")
dir.create(markers_output_dir, recursive = TRUE, showWarnings = FALSE)

write.table(markers_CellType1_top, 
            file = paste0(markers_output_dir, integration_name, "_markers_CellType1_top", top_markers_per_celltype, ".txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(Gene = union_markers_CellType1), 
            file = paste0(markers_output_dir, integration_name, "_union_markers_CellType1.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## Read union marker gene sets from saved files ##
union_markers_CellType1 <- read.table(paste0(markers_output_dir, integration_name, "_union_markers_CellType1.txt"),
                                       sep = "\t", 
                                       header = FALSE, 
                                       stringsAsFactors = FALSE)$V1


#########################################################
#### Extract Expression Matrix using DEGs ##
# Get expression matrix as sparse matrix (keep it sparse to save memory)
expression_matrix <- seurat_combined@assays$RNA@layers$data
features <- rownames(seurat_combined)
rownames(expression_matrix) <- features

## Output files using CellType1 DEGs ##
# Filter expression matrix to include only DEG markers for CellType1 (while still sparse)
# Check which markers are actually present in the expression matrix
markers_in_matrix_CellType1 <- union_markers_CellType1[union_markers_CellType1 %in% features]
if (length(markers_in_matrix_CellType1) < length(union_markers_CellType1)) {
  missing_markers <- setdiff(union_markers_CellType1, markers_in_matrix_CellType1)
  warning(length(missing_markers), " CellType1 markers not found in expression matrix")
}

# Filter sparse matrix first, then convert only the subset to data.frame
expression_matrix_CellType1 <- expression_matrix[features %in% union_markers_CellType1, ]
cat("Filtered matrix for CellType1:", dim(expression_matrix_CellType1), "\n")

# Convert only the filtered subset to data.frame (much smaller!)
expression_df_CellType1 <- as.data.frame(expression_matrix_CellType1, check.names = FALSE)
rownames(expression_df_CellType1) <- rownames(expression_matrix_CellType1)

# Prepare data frame with Gene column
expression_df_CellType1$Gene <- rownames(expression_df_CellType1)
expression_df_CellType1 <- expression_df_CellType1[, c(ncol(expression_df_CellType1), 1:(ncol(expression_df_CellType1) - 1))]
# Output file without column names (for appending headers)
output_path_CellType1 <- paste0(output_base_dir, "/", integration_name, "_scRNA_matrix_no_colnames_CellType1.txt")
write.table(expression_df_CellType1, file = output_path_CellType1, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# CellType1 output
celltype_annotations_CellType1 <- as.character(seurat_combined@meta.data$CellType1)
header_CellType1 <- c("Gene", celltype_annotations_CellType1)
formatted_output_path_CellType1 <- paste0(output_base_dir, "/", integration_name, "_scRNA_matrix_DEGs_CellType1.txt")
writeLines(paste(header_CellType1, collapse = "\t"), con = formatted_output_path_CellType1)
file.append(formatted_output_path_CellType1, output_path_CellType1)
