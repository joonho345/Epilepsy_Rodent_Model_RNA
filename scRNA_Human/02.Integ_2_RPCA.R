library(Seurat)
library(Matrix)
library(Cairo)

#########################################################
#### Integration of two human datasets ####
seurat_human_1 <- GSE160189_hippo_int
seurat_human_2 <- GSE186538_hippo_int
integration_name <- 'GSE160189_186538_Integration'

## Find common genes between the two datasets ####
human_features_1 <- rownames(seurat_human_1)
human_features_2 <- rownames(seurat_human_2)
common_genes <- intersect(human_features_1, human_features_2)
cat("Total genes in dataset 1 (GSE160189):", length(human_features_1), "\n")
cat("Total genes in dataset 2 (GSE186538):", length(human_features_2), "\n")
cat("Common genes between datasets:", length(common_genes), "\n")
seurat_human_1 <- subset(seurat_human_1, features = common_genes)
seurat_human_2 <- subset(seurat_human_2, features = common_genes)

# species column
seurat_human_1@meta.data$species <- 'human'
seurat_human_2@meta.data$species <- 'human'
# group column
seurat_human_1@meta.data$group <- 'H_case'
seurat_human_2@meta.data$group <- 'H_control'


#########################################################
#### Integration key: orig.ident ####
human_list_1 <- SplitObject(seurat_human_1, split.by = "orig.ident")
human_list_2 <- SplitObject(seurat_human_2, split.by = "orig.ident")
combined_list_orig_2 <- c(human_list_1, human_list_2)
# Normalize and find variable features for each split
combined_list_orig_2 <- lapply(combined_list_orig_2, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
  return(x)
})

#### Integration of objects RPCA #### 
combined_list_orig_2_RPCA <- combined_list_orig_2
RPCA_reference <- c(1,2)

features <- SelectIntegrationFeatures(object.list = combined_list_orig_2_RPCA)
combined_list_orig_2_RPCA <- lapply(X = combined_list_orig_2_RPCA, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = combined_list_orig_2_RPCA, reference = RPCA_reference, reduction = "rpca", 
                                  dims = 1:30)
seurat_combined_orig_2_RPCA <- IntegrateData(anchorset = anchors, dims = 1:30)
seurat_combined_orig_2_RPCA <- ScaleData(seurat_combined_orig_2_RPCA, verbose = FALSE)
seurat_combined_orig_2_RPCA <- RunPCA(seurat_combined_orig_2_RPCA, verbose = FALSE)
seurat_combined_orig_2_RPCA <- RunUMAP(seurat_combined_orig_2_RPCA, dims = 1:30)


#########################################################
#### Cell type trimming ####
seurat_combined_orig_2_RPCA <- subset(seurat_combined_orig_2_RPCA, subset = CellType != "Imm" & CellType != 'CR')
exn_cells <- c("Pyr1", "Pyr2", "ExN", "Den.Gyr1", "Den.Gyr2", "Den.Gyr3")
inn_cells <- c("In1", "In2", "In3", "InN")
astro_cells <- c("Astro1", "Astro2", "Astro3", "Astro")
micro_cells <- c("Micro1", "Micro2", "Micro3", "Micro")
oligo_cells <- c("Olig1", "Olig2", "Olig3", "Olig4", "Olig5", "Oligo")
opc_cells <- c("OPC1", "OPC2", "OPC3", "OPC4", "OPC_H")
endo_cells <- c("Endo", "Vascular")

seurat_combined_orig_2_RPCA$CellType1 <- NA
seurat_combined_orig_2_RPCA@meta.data[seurat_combined_orig_2_RPCA@meta.data$CellType %in% astro_cells, "CellType1"] <- "Astro_H"
seurat_combined_orig_2_RPCA@meta.data[seurat_combined_orig_2_RPCA@meta.data$CellType %in% endo_cells, "CellType1"] <- "Endo_H"
seurat_combined_orig_2_RPCA@meta.data[seurat_combined_orig_2_RPCA@meta.data$CellType %in% exn_cells, "CellType1"] <- "ExN_H"
seurat_combined_orig_2_RPCA@meta.data[seurat_combined_orig_2_RPCA@meta.data$CellType %in% inn_cells, "CellType1"] <- "InN_H"
seurat_combined_orig_2_RPCA@meta.data[seurat_combined_orig_2_RPCA@meta.data$CellType %in% micro_cells, "CellType1"] <- "Micro_H"
seurat_combined_orig_2_RPCA@meta.data[seurat_combined_orig_2_RPCA@meta.data$CellType %in% oligo_cells, "CellType1"] <- "Oligo_H"
seurat_combined_orig_2_RPCA@meta.data[seurat_combined_orig_2_RPCA@meta.data$CellType %in% opc_cells, "CellType1"] <- "OPC_H"

table(seurat_combined_orig_2_RPCA$CellType1)
table(seurat_combined_orig_2_RPCA$CellType)


#########################################################
#### Normalize RNA assay after subsetting ####
DefaultAssay(seurat_combined_orig_2_RPCA) <- "RNA"
seurat_combined_orig_2_RPCA <- NormalizeData(seurat_combined_orig_2_RPCA, verbose = FALSE)
seurat_combined_orig_2_RPCA <- JoinLayers(seurat_combined_orig_2_RPCA)


#########################################################
#### Visualization ####
output_base_dir <- "/home/joonho345/1_Epilepsy_RNA/scRNA_Human/02.GSE160189_186538_Integration_Seurat_Deconv"
dir.create(paste0(output_base_dir, "/04_01.integration_orig"), recursive = TRUE, showWarnings = FALSE)

CairoPNG(paste0(output_base_dir, '/04_01.integration_orig/umaprpca_celltype1.png'), width = 800, height = 600)
DimPlot(seurat_combined_orig_2_RPCA, reduction = "umap", group.by = "CellType1", shuffle = TRUE, label = TRUE)
dev.off()
CairoPNG(paste0(output_base_dir, '/04_01.integration_orig/umaprpca_group.png'), width = 800, height = 600)
DimPlot(seurat_combined_orig_2_RPCA, reduction = "umap", group.by = "group", shuffle = TRUE, label = TRUE)
dev.off()
markers = c("SLC17A7", "GAD1", "GFAP", "PTPRC", "MOBP", "PDGFRA", "PECAM1")
CairoPNG(paste0(output_base_dir, '/04_01.integration_orig/umaprpca_marker.png'), width = 800, height = 800)
FeaturePlot(seurat_combined_orig_2_RPCA, features=markers, reduction="umap", raster=FALSE) ;
dev.off()


#########################################################
#### Merge the two objects without integration ####
#########################################################
seurat_merged <- merge(seurat_human_1, y = seurat_human_2)
seurat_merged <- NormalizeData(seurat_merged)
seurat_merged <- FindVariableFeatures(seurat_merged)
seurat_merged <- ScaleData(seurat_merged, verbose = FALSE)
seurat_merged <- RunPCA(seurat_merged, npcs = 30, verbose = FALSE)
seurat_merged <- RunUMAP(seurat_merged, dims = 1:30)

CairoPNG(paste0(output_base_dir, '/04_00.bf_integration/umap_group.png'), width = 800, height = 600)
DimPlot(seurat_merged, reduction = "umap", group.by = "group", shuffle = TRUE, label = TRUE)
dev.off()
