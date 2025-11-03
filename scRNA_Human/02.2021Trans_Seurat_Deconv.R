library(Seurat)
library(SingleCellExperiment)
library(Matrix)

setwd('/home/joonho345/3_RNA/scRNA_Human/')
GSE <- '2021Trans'

#### Create Seurat Object ####
load("/home/joonho345/3_RNA/scRNA_Human/Raw_Data/2021Trans/SCE_HPC-n3_tran-etal.rda")
sce.hpc.tran # SingleCellExperiment object

hippo_trans <- as.Seurat(sce.hpc.tran, counts = "counts", data = "logcounts")
hippo_trans # 33538 features across 10268 samples
unique(hippo_trans@meta.data$orig.ident) # Levels: donor1 donor2 donor3
unique(hippo_trans@meta.data$region) # Levels: hpc
unique(hippo_trans@meta.data$cellType)
# Excit_A, Excit_B, Excit_C, Excit_D, Excit_E, Excit_F, Excit_H, Excit_G
# Inhib_A, Inhib_B, Inhib_C, Inhib_D
# Astro_A, Astro_B
# Oligo, OPC, OPC_COP
# Micro
# Tcell, Mural
# drop.doublet, drop.lowNTx_A, drop.lowNTx_B

#### Quality Control ####
## visualization
CairoPNG('/home/joonho345/3_RNA/scRNA_Human/02.2021Trans_Seurat_Deconv/01.QC/2021Trans_VlnPlot.png', width = 800, height = 600)
VlnPlot(hippo_trans, features = c("nFeature_originalexp", "nCount_originalexp"), ncol = 2)
dev.off()

CairoPNG('/home/joonho345/3_RNA/scRNA_Human/02.2021Trans_Seurat_Deconv/01.QC/2021Trans_FeatureScatter.png', width = 800, height = 600)
FeatureScatter(hippo_trans, feature1 = "nCount_originalexp", feature2 = "nFeature_originalexp")
dev.off()

## cell QC - Done
## feature QC - Done
## Normalization - Done


#### Clustering ####
exn_cells <- c("Excit_A", "Excit_B", "Excit_C", "Excit_D", "Excit_E", "Excit_F", "Excit_G", "Excit_H")
inn_cells <- c("Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D")
astro_cells <- c("Astro_A", "Astro_B")
oligo_cells <- c("Oligo", "OPC", "OPC_COP")
micro_cells <- c("Micro")
imm_cells <- c("Tcell", "Mural")
drop_cells <- c("drop.doublet", "drop.lowNTx_A", "drop.lowNTx_B")

## cluster_1 ##
hippo_trans$cluster_1 <- NA
hippo_trans@meta.data[hippo_trans@meta.data$cellType %in% exn_cells, ]$cluster_1 <- "ExN"
hippo_trans@meta.data[hippo_trans@meta.data$cellType %in% inn_cells, ]$cluster_1 <- "InN"
hippo_trans@meta.data[hippo_trans@meta.data$cellType %in% astro_cells, ]$cluster_1 <- "Astro"
hippo_trans@meta.data[hippo_trans@meta.data$cellType %in% oligo_cells, ]$cluster_1 <- "Oligo"
hippo_trans@meta.data[hippo_trans@meta.data$cellType %in% micro_cells, ]$cluster_1 <- "Micro"
hippo_trans@meta.data[hippo_trans@meta.data$cellType %in% imm_cells, ]$cluster_1 <- "Imm"
hippo_trans@meta.data[hippo_trans@meta.data$cellType %in% drop_cells, ]$cluster_1 <- "Drop"

## filter cells ##
hippo_trans <- subset(hippo_trans, subset = cluster_1 != "Imm" & cluster_1 != "Drop")
hippo_trans # 33538 features across 10070 samples
unique(hippo_trans@meta.data$cluster_1)
sum(is.na(hippo_trans@meta.data$cluster_1)) # no NA
# cluster_1: ExN  InN Astro Oligo Micro


#### Visualization of clusters ####
# Clusters
metadata_fields <- c("orig.ident", "region", "cellType", "cluster_1")
reduction_types <- c("PCA_correct", "PCA_opt", "UMAP", "TSNE")
plot_directory <- "/home/joonho345/3_RNA/scRNA_Human/02.2021Trans_Seurat_Deconv/02.Cluster/"

reduc <- "UMAP"
field <- "cluster_1"
plot_path <- paste0(plot_directory, "2021Trans_", reduc, "plot_", field, ".png")

CairoPNG(plot_path, width = 800, height = 600)
DimPlot(hippo_trans, reduction = reduc, group.by = field, raster = FALSE, label = T)
dev.off()

# Marker genes
ExN.markers = c("SLC17A7", "PROX1", "CFAP299", "SYN3", "HGF", "GRIK1")
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.2021Trans_Seurat_Deconv/03.Markers/2021Trans_UMAPplot_ExN.png", 
         width = 800, height = 800)
FeaturePlot(hippo_trans, features=ExN.markers, reduction="UMAP", raster=FALSE) ;
dev.off()

InN.markers = c("GAD1", "PVALB", "SST", "VIP")
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.2021Trans_Seurat_Deconv/03.Markers/2021Trans_UMAPplot_InN.png", 
         width = 800, height = 800)
FeaturePlot(hippo_trans, features=InN.markers, reduction="UMAP", raster=FALSE) ;
dev.off()

Astro.markers = c("AQP4", "GFAP", "AQP1", "APOE")
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.2019CA_Seurat_Deconv/03.Markers/2021Trans_UMAPplot_Astro.png", 
         width = 800, height = 800)
FeaturePlot(CA, features=Astro.markers, reduction="umap", raster=FALSE) ;
dev.off()

Oligo.markers = c("MOBP", "SOX10", "OPALIN", "PDGFRA", "GPR17")
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.2021Trans_Seurat_Deconv/03.Markers/2021Trans_UMAPplot_Oligo.png", 
         width = 800, height = 800)
FeaturePlot(CA, features=Oligo.markers, reduction="umap", raster=FALSE) ;
dev.off()

Micro.markers = c("AIF1", "PTPRC", "SLC1A3")
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.2021Trans_Seurat_Deconv/03.Markers/2021Trans_UMAPplot_Micro.png", 
         width = 800, height = 800)
FeaturePlot(CA, features=Micro.markers, reduction="umap", raster=FALSE) ;
dev.off()


