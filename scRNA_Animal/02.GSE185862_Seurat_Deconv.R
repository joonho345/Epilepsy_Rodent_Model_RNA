library(Seurat)
library(Matrix)
library(Cairo)

setwd('/home/joonho345/3_RNA/scRNA_Animal/')
GSE <- 'GSE185862'


#### Create Seurat Object ####
file_path <- "/home/joonho345/3_RNA/scRNA_Animal/Raw_Data/Animal_GSE185862/Seurat.ss.rda"
load(file_path)
ss.seurat[["percent.mt"]] <- PercentageFeatureSet(ss.seurat, pattern = "^mt-")
print(ss.seurat)
# Cell: 73363 / Gene: 42055


#### Quality Control ####
# visualization
CairoPNG('/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/01.QC/GSE185862_VlnPlot.png', width = 800, height = 600)
VlnPlot(ss.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

CairoPNG('/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/01.QC/GSE185862_FeatureScatter.png', width = 800, height = 600)
plot1 <- FeatureScatter(ss.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ss.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

## cell QC
ss.seurat <- subset(ss.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 15000 & percent.mt < 10)
ss.seurat
# Cell: 73339 / Gene: 42055 

## feature QC
min_cells <- 0.01 * ncol(ss.seurat)  # Example: 1% of the total cells
counts_matrix <- ss.seurat@assays$RNA@counts
genes_to_keep <- rowSums(counts_matrix > 0) > min_cells
length(genes_to_keep) 
# Cell: 73339 / Gene: 42055 

## Normalization
ss.seurat <- NormalizeData(ss.seurat, scale.factor = 10000)
ss.seurat@assays$RNA@counts[1:100,1:100]
ss.seurat@assays$RNA@data[1:100,1:100]
ss.seurat <- ScaleData(ss.seurat, features = rownames(ss.seurat))
ss.seurat@assays$RNA@scale.data[1:100,1:100]

#### (hip_subset) cells from hippocampus (not CTX) ####
unique(ss.seurat@meta.data$region_label)
hip_subset <- subset(ss.seurat, subset = region_label == "HIP")
unique(hip_subset@meta.data$region_label)
print(hip_subset)
# Cell: 6564 / Gene: 42055

#### (hip_subset) Cluster_1 & Cluster_2 ####
# create cluster_1 and cluster_2 matched with Human 2022 Neuron data regarding cell types
unique(hip_subset@meta.data$cluster_label)
unique(hip_subset@meta.data$class_label)
# Define vectors for cluster_1
exn_cells <- c("Glutamatergic")
exn_dg_cells <- c("361_DG", "362_DG", "363_DG", "364_DG")
exn_ca_cells <- c("340_CA1", "338_CA1", "333_CA1-ProS", "344_CA1", "356_CA3-do", "353_CA3-ve", "336_CA1-ve", 
                     "337_CA1", "332_CA1-ProS", "339_CA1", "347_CA1-do", "348_CA1-do", "351_CA3-ve", "360_CA2-IG-FC", 
                     "346_CA1-do", "352_CA3-ve", "355_CA3-ve", "341_CA1", "345_CA1", "357_CA3-do", "354_CA3-ve", 
                     "359_CA2-IG-FC", "342_CA1", "335_CA1-ve", "358_CA3-do", "334_CA1-ve", "343_CA1", "330_CA1-ProS", 
                     "329_CA1-ProS")
exn_sub_cells <- c("318_SUB", "321_SUB")
exn_unassigned_cells <- c("156_L2/3 IT ENTl", "159_L2/3 IT AI", "173_L2/3 IT ProS", "271_NP SUB", 
                          "272_NP SUB", "283_L6 CT CTX", "290_L6 CT CTX", "310_L6b RHP", 
                          "323_ProS", "324_ProS", "350_Mossy")
inn_cells <- c("GABAergic")
oligo_cells <- c("368_Oligo")
astro_cells <- c("378_Astro", "377_Astro", "376_Astro")
micro_cells <- c("387_Micro-PVM", "386_Micro-PVM", "388_Micro-PVM")

hip_subset@meta.data$cluster_2[is.na(hip_subset@meta.data$cluster_2)] <- "Unassigned"

# Assign cluster_1 labels
hip_subset@meta.data$cluster_1 <- NA
hip_subset@meta.data$cluster_1[hip_subset@meta.data$class_label %in% exn_cells] <- "ExN"
hip_subset@meta.data$cluster_1[hip_subset@meta.data$class_label %in% inn_cells] <- "InN"
hip_subset@meta.data$cluster_1[hip_subset@meta.data$cluster_label %in% oligo_cells] <- "Oligo"
hip_subset@meta.data$cluster_1[hip_subset@meta.data$cluster_label %in% astro_cells] <- "Astro"
hip_subset@meta.data$cluster_1[hip_subset@meta.data$cluster_label %in% micro_cells] <- "Micro"
# Assign cluster_2 labels
hip_subset@meta.data$cluster_2 <- NA
hip_subset@meta.data$cluster_2[hip_subset@meta.data$class_label == "Glutamatergic" & 
                                 hip_subset@meta.data$cluster_label %in% exn_dg_cells] <- "ExN_DG"
hip_subset@meta.data$cluster_2[hip_subset@meta.data$class_label == "Glutamatergic" & 
                                 hip_subset@meta.data$cluster_label %in% exn_ca_cells] <- "ExN_CA"
hip_subset@meta.data$cluster_2[hip_subset@meta.data$class_label == "Glutamatergic" & 
                                 hip_subset@meta.data$cluster_label %in% exn_sub_cells] <- "ExN_SUB"
hip_subset@meta.data$cluster_2[hip_subset@meta.data$class_label == "Glutamatergic" & 
                                 hip_subset@meta.data$cluster_label %in% exn_unassigned_cells] <- "ExN_NA"
hip_subset@meta.data$cluster_2[hip_subset@meta.data$class_label %in% inn_cells] <- "InN"
hip_subset@meta.data$cluster_2[hip_subset@meta.data$cluster_label %in% oligo_cells] <- "Oligo"
hip_subset@meta.data$cluster_2[hip_subset@meta.data$cluster_label %in% astro_cells] <- "Astro"
hip_subset@meta.data$cluster_2[hip_subset@meta.data$cluster_label %in% micro_cells] <- "Micro"

table(hip_subset@meta.data$cluster_label)
table(hip_subset@meta.data$cluster_1)
table(hip_subset@meta.data$cluster_2)


############### (hip_subset) Clustering & Cell annotation ############### 
#### FindVariableFeatures ####
hip_subset <- FindVariableFeatures(hip_subset, nfeatures = 2000)
CairoPNG("/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/01.QC/GSE185862_VariableFeaturePlot.png", 
         width = 800, height = 600)
VariableFeaturePlot(hip_subset)
dev.off()

#### Clustering ####
## PCA ##
hip_subset <- RunPCA(hip_subset, features = VariableFeatures(hip_subset))
# Decide the number of PC dimensions
# ElbowPlot
CairoPNG("/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/02.Cluster/GSE185862_ElbowPlot.png", 
         width = 800, height = 600)
ElbowPlot(hip_subset, ndims = 50)
dev.off()
# JackStrawPlot
#hip_subset <- JackStraw(hip_subset, dims = 50)
#hip_subset <- ScoreJackStraw(hip_subset, dims = 1:50)
#CairoPNG("/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/02.Cluster/GSE185862_JackStrawPlot.png", 
#         width = 800, height = 600)
#JackStrawPlot(hip_subset, dims = 1:50)
#dev.off()
# Set the number of dimensions
num_dims <- 22
message("Number of dimensions chosen: ", num_dims)

## Clustering ##
hip_subset <- FindNeighbors(hip_subset, dims = 1:num_dims)
hip_subset@graphs$RNA_nn@x
hip_subset@graphs$RNA_snn@x
hip_subset <- FindClusters(hip_subset, resolution = 0.4) # 0.4 ~ 1.2
hip_subset@meta.data$seurat_clusters

## UMAP ##
hip_subset <- RunUMAP(hip_subset, reduction="pca", dims = 1:num_dims)

## TSNE ##
hip_subset <- RunTSNE(hip_subset, reduction="pca", dims = 1:num_dims)

## Visualization ##
metadata_fields <- c("region_label", "cluster_label", "cluster_1", "cluster_2", "seurat_clusters")
reduction_types <- c("pca", "umap", "tsne")
plot_directory <- "/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/02.Cluster/"

reduc <- "umap"
field <- "cluster_2"
plot_path <- paste0(plot_directory, "GSE185862_", toupper(reduc), "plot_", field, ".png")

CairoPNG(plot_path, width = 800, height = 600)
DimPlot(hip_subset, reduction = reduc, group.by = field, raster = FALSE, label = T)
dev.off()

#### Cell annotation ####
## Data-derived Marker genes ##
cluster2.markers <- FindMarkers(hip_subset, ident.1 = 2, min.pct = 0.25)
hip_subset.markers <- FindAllMarkers(hip_subset, only.pos=TRUE)
hip_subset.markers <- FindAllMarkers(hip_subset, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25, test.use="MAST") ; 
head(hip_subset.markers, 30)
# pct.1 -> the proportion of cells expressing the marker gene in cluster 2
# pct.2 -> the proportion of cells expression the marker gene in the rest of clusters
hip_subset.markers.top20 <- hip_subset.markers %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(n = 20, wt=avg_log2FC) # filtering top genes using log2FC
hip_subset.markers.top10 = hip_subset.markers %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(n = 10, wt=avg_log2FC)
hip_subset.markers.top5 = hip_subset.markers %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(n = 5, wt=avg_log2FC)
#DoHeatmap(hip_subset, features = hip_subset.markers.top20$gene)
#DoHeatmap(hip_subset, features = hip_subset.markers.top10$gene)
#DoHeatmap(hip_subset, features = hip_subset.markers.top5$gene)

# dataframe of top20 genes for each clusters (hip_subset.markers.top20s)
mySeuratClusters = unique(hip_subset.markers.top20$cluster)
for(c in 1:length(mySeuratClusters)){
  hip_subset.markers.top20.c <- data.frame(
    cluster <- hip_subset.markers.top20[hip_subset.markers.top20$cluster %in% mySeuratClusters[c], "gene"]) ;
  colnames(hip_subset.markers.top20.c) <- mySeuratClusters[c]
  if(c == 1){hip_subset.markers.top20s <- hip_subset.markers.top20.c} else {
    hip_subset.markers.top20s = cbind(hip_subset.markers.top20s, hip_subset.markers.top20.c)}
}

## Previously known Marker genes ##
ExN.markers = c("Slc17a7", "Prox1", "Cfap299", "Syn3", "Hgf", "Grik1")
CairoPNG("/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/03.Markers/GSE185862_UMAPplot_ExN.png", 
         width = 800, height = 800)
FeaturePlot(hip_subset, features=ExN.markers, reduction="umap", raster=FALSE) ;
dev.off()

InN.markers = c("Gad1", "Pvalb", "Sst", "Vip")
CairoPNG("/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/03.Markers/GSE185862_UMAPplot_InN.png", 
         width = 800, height = 800)
FeaturePlot(hip_subset, features=InN.markers, reduction="umap", raster=FALSE) ;
dev.off()

markers = c("Slc17a7", "Gad1", "Ptprc", "Gfap")
CairoPNG("/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/03.Markers/GSE185862_UMAPplot_Whole.png", 
         width = 800, height = 600)
FeaturePlot(hip_subset, features=markers, reduction="umap", raster=FALSE) ;
dev.off()

## Assigning cell type identity to clusters ##


############### Save ############### 
saveRDS(hip_subset, file = "/data/project/HS/scRNA_Animal/GSE185862/GSE185862_02SeuratDeconvSig_object.rds")
hip_subset <- readRDS(file = "/data/project/HS/scRNA_Animal/GSE185862/GSE185862_02SeuratDeconvSig_object.rds")


####################################
#### (ss.seurat) visualization for all (including CTX) ###
## FindVariableFeatures
ss.seurat <- FindVariableFeatures(ss.seurat, nfeatures = 2000)
## PCA
ss.seurat <- RunPCA(ss.seurat, features = VariableFeatures(hip_subset))
num_dims <- 22
## UMAP
ss.seurat <- RunUMAP(ss.seurat, reduction="pca", dims = 1:num_dims)
## TSNE 
ss.seurat <- RunTSNE(ss.seurat, reduction="pca", dims = 1:num_dims)
## modified region label
ss.seurat@meta.data$region_label[is.na(ss.seurat@meta.data$region_label)] <- "Others"
ss.seurat@meta.data$region_label_modified <- ifelse(ss.seurat@meta.data$region_label == "HIP", "HIP", "Others")
## Visualization
metadata_fields <- c("region_label", "region_label_modified", "cluster_label")
reduction_types <- c("pca", "umap", "tsne")
plot_directory <- "/home/joonho345/3_RNA/scRNA_Animal/02.GSE185862_Seurat_Deconv/02.Cluster/"
reduc <- "umap"
field <- "region_label"
plot_path <- paste0(plot_directory, "ALL_GSE185862_", toupper(reduc), "plot_", field, ".png")
CairoPNG(plot_path, width = 800, height = 600)
DimPlot(ss.seurat, reduction = reduc, group.by = field, raster = FALSE, label = T)
dev.off()


####################################
############### (hip_subset_1) Input for Cell Deconvolution ############### 
## XX - Get the cell types of interest 
## XX - for MuSiC ALL genes

#### for CIBERSORT with 1000 genes ####
## 02.Seurat_Deconv_M_Sig.R 
