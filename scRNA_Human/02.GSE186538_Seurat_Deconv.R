library(Seurat)
library(Matrix)

setwd('/home/joonho345/3_RNA/scRNA_Human/')
GSE <- 'GSE186538'

#### Create Seurat Object ####
metadata_path <- paste0("Raw_Data/", GSE, "/", GSE, "_Human_cell_meta.txt")
genes_path <- paste0("Raw_Data/", GSE, "/", GSE, "_Human_genes.txt")
counts_path <- paste0("Raw_Data/", GSE, "/", GSE, "_Human_counts.mtx")

metadata <- read.table(metadata_path, header = TRUE, sep = "\t", row.names = 1)
genes <- read.table(genes_path, header = FALSE, stringsAsFactors = FALSE)
colnames(genes) <- c("Gene")
counts <- readMM(counts_path)

rownames(counts) <- genes$Gene
colnames(counts) <- rownames(metadata)
human_hippo <- CreateSeuratObject(counts = counts, meta.data = metadata)
human_hippo[["percent.mt"]] <- PercentageFeatureSet(human_hippo, pattern = "^MT-")

human_hippo
head(human_hippo@meta.data)
# Cell: 219058 / Gene: 33939

#### Quality Control ####
## visualization
CairoPNG('/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/01.QC/GSE186538_VlnPlot.png', width = 800, height = 600)
VlnPlot(human_hippo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

CairoPNG('/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/01.QC/GSE186538_FeatureScatter.png', width = 800, height = 600)
plot1 <- FeatureScatter(human_hippo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(human_hippo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

library(SeuratWrappers)
miqc.obj <- RunMiQC(human_hippo, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", model.slot = "flexmix_model")
head(miqc.obj@meta.data)
miqc.obj@meta.data$miQC.keep
miqc.obj@meta.data$miQC.probability
CairoPNG('/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/01.QC/GSE186538_PlotMiQC.png', width = 800, height = 600)
PlotMiQC(miqc.obj, color.by = "miQC.keep")
dev.off()
miqc.obj@meta.data %>% 
  group_by(miQC.keep) %>% 
  summarise(min_nFeature_RNA = min(nFeature_RNA), max_nFeature_RNA = max(nFeature_RNA),
            min_percent.mt = min(percent.mt), max_percent.mt = max(percent.mt))

## cell QC
human_hippo <- subset(human_hippo, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
human_hippo # Cell: 200481 / Gene: 33939
head(human_hippo@meta.data)
# Cell: 200481 / Gene: 33939

## feature QC
min_cells <- 0.01 * ncol(human_hippo)  # Example: 1% of the total cells
counts_matrix <- human_hippo@assays$RNA@layers$counts
genes_to_keep <- rowSums(counts_matrix > 0) > min_cells
length(genes_to_keep) # Gene: 33939

#### Normalization ####
human_hippo <- NormalizeData(human_hippo, scale.factor = 10000)
human_hippo@assays$RNA@layers$counts[1:100,1:100]
human_hippo@assays$RNA@layers$data[1:100,1:100]
human_hippo <- ScaleData(human_hippo, features = rownames(human_hippo))
human_hippo@assays$RNA@layers$scale.data[1:100,1:100]

#### Cluster_1 & Cluster_2 ####
exn_cells <- c(
  "CA3 CFAP299 SYN3", "CA1 dorsal GRIK1 GRM3", "DG MC ARHGAP24 DLC1", "EC L6 TLE4 SULF1",
  "CA1 ventral ACVR1C SYT13", "EC L3 PCP4 ADARB2", "EC L2 CUX2 CALB1", "SUB proximal ROBO1 COL5A2",
  "EC L2 CUX2 IL1RAPL2", "CA2 CFAP299 HGF", "EC L2 RELN BCL11B", "DG GC PROX1 SGCZ",
  "DG GC PROX1 PDLIM5", "SUB proximal ROBO1 SEMA3E", "SUB distal FN1 NTNG1", "EC L2 CUX2 PDGFD",
  "EC L5 RORB TLL1", "EC L5 RORB TPBG", "EC L56 TLE4 NXPH2", "EC L6 THEMIS CDH13",
  "EC L5 BCL11B ADRA1A", "EC L2 CUX2 LAMA3", "EC L6b TLE4 CCN2", "EC L6 THEMIS RGS12",
  "EC L2 RELN BMPR1B")
exn_dg_cells <- c("DG MC ARHGAP24 DLC1", "DG GC PROX1 SGCZ", "DG GC PROX1 PDLIM5")
exn_ca_cells <- c("CA3 CFAP299 SYN3", "CA1 dorsal GRIK1 GRM3", "CA1 ventral ACVR1C SYT13", "CA2 CFAP299 HGF")
exn_sub_cells <- c("SUB proximal ROBO1 COL5A2", "SUB proximal ROBO1 SEMA3E", "SUB distal FN1 NTNG1")
exn_ec_cells <- c(
  "EC L6 TLE4 SULF1", "EC L3 PCP4 ADARB2", "EC L2 CUX2 CALB1", "EC L2 CUX2 IL1RAPL2",
  "EC L2 RELN BCL11B", "EC L2 CUX2 PDGFD", "EC L5 RORB TLL1", "EC L5 RORB TPBG",
  "EC L56 TLE4 NXPH2", "EC L6 THEMIS CDH13", "EC L5 BCL11B ADRA1A", "EC L2 CUX2 LAMA3",
  "EC L6b TLE4 CCN2", "EC L6 THEMIS RGS12", "EC L2 RELN BMPR1B")
inn_cells <- c(
  "InN LAMP5 KIT", "InN VIP SCTR", "InN SST EPB41L4A", "InN VIP ABI3BP", "InN PVALB PLCL1",
  "InN SST ADAMTS12", "InN LAMP5 NMBR", "InN NR2F2 DDR2", "InN PVALB PLEKHH2", "InN SST OTOF",
  "InN LAMP5 CHST9", "InN VIP NOX4", "InN VIP SCML4", "InN VIP CHRNA2", "InN NR2F2 SLC17A8",
  "InN PVALB MEPE", "InN NR2F2 ANO2", "InN NR2F2 PTPRK", "InN LHX6 AC008415.1", "InN VIP PENK",
  "InN NR2F2 MIR4300HG", "InN SST NPY", "InN MEIS2 SHISAL2B")
cr_cells <- c("CR RELN NDNF")
oligo_cells <- c(
  "Oligo OPALIN LAMA2", "Oligo CPXM2 KANK4", "Oligo OPALIN SLC5A11", "Oligo OPALIN LINC01098",
  "OPC PDGFRA GRIA4", "OPC PDGFRA EGR1", "COP GPR17 ADAM33")
vascular_cells <- c(
  "aSMC ACTA2 TAGLN", "vSMC ABCC9 P2RY14", "PC CLDN5 ABCC9", "Endo CLDN5 VWF", 
  "VLMC COL1A1 COL1A2", "aEndo DKK2 FBLN5")
astro_cells <- c("Astro AQP4 GFAP", "Astro AQP4 CHRDL1")
micro_cells <- c("Micro C1QB P2RY12", "Micro C1QB CD83")
imm_cells <- c("Macro F13A1 COLEC12", "T SKAP1 CD247", "Myeloid LSP1 LYZ")

## cluster_1 ##
human_hippo$cluster_1 <- NA  # Initialize with NA
human_hippo@meta.data[human_hippo@meta.data$cluster %in% exn_cells, ]$cluster_1 <- "ExN"
human_hippo@meta.data[human_hippo@meta.data$cluster %in% inn_cells, ]$cluster_1 <- "InN"
human_hippo@meta.data[human_hippo@meta.data$cluster %in% cr_cells, ]$cluster_1 <- "CR"
human_hippo@meta.data[human_hippo@meta.data$cluster %in% oligo_cells, ]$cluster_1 <- "Oligo"
human_hippo@meta.data[human_hippo@meta.data$cluster %in% vascular_cells, ]$cluster_1 <- "Vascular"
human_hippo@meta.data[human_hippo@meta.data$cluster %in% astro_cells, ]$cluster_1 <- "Astro"
human_hippo@meta.data[human_hippo@meta.data$cluster %in% micro_cells, ]$cluster_1 <- "Micro"
human_hippo@meta.data[human_hippo@meta.data$cluster %in% imm_cells, ]$cluster_1 <- "Imm"

## cluster_2 ##
human_hippo$cluster_2 <- NA  # Initialize with NA
human_hippo@meta.data[human_hippo@meta.data$cluster %in% exn_dg_cells, ]$cluster_2 <- "ExN_DG"
human_hippo@meta.data[human_hippo@meta.data$cluster %in% exn_ca_cells, ]$cluster_2 <- "ExN_CA"
human_hippo@meta.data[human_hippo@meta.data$cluster %in% exn_sub_cells, ]$cluster_2 <- "ExN_SUB"
human_hippo@meta.data[human_hippo@meta.data$cluster %in% exn_ec_cells, ]$cluster_2 <- "ExN_EC"
human_hippo@meta.data[human_hippo@meta.data$cluster %in% inn_cells, ]$cluster_2 <- "InN"
human_hippo@meta.data[human_hippo@meta.data$cluster %in% cr_cells, ]$cluster_2 <- "CR"
human_hippo@meta.data[human_hippo@meta.data$cluster %in% oligo_cells, ]$cluster_2 <- "Oligo"
human_hippo@meta.data[human_hippo@meta.data$cluster %in% vascular_cells, ]$cluster_2 <- "Vascular"
human_hippo@meta.data[human_hippo@meta.data$cluster %in% astro_cells, ]$cluster_2 <- "Astro"
human_hippo@meta.data[human_hippo@meta.data$cluster %in% micro_cells, ]$cluster_2 <- "Micro"
human_hippo@meta.data[human_hippo@meta.data$cluster %in% imm_cells, ]$cluster_2 <- "Imm"

## filter cells ##
human_hippo <- subset(human_hippo, subset = cluster_1 != "Imm" & cluster_1 != "Vascular" & cluster_1 != "CR")
human_hippo # 33939 features across 193808 samples
unique(human_hippo@meta.data$cluster_1)
unique(human_hippo@meta.data$cluster_2)
sum(is.na(human_hippo@meta.data$cluster_1)) # no NA
sum(is.na(human_hippo@meta.data$cluster_1)) # no NA
# cluster_1: ExN  InN Astro Oligo Micro
# cluster_1: ExN_DG ExN_CA  ExN_SUB ExN_EC  InN Astro Oligo Micro


############### Clustering & Cell annotation ############### 
#### FindVariableFeatures ####
human_hippo <- FindVariableFeatures(human_hippo, nfeatures = 2000)
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/01.QC/GSE186538_VariableFeaturePlot.png", 
         width = 800, height = 600)
VariableFeaturePlot(human_hippo)
dev.off()

#### Clustering ####
## PCA ##
human_hippo <- RunPCA(human_hippo, features = VariableFeatures(human_hippo))
# Decide the number of PC dimensions
# ElbowPlot
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/02.Cluster/GSE186538_ElbowPlot.png", 
         width = 800, height = 600)
ElbowPlot(human_hippo, ndims = 50)
dev.off()
# JackStrawPlot
#human_hippo <- JackStraw(human_hippo, dims = 50)
#human_hippo <- ScoreJackStraw(human_hippo, dims = 1:50)
#CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/GSE186538_02Cluster_JackStrawPlot.png", 
#         width = 800, height = 600)
#JackStrawPlot(human_hippo, dims = 1:50)
#dev.off()
# Set the number of dimensions
num_dims <- 22
message("Number of dimensions chosen: ", num_dims)

## Clustering ##
human_hippo <- FindNeighbors(human_hippo, dims = 1:num_dims)
human_hippo@graphs$RNA_nn@x
human_hippo@graphs$RNA_snn@x
human_hippo <- FindClusters(human_hippo, resolution = 0.4) # 0.4 ~ 1.2
human_hippo@meta.data$seurat_clusters

## UMAP ##
human_hippo <- RunUMAP(human_hippo, reduction="pca", dims = 1:num_dims)

## TSNE ##
human_hippo <- RunTSNE(human_hippo, reduction="pca", dims = 1:num_dims)

## Visualization ##
metadata_fields <- c("orig.ident", "cluster", "cluster_1", "cluster_2", "seurat_clusters")
reduction_types <- c("pca", "umap", "tsne")
plot_directory <- "/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/02.Cluster/"

reduc <- "umap"
field <- "cluster"
plot_path <- paste0(plot_directory, "GSE186538_", toupper(reduc), "plot_", field, ".png")

CairoPNG(plot_path, width = 1200, height = 600)
DimPlot(human_hippo, reduction = reduc, group.by = field, raster = FALSE, label = T)
dev.off()

#### Cell annotation ####
## Data-derived Marker genes ##
cluster2.markers <- FindMarkers(human_hippo, ident.1 = 2, min.pct = 0.25)
human_hippo.markers <- FindAllMarkers(human_hippo, only.pos=TRUE)
human_hippo.markers <- FindAllMarkers(human_hippo, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25, test.use="MAST") ; 
head(human_hippo.markers, 30)
# pct.1 -> the proportion of cells expressing the marker gene in cluster 2
# pct.2 -> the proportion of cells expression the marker gene in the rest of clusters
human_hippo.markers.top20 <- human_hippo.markers %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(n = 20, wt=avg_log2FC) # filtering top genes using log2FC
human_hippo.markers.top10 = human_hippo.markers %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(n = 10, wt=avg_log2FC)
human_hippo.markers.top5 = human_hippo.markers %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(n = 5, wt=avg_log2FC)
#DoHeatmap(human_hippo, features = human_hippo.markers.top20$gene)
#DoHeatmap(human_hippo, features = human_hippo.markers.top10$gene)
#DoHeatmap(human_hippo, features = human_hippo.markers.top5$gene)

# dataframe of top20 genes for each clusters (human_hippo.markers.top20s)
mySeuratClusters = unique(human_hippo.markers.top20$cluster)
for(c in 1:length(mySeuratClusters)){
  human_hippo.markers.top20.c <- data.frame(
    cluster <- human_hippo.markers.top20[human_hippo.markers.top20$cluster %in% mySeuratClusters[c], "gene"]) ;
  colnames(human_hippo.markers.top20.c) <- mySeuratClusters[c]
  if(c == 1){human_hippo.markers.top20s <- human_hippo.markers.top20.c} else {
    human_hippo.markers.top20s = cbind(human_hippo.markers.top20s, human_hippo.markers.top20.c)}
}

## Previously known Marker genes ##
ExN.markers = c("SLC17A7", "PROX1", "CFAP299", "SYN3", "HGF", "GRIK1")
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/03.Markers/GSE186538_UMAPplot_ExN.png", 
         width = 800, height = 800)
FeaturePlot(human_hippo, features=ExN.markers, reduction="umap", raster=FALSE) ;
dev.off()

InN.markers = c("GAD1", "PVALB", "SST", "VIP")
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/03.Markers/GSE186538_UMAPplot_InN.png", 
         width = 800, height = 800)
FeaturePlot(human_hippo, features=InN.markers, reduction="umap", raster=FALSE) ;
dev.off()

Astro.markers = c("AQP4", "GFAP", "AQP1", "APOE")
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/03.Markers/GSE186538_UMAPplot_Astro.png", 
         width = 800, height = 800)
FeaturePlot(CA, features=Astro.markers, reduction="umap", raster=FALSE) ;
dev.off()

Oligo.markers = c("MOBP", "SOX10", "OPALIN", "PDGFRA", "GPR17")
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/03.Markers/GSE186538_UMAPplot_Oligo.png", 
         width = 800, height = 800)
FeaturePlot(CA, features=Oligo.markers, reduction="umap", raster=FALSE) ;
dev.off()

Micro.markers = c("AIF1", "PTPRC", "SLC1A3")
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/03.Markers/GSE186538_UMAPplot_Micro.png", 
         width = 800, height = 800)
FeaturePlot(CA, features=Micro.markers, reduction="umap", raster=FALSE) ;
dev.off()

markers = c("SLC17A7", "GAD1", "PTPRC", "GFAP")
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.GSE186538_Seurat_Deconv/03.Markers/GSE186538_UMAPplot_Whole.png", 
         width = 800, height = 600)
FeaturePlot(human_hippo, features=markers, reduction="umap", raster=FALSE) ;
dev.off()


## Assigning cell type identity to clusters ##


############### Save ############### 
saveRDS(human_hippo, file = "/data/project/HS/scRNA_Human/GSE186538/GSE186538_02SeuratDeconvSig_object.rds")
human_hippo <- readRDS(file = "/data/project/HS/scRNA_Human/GSE186538/GSE186538_02SeuratDeconvSig_object.rds")
