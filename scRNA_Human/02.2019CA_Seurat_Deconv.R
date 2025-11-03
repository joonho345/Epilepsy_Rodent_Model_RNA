####### scRNAseq #######
####### NatCom: CA for Human #######
####### NatCom: TS for Mouse #######
# https://github.com/Voineagulab/BrainCellularComposition/blob/main/Preprocessing_a_SingleNucleusData.R
# https://github.com/Voineagulab/BrainCellularComposition/blob/main/Preprocessing_b_SignatureGeneration.R

library(Seurat)
library(Matrix)

setwd('/home/joonho345/3_RNA/scRNA_Human/')
GSE <- '2019CA'

#### Functions ####
# for seurat preprocessing
min.cells <- 3 # during the initial load, a gene is excluded if in < 3 cells 
min.features <- 200 # during the initial load, a barcode is excluded < 200 features are expressed
min.depth <- 1000 # a barcode is excluded if nCount_RNA < this value
max.depth.percentile <- 0.995 # a barcode is excluded if nCount_RNA > this percentile within the dataset
max.mito <- 5
# preprocessing options
downsample <- FALSE
downsample.n <- NA; if (downsample) downsample.n <- NA
use.SCTransform <- FALSE
# Function for downsampling the dataset to a set number of barcodes
downsample.fun <- function(x, n = downsample.n) {
  if (ncol(x) <= downsample.n) {
    print("No downsampling performed (Reason: number of cells in dataset is already less than or equal to the downsampling number)")
  } else {
    sample <- sample(colnames(x), size = n, replace = FALSE)
    x <- subset(x, cells = sample)  
  }
  return(x)
}
# General function for preprocessing sn data (normalise, filters, and scales)
get.max.depth <- function(x) {
  max.depth <- quantile(x@meta.data$nCount_RNA, probs = max.depth.percentile)
}
preprocess.fun <- function(x, run.downsample = downsample, SCTransform = use.SCTransform, max.depth = max.depth) {
  x[["percent.mito"]] <- PercentageFeatureSet(object = x, pattern = "^MT-")
  x <- subset(x = x, subset = (nCount_RNA > min.depth) & (nCount_RNA < max.depth) & (percent.mito < max.mito))
  if (run.downsample) { x <- downsample.fun(x) }
  x <- NormalizeData(object = x, normalization.method = "LogNormalize", scale.factor = 10000) # standard parameters for Seurat
  x <- FindVariableFeatures(object = x, selection.method = "vst", nfeatures = 2000)
  if (use.SCTransform) {
    x <- SCTransform(object = x, vars.to.regress = c("nCount_RNA", "percent.mito")) 
  }
  return(x)
}


#### Create Seurat Object (NatCom: CA for Human) ####
# Download this and unzip: https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-smart-seq  
dat <- read.csv("/home/joonho345/3_RNA/scRNA_Human/Raw_Data/2019CA/NatCom/CA/human_MTG_2018-06-14_exon-matrix.csv") 
meta <- read.csv("/home/joonho345/3_RNA/scRNA_Human/Raw_Data/2019CA/NatCom/CA/human_MTG_2018-06-14_genes-rows.csv")
dat <- dat[,-1] # remove an annotation column
rownames(dat) <- meta$gene
# create Seurat object
CA <- CreateSeuratObject(counts = dat,
                         min.cells = round(ncol(dat) / 100),
                         min.features = min.features,
                         project = "HCA")
# augment metadata
meta <- read.csv("/home/joonho345/3_RNA/scRNA_Human/Raw_Data/2019CA/NatCom/CA/human_MTG_2018-06-14_samples-columns.csv")
rownames(meta) <- meta$sample_name
meta <- meta[colnames(CA),]
CA$Individual <- meta$donor # samples of CA
CA$orig.celltype <- meta$cluster # cell type of CA
# Remove cells with no class
keep <- which(!(CA$orig.celltype == "no class"))
CA <- subset(CA, cells = keep)
# Preprocess
max.depth <- get.max.depth(CA)
CA <- preprocess.fun(CA, max.depth = max.depth)
# types
unique(CA@meta.data$orig.ident) # Levels: F1S4 F2S4
unique(CA@meta.data$orig.celltype)


#### Quality Control ####
## visualization
CairoPNG('/home/joonho345/3_RNA/scRNA_Human/02.2019CA_Seurat_Deconv/01.QC/2019CA_VlnPlot.png', width = 800, height = 600)
VlnPlot(CA, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
dev.off()

CairoPNG('/home/joonho345/3_RNA/scRNA_Human/02.2019CA_Seurat_Deconv/01.QC/2019CA_FeatureScatter.png', width = 800, height = 600)
FeatureScatter(CA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

## cell QC - Done
## feature QC - Done
## Normalization (NormalizeData - Done)
CA <- ScaleData(CA, features = rownames(CA))


#### Clustering ####
## filter cells ##
CA <- subset(CA, cells = grep("Endo", CA$orig.celltype, invert = TRUE))
CA # 27903 features across 15515 samples

## cluster_1 ##
rename.ct <- function(x, orig.label, new.label) {
  if (!("cluster_1" %in% colnames(x@meta.data))) x$cluster_1 <- "."
  g <- grep(orig.label, x$orig.celltype)
  x$cluster_1[g] <- new.label
  return(x)
}
CA <- rename.ct(CA, "Astro", "Astro")
CA <- rename.ct(CA, "Inh", "InN")
CA <- rename.ct(CA, "Exc", "ExN")
CA <- rename.ct(CA, "Micro", "Micro")
CA <- rename.ct(CA, "Oligo", "Oligo")
CA <- rename.ct(CA, "OPC", "Oligo")
unique(CA@meta.data$cluster_1)
sum(is.na(CA@meta.data$cluster_1)) # no NA
# cluster_1: ExN  InN Astro Oligo Micro

############### Clustering & Cell annotation ############### 
#### FindVariableFeatures ####
CA <- FindVariableFeatures(CA, nfeatures = 2000)
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.2019CA_Seurat_Deconv/01.QC/2019CA_VariableFeaturePlot.png", 
         width = 800, height = 600)
VariableFeaturePlot(CA)
dev.off()

#### Clustering ####
## PCA ##
CA <- RunPCA(CA, features = VariableFeatures(CA))
# Decide the number of PC dimensions
# ElbowPlot
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.2019CA_Seurat_Deconv/02.Cluster/2019CA_ElbowPlot.png", 
         width = 800, height = 600)
ElbowPlot(CA, ndims = 50)
dev.off()
# JackStrawPlot
#CA <- JackStraw(CA, dims = 50)
#CA <- ScoreJackStraw(CA, dims = 1:50)
#CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.2019CA_Seurat_Deconv/02.Cluster/2019CA_JackStrawPlot.png", 
#         width = 800, height = 600)
#JackStrawPlot(CA, dims = 1:50)
#dev.off()
# Set the number of dimensions
num_dims <- 30
message("Number of dimensions chosen: ", num_dims)

## UMAP ##
CA <- RunUMAP(CA, reduction="pca", dims = 1:num_dims)

## TSNE ##
CA <- RunTSNE(CA, reduction="pca", dims = 1:num_dims)


#### Visualization of clusters ####
# Clusters
metadata_fields <- c("orig.ident", "Individual", "orig.celltype", "cluster_1")
reduction_types <- c("pca", "umap", "tsne")
plot_directory <- "/home/joonho345/3_RNA/scRNA_Human/02.2019CA_Seurat_Deconv/02.Cluster/"

reduc <- "umap"
field <- "cluster_1"
plot_path <- paste0(plot_directory, "2019CA_", toupper(reduc), "plot_", field, ".png")

CairoPNG(plot_path, width = 800, height = 600)
DimPlot(CA, reduction = reduc, group.by = field, raster = FALSE, label = T)
dev.off()


# Marker genes
ExN.markers = c("SLC17A7", "PROX1", "CFAP299", "SYN3", "HGF", "GRIK1")
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.2019CA_Seurat_Deconv/03.Markers/2019CA_UMAPplot_ExN.png", 
         width = 800, height = 800)
FeaturePlot(CA, features=ExN.markers, reduction="umap", raster=FALSE) ;
dev.off()

InN.markers = c("GAD1", "PVALB", "SST", "VIP")
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.2019CA_Seurat_Deconv/03.Markers/2019CA_UMAPplot_InN.png", 
         width = 800, height = 800)
FeaturePlot(CA, features=InN.markers, reduction="umap", raster=FALSE) ;
dev.off()

Astro.markers = c("AQP4", "GFAP", "AQP1", "APOE")
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.2019CA_Seurat_Deconv/03.Markers/2019CA_UMAPplot_Astro.png", 
         width = 800, height = 800)
FeaturePlot(CA, features=Astro.markers, reduction="umap", raster=FALSE) ;
dev.off()

Oligo.markers = c("MOBP", "SOX10", "OPALIN", "PDGFRA", "GPR17")
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.2019CA_Seurat_Deconv/03.Markers/2019CA_UMAPplot_Oligo.png", 
         width = 800, height = 800)
FeaturePlot(CA, features=Oligo.markers, reduction="umap", raster=FALSE) ;
dev.off()

Micro.markers = c("AIF1", "PTPRC", "SLC1A3")
CairoPNG("/home/joonho345/3_RNA/scRNA_Human/02.2019CA_Seurat_Deconv/03.Markers/2019CA_UMAPplot_Micro.png", 
         width = 800, height = 800)
FeaturePlot(CA, features=Micro.markers, reduction="umap", raster=FALSE) ;
dev.off()




############################################################
### NatCom: TS for Mouse ###
dat <- read.csv("/home/joonho345/3_RNA/RNA_Human/script/05.Deconvolution/NatCom/TS/GSE115746_cells_exon_counts.csv")
# meta <- read.csv("Raw/Tasic2018_ST10_Full_Metadata.csv")
meta <- readxl::read_xlsx("/home/joonho345/3_RNA/RNA_Human/script/05.Deconvolution/NatCom/TS/Supplementary_Table_10_Full_Metadata.xlsx", sheet = 1)
meta <- as.data.frame(meta)
rownames(meta) <- meta$sample_name
## Convert to human EnsID
# move gene names to rows
rownames(dat) <- dat$X
# filter to homologous genes
homologues <- read.delim("/home/joonho345/3_RNA/RNA_Human/script/05.Deconvolution/NatCom/TS//HOM_MouseHumanSequence.txt")
dat <- dat[which(rownames(dat) %in% homologues$Symbol),]
# replace mouse gene symbol with human gene symbol
m <- match(rownames(dat), homologues$Symbol)
m <- m + 1
dat$Symbol <- homologues$Symbol[m]
# removes NAs
remove <- which(is.na(dat$Symbol)) 
dat <- dat[-remove,]
# removes duplicated genes
remove <- dat$Symbol[which(duplicated(dat$Symbol))] 
dat <- dat[-which(dat$Symbol %in% remove),]
rownames(dat) <- dat$Symbol
dat <- dat[,-which(colnames(dat) == "Symbol")]
## First pass filtering of cells
# restrict to cells with metadata
keep <- intersect(colnames(dat), meta$sample_name)
dat <- dat[,keep]
meta <- meta[keep,]
# restrict to ALM (anterior lateral motor cortex), as that's a frontal cortical areas
keep <- which(meta$brain_region == "ALM")
dat <- dat[,keep]
meta <- meta[keep,]
# filter out classless cells or those classified as low quality
remove <- which(meta$class %in% c("Low Quality", "No Class"))
dat <- dat[,-remove]
meta <- meta[-remove,]
## Make into Seurat object!
TS <- CreateSeuratObject(counts = dat,
                         min.cells = round(ncol(dat) / 100),
                         min.features = min.features,
                         project = "Tasic")
## Augment the metadata
TS$Individual <- as.character(meta$donor)
TS$orig.celltype <- as.character(meta$cluster)
TS$Class <- as.character(meta$class)
TS$Subclass <- as.character(meta$subclass)
## Remove celltypes
keep <- grep("Peri|^SMC|^VLMC", TS$orig.celltype, invert = TRUE) # this removes pericytes, smooth muscle cells, and VLMCs
TS <- subset(TS, cells = keep)  
## Preprocess and normalise
max.depth <- get.max.depth(TS)
TS <- preprocess.fun(TS, max.depth = max.depth)
# Seurat object of TS
TS
##############################


