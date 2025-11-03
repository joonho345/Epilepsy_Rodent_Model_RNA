library(Seurat)
library(Matrix)
library(Cairo)

setwd('/home/joonho345/3_RNA/scRNA_Animal/')
GSE <- 'GSE185862'

#########################################
####### 
#######
#########################################

#### Create Seurat Object ####
file_path <- "/home/joonho345/3_RNA/scRNA_Animal/Raw_Data/Animal_GSE185862/Seurat.ss.rda"
load(file_path)
ss.seurat[["percent.mt"]] <- PercentageFeatureSet(ss.seurat, pattern = "^mt-")
print(ss.seurat)
# Cell: 73363 / Gene: 42055

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

#### (hip_subset) cells from hippocampus (not CTX) ####
unique(ss.seurat@meta.data$region_label)
hip_subset <- subset(ss.seurat, subset = region_label == "HIP")
unique(hip_subset@meta.data$region_label)
print(hip_subset)
# Cell: 6564 / Gene: 42055


##########################################################
#### Trimming meta.data ... leave orig.ident & CellType ####
GSE185862_hippo_M_int <- ss.seurat
GSE185862_hippo_M_int@meta.data$orig.ident <- GSE185862_hippo_M_int@meta.data$donor_sex_id
GSE185862_hippo_M_int@meta.data$CellType <- GSE185862_hippo_M_int@meta.data$subclass_label
GSE185862_hippo_M_int@meta.data$species <- 'mouse'
columns_to_remove <- c(
  "sample_name", "exp_component_name", "platform_label",
  "cluster_color", "cluster_id", "cluster_label", "class_color", "class_id", "class_label",
  "subclass_color", "subclass_id", "subclass_label", "donor_sex_color", "donor_sex_id", "donor_sex_label",
  "region_color", "region_id", "region_label", "cell_set_accession_color", "cell_set_accession_id",
  "cell_set_accession_label", "cell_set_alt_alias_color", "cell_set_alt_alias_id", "cell_set_alt_alias_label",
  "cell_set_designation_color", "cell_set_designation_id", "cell_set_designation_label",
  "neighborhood_label", "neighborhood_id", "neighborhood_color", "percent.mt"
)
GSE185862_hippo_M_int@meta.data <- GSE185862_hippo_M_int@meta.data[, !colnames(GSE185862_hippo_M_int@meta.data) %in% columns_to_remove]

head(GSE185862_hippo_M_int@meta.data)
table(GSE185862_hippo_M_int@meta.data$orig.ident)
table(GSE185862_hippo_M_int@meta.data$CellType)
table(GSE185862_hippo_M_int@meta.data$neighborhood_id)


##########################################################
ss.seurat <- NormalizeData(ss.seurat, scale.factor = 10000)
ss.seurat <- ScaleData(ss.seurat, features = rownames(ss.seurat))
ss.seurat <- FindVariableFeatures(ss.seurat, nfeatures = 2000)
ss.seurat <- RunPCA(ss.seurat, features = VariableFeatures(ss.seurat))
ss.seurat <- RunUMAP(ss.seurat, reduction="pca", dims = 1:30)

CairoPNG('/home/joonho345/3_RNA/scRNA_Animal/03.GSE185862_Seurat/01.02umap_donor_sex_id.png', width = 800, height = 600)
DimPlot(ss.seurat, reduction = "umap", group.by = "donor_sex_id", shuffle = TRUE, label = TRUE)
dev.off()
