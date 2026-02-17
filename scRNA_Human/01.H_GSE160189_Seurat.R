library(Seurat)
library(Matrix)

setwd('/home/joonho345/1_Epilepsy_RNA/scRNA_Human/')
GSE <- 'TLE_GSE160189'

#########################################################
#### Create Seurat Object ####
raw_counts <- read.csv("/data/project/1_Epilepsy_RNA/scRNA_Human/Raw_Data/TLE_GSE160189/GSE160189_Hippo_Counts.csv", row.names = 1)
metadata <- read.csv("/data/project/1_Epilepsy_RNA/scRNA_Human/Raw_Data/TLE_GSE160189/HippoeSeq_Fig1_All_Meta.txt", sep = "\t", row.names = 1)
GSE160189_hippo <- CreateSeuratObject(counts = raw_counts, meta.data = metadata, project = GSE)

print(GSE160189_hippo) # 17180 features across 131325 samples
head(GSE160189_hippo@meta.data)
table(GSE160189_hippo@meta.data$orig.ident)
table(GSE160189_hippo@meta.data$group)
table(GSE160189_hippo@meta.data$batch)
table(GSE160189_hippo@meta.data$donor)


#########################################################
#### Quality Control ####
## cell QC (already done by 2021 Neuron)
# nCount_RNA< 10000 & nFeature_RNA >300 & percent.mt<5

## feature QC (already done by 2021 Neuron)
# remove mGenes,xGenes,yGenes


#########################################################
#### Cell Annotation ####
table(GSE160189_hippo@meta.data$CellType)
Den.Gyr1 Den.Gyr2 Den.Gyr3
Pyr1 Pyr2
In1 In2 In3
Astro1 Astro2 Astro3
Micro1 Micro2 Micro3
Olig1 Olig2 Olig3 Olig4 Olig5 
OPC1 OPC2 OPC3 OPC4
Endo
NA


##########################################################
#### Trimming meta.data for integration object ####
GSE160189_hippo_int <- GSE160189_hippo
GSE160189_hippo_int@meta.data$orig.ident <- GSE160189_hippo_int@meta.data$donor
GSE160189_hippo_int@meta.data$species <- 'human'
columns_to_remove <- c("percent.mt", "group", "age", "sex", "dur", "batch", "version", "donor", "UMAP_1", "UMAP_2")
GSE160189_hippo_int@meta.data <- GSE160189_hippo_int@meta.data[, !colnames(GSE160189_hippo_int@meta.data) %in% columns_to_remove]

head(GSE160189_hippo_int@meta.data)
table(GSE160189_hippo_int@meta.data$orig.ident)
table(GSE160189_hippo_int@meta.data$CellType)

