library(Seurat)
library(Matrix)

setwd('/home/joonho345/1_Epilepsy_RNA/scRNA_Human/')
GSE <- 'GSE186538'

#########################################################
#### Create Seurat Object ####
metadata_path <- "/data/project/1_Epilepsy_RNA/scRNA_Human/Raw_Data/GSE186538/GSE186538_Human_cell_meta.txt"
genes_path <- "/data/project/1_Epilepsy_RNA/scRNA_Human/Raw_Data/GSE186538/GSE186538_Human_genes.txt"
counts_path <- "/data/project/1_Epilepsy_RNA/scRNA_Human/Raw_Data/GSE186538/GSE186538_Human_counts.mtx"

metadata <- read.table(metadata_path, header = TRUE, sep = "\t", row.names = 1)
genes <- read.table(genes_path, header = FALSE, stringsAsFactors = FALSE)
colnames(genes) <- c("Gene")
counts <- readMM(counts_path)

rownames(counts) <- genes$Gene
colnames(counts) <- rownames(metadata)
GSE186538_hippo <- CreateSeuratObject(counts = counts, meta.data = metadata)
GSE186538_hippo[["percent.mt"]] <- PercentageFeatureSet(GSE186538_hippo, pattern = "^MT-")

GSE186538_hippo
head(GSE186538_hippo@meta.data)


#########################################################
#### Quality Control ####
## cell QC
GSE186538_hippo <- subset(GSE186538_hippo, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
GSE186538_hippo # Cell: 200481 / Gene: 33939

## feature QC
min_cells <- 0.01 * ncol(GSE186538_hippo)
counts_matrix <- GSE186538_hippo@assays$RNA@layers$counts
genes_to_keep <- rowSums(counts_matrix > 0) > min_cells
length(genes_to_keep) # Gene: 33939


#########################################################
#### Cell Annotation ####
exn_cells <- c(
  "CA3 CFAP299 SYN3", "CA1 dorsal GRIK1 GRM3", "DG MC ARHGAP24 DLC1", "EC L6 TLE4 SULF1",
  "CA1 ventral ACVR1C SYT13", "EC L3 PCP4 ADARB2", "EC L2 CUX2 CALB1", "SUB proximal ROBO1 COL5A2",
  "EC L2 CUX2 IL1RAPL2", "CA2 CFAP299 HGF", "EC L2 RELN BCL11B", "DG GC PROX1 SGCZ",
  "DG GC PROX1 PDLIM5", "SUB proximal ROBO1 SEMA3E", "SUB distal FN1 NTNG1", "EC L2 CUX2 PDGFD",
  "EC L5 RORB TLL1", "EC L5 RORB TPBG", "EC L56 TLE4 NXPH2", "EC L6 THEMIS CDH13",
  "EC L5 BCL11B ADRA1A", "EC L2 CUX2 LAMA3", "EC L6b TLE4 CCN2", "EC L6 THEMIS RGS12",
  "EC L2 RELN BMPR1B")
inn_cells <- c(
  "InN LAMP5 KIT", "InN VIP SCTR", "InN SST EPB41L4A", "InN VIP ABI3BP", "InN PVALB PLCL1",
  "InN SST ADAMTS12", "InN LAMP5 NMBR", "InN NR2F2 DDR2", "InN PVALB PLEKHH2", "InN SST OTOF",
  "InN LAMP5 CHST9", "InN VIP NOX4", "InN VIP SCML4", "InN VIP CHRNA2", "InN NR2F2 SLC17A8",
  "InN PVALB MEPE", "InN NR2F2 ANO2", "InN NR2F2 PTPRK", "InN LHX6 AC008415.1", "InN VIP PENK",
  "InN NR2F2 MIR4300HG", "InN SST NPY", "InN MEIS2 SHISAL2B")
cr_cells <- c("CR RELN NDNF")
oligo_cells <- c("Oligo OPALIN LAMA2", "Oligo CPXM2 KANK4", "Oligo OPALIN SLC5A11", "Oligo OPALIN LINC01098")
OPC_cells <- c("OPC PDGFRA GRIA4", "OPC PDGFRA EGR1", "COP GPR17 ADAM33")
vascular_cells <- c(
  "aSMC ACTA2 TAGLN", "vSMC ABCC9 P2RY14", "PC CLDN5 ABCC9", "Endo CLDN5 VWF", 
  "VLMC COL1A1 COL1A2", "aEndo DKK2 FBLN5")
astro_cells <- c("Astro AQP4 GFAP", "Astro AQP4 CHRDL1")
micro_cells <- c("Micro C1QB P2RY12", "Micro C1QB CD83")
imm_cells <- c("Macro F13A1 COLEC12", "T SKAP1 CD247", "Myeloid LSP1 LYZ")

## cluster_1 ##
GSE186538_hippo$cluster_1 <- NA  # Initialize with NA
GSE186538_hippo@meta.data[GSE186538_hippo@meta.data$cluster %in% exn_cells, ]$cluster_1 <- "ExN"
GSE186538_hippo@meta.data[GSE186538_hippo@meta.data$cluster %in% inn_cells, ]$cluster_1 <- "InN"
GSE186538_hippo@meta.data[GSE186538_hippo@meta.data$cluster %in% cr_cells, ]$cluster_1 <- "CR"
GSE186538_hippo@meta.data[GSE186538_hippo@meta.data$cluster %in% oligo_cells, ]$cluster_1 <- "Oligo"
GSE186538_hippo@meta.data[GSE186538_hippo@meta.data$cluster %in% vascular_cells, ]$cluster_1 <- "Vascular"
GSE186538_hippo@meta.data[GSE186538_hippo@meta.data$cluster %in% astro_cells, ]$cluster_1 <- "Astro"
GSE186538_hippo@meta.data[GSE186538_hippo@meta.data$cluster %in% micro_cells, ]$cluster_1 <- "Micro"
GSE186538_hippo@meta.data[GSE186538_hippo@meta.data$cluster %in% imm_cells, ]$cluster_1 <- "Imm"
GSE186538_hippo@meta.data[GSE186538_hippo@meta.data$cluster %in% OPC_cells, ]$cluster_1 <- "OPC_H"


#########################################################
#### Trimming meta.data for integration object ####
head(GSE186538_hippo@meta.data)
GSE186538_hippo_int <- GSE186538_hippo
GSE186538_hippo_int@meta.data$CellType <- GSE186538_hippo_int@meta.data$cluster_1
GSE186538_hippo_int@meta.data$species <- 'human'
columns_to_remove <- c("percent.mt", "samplename", "region", "cluster", "cluster_1")
GSE186538_hippo_int@meta.data <- GSE186538_hippo_int@meta.data[, !colnames(GSE186538_hippo_int@meta.data) %in% columns_to_remove]

head(GSE186538_hippo_int@meta.data)
table(GSE186538_hippo_int@meta.data$orig.ident)
table(GSE186538_hippo_int@meta.data$CellType)
