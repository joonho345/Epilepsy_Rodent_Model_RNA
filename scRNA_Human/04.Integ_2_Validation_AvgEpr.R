library(Seurat)
library(Matrix)

setwd("/home/joonho345/1_Epilepsy_RNA/scRNA_Human/")


#########################################################
#### Load Seurat Object (Deconv) ####
integration_name <- "GSE160189_186538_Integration"
seurat_rds <- "/data/project/1_Epilepsy_RNA/scRNA_Human/Deconv_objects/GSE160189_186538_orig_2_RPCA_1.rds"
seurat_combined <- readRDS(seurat_rds)

# group.by column for averaging
group_by <- "CellType1"

# Output base
output_base_dir <- "/home/joonho345/1_Epilepsy_RNA/scRNA_Human/02.GSE160189_186538_Integration_Seurat_Deconv"

# Integrated DEG gene lists (6 txt files: Up/Down x 3 unions)
deg_list_dir <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/07.DESeq/03.DESeq_sig_genelist"
Target_group_name <- "1"
deg_files <- c(
  file.path(deg_list_dir, paste0(Target_group_name, "_Sig_GeneList_NT_IC_UP.txt")),
  file.path(deg_list_dir, paste0(Target_group_name, "_Sig_GeneList_NT_IC_DOWN.txt")),
  file.path(deg_list_dir, paste0(Target_group_name, "_Sig_GeneList_NI_GO_UP.txt")),
  file.path(deg_list_dir, paste0(Target_group_name, "_Sig_GeneList_NI_GO_DOWN.txt")),
  file.path(deg_list_dir, paste0(Target_group_name, "_Sig_GeneList_NG_GG_GO_UP.txt")),
  file.path(deg_list_dir, paste0(Target_group_name, "_Sig_GeneList_NG_GG_GO_DOWN.txt"))
)

#########################################################
#### Output directory ####
out_dir <- paste0(output_base_dir, "/06.Average_Expr/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#########################################################
#### Pick genes (from 6 DEG txt files) ####
read_genelist <- function(path) {
  if (!file.exists(path)) {
    warning("Missing DEG list file: ", path)
    return(character())
  }
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  unique(x)
}

genes_raw <- unique(unlist(lapply(deg_files, read_genelist), use.names = FALSE))
genes <- genes_raw[genes_raw %in% rownames(seurat_combined)]

#########################################################
#### Average expression (RNA data) ####
DefaultAssay(seurat_combined) <- "RNA"

# Ensure layers are merged if Seurat v5
if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
  try(seurat_combined <- JoinLayers(seurat_combined, assay = "RNA"), silent = TRUE)
}

avg_expr <- AverageExpression(
  seurat_combined,
  features = genes,
  group.by = group_by,
  assays = "RNA",
  verbose = FALSE
)

df_RNA <- as.data.frame(avg_expr$RNA, check.names = FALSE)

#########################################################
#### Ratio (within gene across cell types) ####
df_RNA_ratios <- df_RNA
for (gene in rownames(df_RNA)) {
  total_expression <- sum(df_RNA[gene, ], na.rm = TRUE)
  if (is.finite(total_expression) && total_expression > 0) {
    df_RNA_ratios[gene, ] <- df_RNA[gene, ] / total_expression
  } else {
    df_RNA_ratios[gene, ] <- NA
  }
}

#########################################################
#### Write outputs ####
write.table(
  df_RNA,
  file = paste0(out_dir, "/", integration_name, "_AverageExpr_RNA_by_", group_by, ".tsv"),
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

write.table(
  df_RNA_ratios,
  file = paste0(out_dir, "/", integration_name, "_AverageExpr_RNA_ratio_by_", group_by, ".tsv"),
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)
