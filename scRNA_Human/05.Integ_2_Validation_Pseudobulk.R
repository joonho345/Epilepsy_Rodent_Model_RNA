library(Seurat)
library(Matrix)
library(DESeq2)

setwd("/home/joonho345/1_Epilepsy_RNA/scRNA_Human/")

#########################################################
#### Load Seurat Object (Deconv) ####
integration_name <- "GSE160189_186538_Integration"
seurat_rds <- "/data/project/1_Epilepsy_RNA/scRNA_Human/Deconv_objects/GSE160189_186538_orig_2_RPCA_1.rds"
seurat_combined <- readRDS(seurat_rds)

# Output directory
output_base_dir <- "/home/joonho345/1_Epilepsy_RNA/scRNA_Human/02.GSE160189_186538_Integration_Seurat_Deconv"
out_dir <- paste0(output_base_dir, "/07.Pseudobulk_DEG_DESeq2/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Set default assay to RNA for counts extraction
DefaultAssay(seurat_combined) <- "RNA"

#########################################################
#### Helper: aggregate counts per sample per condition ####
aggregate_sample_condition <- function(sample_name, condition_name, meta_df, count_mat) {
  sample_condition_cells <- meta_df$cell_id[
    meta_df$Sample == sample_name & meta_df$group == condition_name
  ]
  if (length(sample_condition_cells) == 0) {
    return(NULL)
  }
  Matrix::rowSums(count_mat[, sample_condition_cells, drop = FALSE])
}

#########################################################
#### Pseudobulk DEG (DESeq2 only) ####
run_pseudobulk_deseq2 <- function(seurat_obj, celltype_name, output_dir) {
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("Processing cell type:", celltype_name, "\n")
  cat(rep("=", 80), "\n", sep = "")

  # Subset to cell type (CellType1 format: ExN_H, InN_H, etc.)
  celltype_filter <- paste0(celltype_name, "_H")
  cells_keep <- WhichCells(seurat_obj, expression = CellType1 == celltype_filter)

  if (length(cells_keep) == 0) {
    cat("No cells found for", celltype_filter, ". Skipping...\n")
    return(NULL)
  }

  pb_obj <- subset(seurat_obj, cells = cells_keep)
  cat("Cells in subset:", ncol(pb_obj), "\n")

  # Check group distribution
  group_table <- table(pb_obj$group)
  cat("Group distribution:\n")
  print(group_table)

  # Need both groups
  if (!("H_case" %in% names(group_table)) || !("H_control" %in% names(group_table))) {
    cat("Warning: Missing H_case or H_control. Skipping DESeq2.\n")
    return(NULL)
  }

  # Prepare metadata
  meta <- pb_obj@meta.data[, c("group", "orig.ident")]
  meta$cell_id <- rownames(meta)
  meta$Sample <- meta$orig.ident

  # Access counts layer (handle Seurat v5 multiple layers)
  rna_assay <- pb_obj[["RNA"]]
  counts_layers <- Layers(rna_assay, search = "counts")

  if (length(counts_layers) > 1) {
    cat("Found", length(counts_layers), "counts layers. Joining layers...\n")
    pb_obj <- JoinLayers(pb_obj, assay = "RNA", layers = "counts")
  }
  counts <- LayerData(pb_obj, assay = "RNA", layer = "counts")

  # Sample-condition combinations
  sample_condition_combos <- unique(meta[, c("Sample", "group")])
  colnames(sample_condition_combos) <- c("sample", "condition")
  cat("\nSample-condition combinations:\n")
  print(sample_condition_combos)

  samples_case <- unique(meta$Sample[meta$group == "H_case"])
  samples_control <- unique(meta$Sample[meta$group == "H_control"])
  cat("Samples in H_case:", length(samples_case), "-", paste(samples_case, collapse = ", "), "\n")
  cat("Samples in H_control:", length(samples_control), "-", paste(samples_control, collapse = ", "), "\n")

  if (length(samples_case) < 1 || length(samples_control) < 1) {
    cat("Error: Not enough samples per condition (need >=1 per group).\n")
    return(NULL)
  }

  # Aggregate counts for each sample-condition combo
  pb_list <- list()
  pb_coldata_rows <- list()

  for (i in 1:nrow(sample_condition_combos)) {
    samp <- sample_condition_combos$sample[i]
    cond <- sample_condition_combos$condition[i]
    col_name <- paste0(samp, "_", cond)

    aggregated_counts <- aggregate_sample_condition(samp, cond, meta, counts)
    if (!is.null(aggregated_counts)) {
      pb_list[[col_name]] <- aggregated_counts
      pb_coldata_rows[[col_name]] <- data.frame(
        sample = samp,
        condition = cond,
        row.names = col_name,
        stringsAsFactors = FALSE
      )
    }
  }

  pb_mat <- do.call(cbind, pb_list)
  coldata <- do.call(rbind, pb_coldata_rows)

  cat("\nPseudobulk matrix:", nrow(pb_mat), "genes x", ncol(pb_mat), "sample-condition combinations\n")

  # Filter low-count genes
  keep_genes <- rowSums(pb_mat) > 10
  pb_mat_filtered <- pb_mat[keep_genes, ]
  cat("Genes after filtering (total count > 10):", nrow(pb_mat_filtered), "\n")

  if (nrow(pb_mat_filtered) == 0) {
    cat("Error: No genes passed filtering. Skipping DESeq2.\n")
    return(NULL)
  }

  #########################################################
  #### DESeq2 ####
  cat("\n--- Running DESeq2 ---\n")
  coldata$condition <- factor(coldata$condition, levels = c("H_control", "H_case"))

  dds <- DESeqDataSetFromMatrix(
    countData = round(pb_mat_filtered),
    colData = coldata,
    design = ~ condition
  )

  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", "H_case", "H_control"))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]
  res_df <- res_df[order(-res_df$log2FoldChange), ]

  out_file <- file.path(output_dir, paste0("DESeq2_", celltype_name, "_case_vs_control.tsv"))
  write.table(res_df, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("DESeq2 results saved to:", out_file, "\n")
  cat("  Significant genes (padj < 0.05):", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")

  return(TRUE)
}

#########################################################
#### Run analysis for each major cell type ####
celltypes <- c("ExN", "InN", "Astro", "Micro", "Oligo", "OPC", "Endo")

results_summary <- data.frame(
  celltype = character(),
  cells = integer(),
  samples_case = integer(),
  samples_control = integer(),
  genes_tested = integer(),
  deseq_sig = integer(),
  stringsAsFactors = FALSE
)

for (celltype in celltypes) {
  result <- tryCatch(
    run_pseudobulk_deseq2(seurat_combined, celltype, out_dir),
    error = function(e) {
      cat("Error in DESeq2 analysis for", celltype, ":", e$message, "\n")
      NULL
    }
  )

  celltype_filter <- paste0(celltype, "_H")
  cells_subset <- WhichCells(seurat_combined, expression = CellType1 == celltype_filter)
  subset_obj <- if (length(cells_subset) > 0) subset(seurat_combined, cells = cells_subset) else NULL

  samples_case <- if (!is.null(subset_obj)) length(unique(subset_obj$orig.ident[subset_obj$group == "H_case"])) else 0
  samples_control <- if (!is.null(subset_obj)) length(unique(subset_obj$orig.ident[subset_obj$group == "H_control"])) else 0

  deseq_file <- file.path(out_dir, paste0("DESeq2_", celltype, "_case_vs_control.tsv"))
  genes_tested <- 0
  deseq_sig <- 0
  if (file.exists(deseq_file)) {
    deseq_res <- read.table(deseq_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    genes_tested <- nrow(deseq_res)
    deseq_sig <- sum(deseq_res$padj < 0.05, na.rm = TRUE)
  }

  results_summary <- rbind(
    results_summary,
    data.frame(
      celltype = celltype,
      cells = length(cells_subset),
      samples_case = samples_case,
      samples_control = samples_control,
      genes_tested = genes_tested,
      deseq_sig = deseq_sig
    )
  )
}

