library(ggplot2)

#### Define TPM matrices to process ####
tpm_matrices <- c(
  "merged_matrix_TPM_M",
  "adjusted_merged_matrix_TPM_M"
)

#### Read coldata (once, before loop) ####
coldata_df <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/02.Quantification/filtered_coldata_M.txt",
                         sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE, 
                         fill = TRUE, quote = "", comment.char = "")

#### Define colors (once, before loop) ####
cohort_colors <- c(
  "GSE72402"  = "#FDE3A7", # Lighter orange
  "GSE99577"  = "#FAB78B", # Lighter red-orange
  "GSE148028" = "#F58A72",  # Lighter deep red
  "GSE198498" = "#C48EDC",  # Lighter purple
  "GSE205373" = "#C8749F",  # Lighter magenta
  "GSE213393" = "#8CA6E5",  # Lighter blue
  "GSE241219" = "#85C4C9",   # Lighter teal
  "N/A" = "#FFC6D6"  # Lighter pink
)

#### Plot settings ####
axis_text_size <- 12
legend_text_size <- 10
title_text_size <- 14
title_face <- "bold"
font_family <- "Arial"
plot_width <- 8
plot_height <- 6
plot_unit <- 'in'
plot_dpi <- 300

# Create output directory if it doesn't exist
output_dir <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/04.Clustering/01.TPM_PCA"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#### Loop through each TPM matrix ####
cat("=== Processing TPM matrices for PCA ===\n\n")

for (tpm_name in tpm_matrices) {
  cat("Processing:", tpm_name, "\n")
  
  # Read TPM normalized matrix
  tpm_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/03.Normalization/", tpm_name, ".txt")
  
  if (!file.exists(tpm_file)) {
    cat("  ⚠ Warning: File not found:", tpm_file, "\n")
    cat("  Skipping...\n\n")
    next
  }
  
  tpm_matrix <- read.table(tpm_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  
  # Ensure samples match
  common_samples <- intersect(rownames(coldata_df), colnames(tpm_matrix))
  tpm_matrix <- tpm_matrix[, common_samples]
  coldata_subset <- coldata_df[common_samples, , drop = FALSE]
  
  cat("  TPM matrix dimensions:", dim(tpm_matrix), "\n")
  cat("  Number of samples:", length(common_samples), "\n")
  
  #### PCA analysis using prcomp ####
  # Use matrix as-is (genes as rows, samples as columns) - matching old script structure
  pca_data <- as.matrix(tpm_matrix)
  
  # Perform PCA
  pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
  
  # Calculate variance explained
  summary_df <- summary(pca_result)$importance
  variance_explained_pct <- round(summary_df["Proportion of Variance", ] * 100, 2)
  
  # Create PCA data frame using $rotation (sample coordinates when genes are rows)
  pca_df <- as.data.frame(pca_result$rotation)
  pca_df$Sample <- rownames(pca_df)
  
  # Add metadata
  pca_df$GEO <- coldata_subset[pca_df$Sample, "GEO"]
  pca_df$Treatment_Lateral <- coldata_subset[pca_df$Sample, "Treatment_Lateral"]
  pca_df$Sequencing <- coldata_subset[pca_df$Sample, "Sequencing"]
  
  #### Create PCA plot ####
  pca_plot <- ggplot(data = pca_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = GEO), size = 4.0, alpha = 0.8) +
    scale_colour_manual(values = cohort_colors) +
    labs(x = paste0("PC1 (", variance_explained_pct[1], "%)"),
         y = paste0("PC2 (", variance_explained_pct[2], "%)")) +
         #title = paste0("PCA - TPM Normalized (", tpm_name, ")")) +
    theme_minimal() +
    theme(
      #plot.title = element_text(size = title_text_size + 2, family = font_family, face = title_face, hjust = 0.5),
      axis.title = element_text(size = title_text_size, family = font_family, face = title_face),
      axis.text.x = element_text(size = axis_text_size, family = font_family),
      axis.text.y = element_text(size = axis_text_size, family = font_family), 
      #legend.title = element_text(size = legend_text_size, family = font_family),
      #legend.text = element_text(size = legend_text_size, family = font_family), 
      legend.position = "none",
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA), 
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  #### Save plot ####
  file_name <- paste0(output_dir, "/TPM_PCA_M_", tpm_name, ".png")
  ggsave(file_name, plot = pca_plot, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)
  
  cat("  ✓ PCA plot saved:", file_name, "\n")
  cat("  PC1 variance explained:", variance_explained_pct[1], "%\n")
  cat("  PC2 variance explained:", variance_explained_pct[2], "%\n")
  cat("\n")
}

cat("=== All TPM PCA plots completed for Mouse ===\n")
