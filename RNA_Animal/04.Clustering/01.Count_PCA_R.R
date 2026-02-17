library(ggplot2)

#### Define count matrices to process ####
count_matrices <- c(
  "merged_matrix_R",
  "adjusted_merged_matrix_R"
)

#### Read coldata (once, before loop) ####
coldata_df <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/02.Quantification/filtered_coldata_R.txt",
                         sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE, 
                         fill = TRUE, quote = "", comment.char = "")

#### Define colors (once, before loop) ####
cohort_colors <- c(
  "GSE50079"  = "#FDE3A7", # Lighter orange
  "GSE75120"  = "#FAB78B", # Lighter red-orange
  "GSE75402"  = "#F58A72",  # Lighter deep red
  "GSE136913" = "#C48EDC",  # Lighter purple
  "GSE137473" = "#C8749F",  # Lighter magenta
  "GSE143555" = "#8CA6E5"   # Lighter blue
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
output_dir <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/04.Clustering/01.Count_PCA"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (count_name in count_matrices) {
  cat("Processing:", count_name, "\n")
  
  # Read raw count matrix
  count_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/02.Quantification/", count_name, ".txt")
  
  if (!file.exists(count_file)) {
    cat("  ⚠ Warning: File not found:", count_file, "\n")
    cat("  Skipping...\n\n")
    next
  }
  
  count_matrix <- read.table(count_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  
  # Ensure samples match
  common_samples <- intersect(rownames(coldata_df), colnames(count_matrix))
  count_matrix <- count_matrix[, common_samples]
  coldata_subset <- coldata_df[common_samples, , drop = FALSE]
  
  cat("  Count matrix dimensions:", dim(count_matrix), "\n")
  cat("  Number of samples:", length(common_samples), "\n")
  
  #### PCA analysis using prcomp ####
  # Use matrix as-is (genes as rows, samples as columns) - matching old script structure
  pca_data <- as.matrix(count_matrix)
  
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
         #title = paste0("PCA - Raw Counts (", count_name, ")")) +
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
  
  #### Save plot without legend ####
  file_name <- paste0(output_dir, "/Count_PCA_R_", count_name, ".png")
  ggsave(file_name, plot = pca_plot, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)
  
  cat("  ✓ PCA plot saved:", file_name, "\n")
  
  #### Create PCA plot with legend ####
  pca_plot_legend <- ggplot(data = pca_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = GEO), size = 4.0, alpha = 0.8) +
    scale_colour_manual(values = cohort_colors) +
    labs(x = paste0("PC1 (", variance_explained_pct[1], "%)"),
         y = paste0("PC2 (", variance_explained_pct[2], "%)")) +
         #title = paste0("PCA - Raw Counts (", count_name, ")")) +
    theme_minimal() +
    theme(
      #plot.title = element_text(size = title_text_size + 2, family = font_family, face = title_face, hjust = 0.5),
      axis.title = element_text(size = title_text_size, family = font_family, face = title_face),
      axis.text.x = element_text(size = axis_text_size, family = font_family),
      axis.text.y = element_text(size = axis_text_size, family = font_family), 
      legend.title = element_text(size = legend_text_size, family = font_family),
      legend.text = element_text(size = legend_text_size, family = font_family), 
      legend.position = "right",
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA), 
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  #### Save plot with legend ####
  file_name_legend <- paste0(output_dir, "/Count_PCA_R_", count_name, "_legend.png")
  ggsave(file_name_legend, plot = pca_plot_legend, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)
  
  cat("  ✓ PCA plot with legend saved:", file_name_legend, "\n")
  cat("  PC1 variance explained:", paste0("PC1 (", variance_explained_pct[1], "%)"), "\n")
  cat("  PC2 variance explained:", paste0("PC2 (", variance_explained_pct[2], "%)"), "\n")
  cat("\n")
}

cat("=== All raw count PCA plots completed for Rat ===\n")
