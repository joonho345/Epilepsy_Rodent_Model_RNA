library(ggplot2)

#### Define count matrices to process ####
count_matrices <- c(
  "merged_matrix",
  "adjusted_merged_matrix_1"
)

#### Read coldata (once, before loop) ####
coldata_df <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/02.Quantification/filtered_coldata.txt",
                         sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE, 
                         fill = TRUE, quote = "", comment.char = "")

#### Define colors (once, before loop) ####
cohort_colors <- c(
  "PRJNA1073977" = "#FDE3A7", # Lighter orange
  "PRJNA1077986" = "#FAB78B", # Lighter red-orange
  "PRJNA280563" = "#F58A72",  # Lighter deep red
  "PRJNA290212" = "#C48EDC",  # Lighter purple
  "PRJNA322318" = "#C8749F",  # Lighter magenta
  "PRJNA373909" = "#8CA6E5",  # Lighter blue
  "PRJNA525671" = "#85C4C9",  # Lighter teal
  "PRJNA556159" = "#FFC6D6",  # Lighter pink
  "PRJNA600414" = "#9A74B8",  # Lighter dark purple
  "PRJNA753504" = "#93DCB1",  # Lighter light green
  "PRJNA773413" = "#D8F296",  # Lighter lime green
  "PRJNA787219" = "#F5ED92"   # Lighter yellow
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

for (count_name in count_matrices) {
  # Read raw count matrix
  count_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/02.Quantification/", count_name, ".txt")
  count_matrix <- read.table(count_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  
  # Ensure samples match
  common_samples <- intersect(rownames(coldata_df), colnames(count_matrix))
  count_matrix <- count_matrix[, common_samples]
  coldata_subset <- coldata_df[common_samples, , drop = FALSE]

  #### PCA analysis using prcomp ####
  pca_data <- as.matrix(count_matrix)
  pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
  summary_df <- summary(pca_result)$importance
  variance_explained_pct <- round(summary_df["Proportion of Variance", ] * 100, 2)
  
  # Create PCA data frame using $rotation (sample coordinates when genes are rows)
  pca_df <- as.data.frame(pca_result$rotation)
  pca_df$Sample <- rownames(pca_df)
  
  # Add metadata
  pca_df$PRJNA <- coldata_subset[pca_df$Sample, "PRJNA"]
  pca_df$Diagnosis <- coldata_subset[pca_df$Sample, "Diagnosis"]
  pca_df$Brain_Location <- coldata_subset[pca_df$Sample, "Brain_Location"]
  
  #### Create PCA plot ####
  pca_plot <- ggplot(data = pca_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = PRJNA), size = 4.0, alpha = 0.8) +
    scale_colour_manual(values = cohort_colors) +
    labs(x = paste0("PC1 (", variance_explained_pct[1], "%)"),
         y = paste0("PC2 (", variance_explained_pct[2], "%)")) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = title_text_size, family = font_family, face = title_face),
      axis.text.x = element_text(size = axis_text_size, family = font_family),
      axis.text.y = element_text(size = axis_text_size, family = font_family), 
      legend.position = "none",
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA), 
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  #### Save plot ####
  file_name <- paste0(output_dir, "/Count_PCA_all_", count_name, ".png")
  ggsave(file_name, plot = pca_plot, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)
}

