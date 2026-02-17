library(ggplot2)
library(dplyr)
library(svglite)

#### Read outlier files ####
outlier_dir <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/06.Selection/01.TPM_Outliers"

# Try to read the combined all outliers file first (most efficient)
all_outliers_file <- paste0(outlier_dir, "/01.ALL_Outliers_All_Genes.txt")
if (file.exists(all_outliers_file) && file.info(all_outliers_file)$size > 0) {
  tryCatch({
    ALL_outliers <- trimws(read.table(all_outliers_file, stringsAsFactors = FALSE)$V1)
    ALL_outliers <- ALL_outliers[ALL_outliers != ""]  # Remove empty strings
  }, error = function(e) {
    ALL_outliers <- character(0)
  })
} else {
  # Fallback: read individual gene files and combine
  genes <- c("GFAP", "NOTCH1", "AQP4")
  mtle_outliers <- list()
  nl_outliers <- list()
  
  for (gene in genes) {
    mtle_file <- paste0(outlier_dir, "/01.MTLE_Outliers_", gene, ".txt")
    nl_file <- paste0(outlier_dir, "/01.NL_Outliers_", gene, ".txt")
    
    if (file.exists(mtle_file) && file.info(mtle_file)$size > 0) {
      tryCatch({
        mtle_outliers[[gene]] <- trimws(read.table(mtle_file, stringsAsFactors = FALSE)$V1)
        mtle_outliers[[gene]] <- mtle_outliers[[gene]][mtle_outliers[[gene]] != ""]
      }, error = function(e) {
        mtle_outliers[[gene]] <- character(0)
      })
    } else {
      mtle_outliers[[gene]] <- character(0)
    }
    
    if (file.exists(nl_file) && file.info(nl_file)$size > 0) {
      tryCatch({
        nl_outliers[[gene]] <- trimws(read.table(nl_file, stringsAsFactors = FALSE)$V1)
        nl_outliers[[gene]] <- nl_outliers[[gene]][nl_outliers[[gene]] != ""]
      }, error = function(e) {
        nl_outliers[[gene]] <- character(0)
      })
    } else {
      nl_outliers[[gene]] <- character(0)
    }
  }
  
  # Combine all outliers (union of all genes)
  ALL_outliers <- unique(c(
    unlist(mtle_outliers),
    unlist(nl_outliers)
  ))
}

# Read combined MTLE and NL outliers files (for reporting)
all_mtle_file <- paste0(outlier_dir, "/01.MTLE_Outliers_All_Genes.txt")
all_nl_file <- paste0(outlier_dir, "/01.NL_Outliers_All_Genes.txt")

all_mtle_outliers <- if (file.exists(all_mtle_file) && file.info(all_mtle_file)$size > 0) {
  tryCatch({
    tmp <- trimws(read.table(all_mtle_file, stringsAsFactors = FALSE)$V1)
    tmp[tmp != ""]
  }, error = function(e) {
    character(0)
  })
} else {
  character(0)
}

all_nl_outliers <- if (file.exists(all_nl_file) && file.info(all_nl_file)$size > 0) {
  tryCatch({
    tmp <- trimws(read.table(all_nl_file, stringsAsFactors = FALSE)$V1)
    tmp[tmp != ""]
  }, error = function(e) {
    character(0)
  })
} else {
  character(0)
}

# Read common outliers (intersection across all genes)
common_mtle_file <- paste0(outlier_dir, "/01.MTLE_Outliers_Common_All_Genes.txt")
common_nl_file <- paste0(outlier_dir, "/01.NL_Outliers_Common_All_Genes.txt")

common_mtle_outliers <- if (file.exists(common_mtle_file) && file.info(common_mtle_file)$size > 0) {
  tryCatch({
    tmp <- trimws(read.table(common_mtle_file, stringsAsFactors = FALSE)$V1)
    tmp[tmp != ""]
  }, error = function(e) {
    character(0)
  })
} else {
  character(0)
}

common_nl_outliers <- if (file.exists(common_nl_file) && file.info(common_nl_file)$size > 0) {
  tryCatch({
    tmp <- trimws(read.table(common_nl_file, stringsAsFactors = FALSE)$V1)
    tmp[tmp != ""]
  }, error = function(e) {
    character(0)
  })
} else {
  character(0)
}


###############################################################
##################### PCA Plot with Outliers ####################
###############################################################

#### Read TPM matrix ####
tpm_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/03.Normalization/adjusted_merged_matrix_1_TPM.txt"
tpm_matrix <- read.table(tpm_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

#### Read coldata ####
coldata_df <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/02.Quantification/filtered_coldata.txt",
                         sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE, 
                         fill = TRUE, quote = "", comment.char = "")

#### Apply filtering logic (matching 00.Subgrouping_H.R lines 59-69) ####
coldata_df_filtered <- coldata_df[coldata_df$Diagnosis %in% c("MTLEHS", "MTLE", "NL"), ]
coldata_df_filtered$Diagnosis_group <- ifelse(coldata_df_filtered$Diagnosis %in% c("MTLEHS", "MTLE"), "MTLEALL", "NL")
coldata_df_filtered$Diagnosis <- NULL
colnames(coldata_df_filtered)[which(names(coldata_df_filtered) == "Diagnosis_group")] <- "Diagnosis"
coldata_df_filtered <- coldata_df_filtered %>% 
  filter(Brain_Location == "HIPPO") %>% 
  mutate(across(everything(), ~ replace(.x, .x == "", "N/A"))) %>%
  mutate(Diagnosis_Sub = ifelse(Diagnosis_Sub == "TLE", "MTLE", Diagnosis_Sub))

coldata_df_filtered$Diagnosis <- factor(coldata_df_filtered$Diagnosis, levels = c("MTLEALL", "NL"))
if ("Diagnosis_Sub" %in% colnames(coldata_df_filtered)) {
  coldata_df_filtered$Diagnosis_Sub <- factor(coldata_df_filtered$Diagnosis_Sub, 
                                               levels = c('MTLEHS', 'MTLE', 'NL'))
}

# Filter to only H_gene_set_W genes
genes_in_matrix <- rownames(tpm_matrix)
genes_to_keep <- intersect(genes_in_matrix, H_gene_set_W)
tpm_matrix_filtered <- tpm_matrix[genes_to_keep, , drop = FALSE]

# Filter samples
filtered_samples <- intersect(rownames(coldata_df_filtered), colnames(tpm_matrix_filtered))
tpm_matrix_filtered <- tpm_matrix_filtered[, filtered_samples, drop = FALSE]
coldata_subset <- coldata_df_filtered[filtered_samples, , drop = FALSE]

# Ensure sample order matches
if (!all(colnames(tpm_matrix_filtered) == rownames(coldata_subset))) {
  coldata_subset <- coldata_subset[colnames(tpm_matrix_filtered), , drop = FALSE]
}

#### Define colors ####
diagnosis_sub_colors <- c(
  "MTLEHS" = "darkseagreen2",
  "MTLE" = "lightsteelblue1",
  "NL" = "antiquewhite2"
)

shape_vector <- c(
  "Whole" = 16,
  "DG" = 15,
  "CA" = 17,
  "SUB" = 18
)

#### Plot settings ####
axis_text_size <- 6
legend_text_size <- 8
title_text_size <- 10
title_face <- "bold"
font_family <- "Arial"
plot_width <- 8
plot_height <- 6
plot_unit <- 'in'
plot_dpi <- 300

# Create output directory
output_dir <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/06.Selection/02.TPM_Outlier_plot"

#### PCA analysis ####
pca_data <- as.matrix(tpm_matrix_filtered)
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)

summary_df <- summary(pca_result)$importance
variance_explained_pct <- round(summary_df["Proportion of Variance", ] * 100, 2)

pca_df <- as.data.frame(pca_result$rotation)
pca_df$Sample <- rownames(pca_df)

# Add metadata
pca_df$Diagnosis <- coldata_subset[pca_df$Sample, "Diagnosis"]
if ("Diagnosis_Sub" %in% colnames(coldata_subset)) {
  pca_df$Diagnosis_Sub <- coldata_subset[pca_df$Sample, "Diagnosis_Sub"]
}
if ("Brain_Location_Sub" %in% colnames(coldata_subset)) {
  pca_df$Brain_Location_Sub <- coldata_subset[pca_df$Sample, "Brain_Location_Sub"]
}

# Mark outliers
pca_df$Outlier <- ifelse(pca_df$Sample %in% ALL_outliers, "Outlier", "Non-Outlier")

#### Create PCA plots with outliers - H_gene_set_W ####
if ("Diagnosis_Sub" %in% colnames(pca_df)) {
  # Separate outlier data for X markers
  pca_df_outlier <- pca_df[pca_df$Outlier == "Outlier", ]
  
  # Plot 1: H_gene_set_W with Diagnosis_Sub only
  # Separate normal and outlier data
  pca_df_normal_diag <- pca_df[pca_df$Outlier == "Non-Outlier", ]
  
  pca_plot_diagnosis <- ggplot(data = pca_df_normal_diag, aes(x = PC1, y = PC2, color = Diagnosis_Sub)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkred", size = 0.3) +
    geom_point(size = 4.0, alpha = 0.8)
  
  # Add outlier points with their colors, then X markers on top
  if (nrow(pca_df_outlier) > 0) {
    pca_plot_diagnosis <- pca_plot_diagnosis +
      geom_point(data = pca_df_outlier,
                 aes(x = PC1, y = PC2, color = Diagnosis_Sub),
                 size = 4.0, alpha = 0.8) +
      geom_point(data = pca_df_outlier,
                 aes(x = PC1, y = PC2),
                 inherit.aes = FALSE,
                 shape = 4,  # X marker
                 size = 3.2,
                 color = "darkred",
                 stroke = 1.2)
  }
  
  pca_plot_diagnosis <- pca_plot_diagnosis +
    scale_color_manual(values = diagnosis_sub_colors, name = "Diagnosis") +
    labs(x = paste0("PC1 (", variance_explained_pct[1], "%)"),
         y = paste0("PC2 (", variance_explained_pct[2], "%)")) +
         #title = paste0("PCA - TPM Normalized (adjusted_merged_matrix_1_TPM) - TLE vs NL (HIPPO, H_gene_set_W) with Outliers")) +
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
  
  file_name_diagnosis <- paste0(output_dir, "/PCA_TPM_TLE_NL_HIPPO_H_gene_set_W_with_Outliers_Diagnosis.png")
  ggsave(file_name_diagnosis, plot = pca_plot_diagnosis, width = plot_width, height = plot_height, 
         dpi = plot_dpi, units = plot_unit)
}


###############################################################
##################### Box Plot with Outliers ####################
###############################################################

#### Define target genes and their settings ####
target_genes_config <- list(
  ACTB = list(
    yinter = c(0, 1000, 2000, 3000, 4000, 5000),
    y_coord = c(0, 5000)
  ),
  GFAP = list(
    yinter = c(0, 2500, 5000, 7500, 10000),
    y_coord = c(0, 10000)
  ),
  NOTCH1 = list(
    yinter = c(0, 25, 50, 75, 100),
    y_coord = c(0, 100)
  ),
  AQP4 = list(
    yinter = c(0, 1000, 2000, 3000, 4000),
    y_coord = c(0, 4000)
  )
)

Group_1 <- "MTLEHS"
Group_2 <- "MTLE"
Group_3 <- "NL"

group_colors <- c("MTLEHS" = "darkseagreen2", "MTLE" = "lightsteelblue1", "NL" = "antiquewhite2")

# Re-read and apply filtering logic (matching 04.TPM_Comparison_3.R)
Target_matrix <- read.table(tpm_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

Target_coldata <- coldata_df[coldata_df$Diagnosis %in% c("MTLEHS", "MTLE", "NL"), ]
Target_coldata$Diagnosis_F <- ifelse(Target_coldata$Diagnosis %in% c("MTLEHS", "MTLE"), "MTLEALL", "NL")
Target_coldata$Diagnosis <- NULL
colnames(Target_coldata)[which(names(Target_coldata) == "Diagnosis_F")] <- "Diagnosis"
Target_coldata <- Target_coldata %>% 
  filter(Brain_Location == "HIPPO") %>% 
  mutate(across(everything(), ~ replace(.x, .x == "", "N/A"))) %>%
  mutate(Diagnosis_Sub = ifelse(Diagnosis_Sub == "TLE", "MTLE", Diagnosis_Sub))

# Filter matrix to match filtered coldata
filtered_samples <- rownames(Target_coldata)
Target_matrix <- Target_matrix[, colnames(Target_matrix) %in% filtered_samples, drop = FALSE]

Target_coldata$sample_name <- rownames(Target_coldata)

#### Loop through target genes ####
for (Target_gene in names(target_genes_config)) {
  if (!Target_gene %in% rownames(Target_matrix)) {
    next
  }
  
  # Get gene-specific settings
  yinter <- target_genes_config[[Target_gene]]$yinter
  y_coord <- target_genes_config[[Target_gene]]$y_coord
  
  # Extract expression data
  gene_expression <- Target_matrix[Target_gene, ]
  gene_expression <- t(gene_expression)
  gene_expression <- as.data.frame(gene_expression)
  gene_expression$sample_name <- rownames(gene_expression)
  colnames(gene_expression)[1] <- "target_gene_expression"
  
  # Merge with metadata
  gene_expression_merged <- merge(gene_expression, Target_coldata, by = "sample_name")
  
  # Use Diagnosis_Sub for grouping (matching 04.TPM_Comparison_3.R)
  gene_expression_filtered <- gene_expression_merged
  if ("Diagnosis_Sub" %in% colnames(gene_expression_filtered)) {
    gene_expression_filtered$Diagnosis_Sub <- factor(gene_expression_filtered$Diagnosis_Sub, 
                                                      levels = c('MTLEHS', 'MTLE', 'NL'))
  }
  
  # Mark outliers
  gene_expression_filtered$highlight <- ifelse(gene_expression_filtered$sample_name %in% ALL_outliers, "highlight", "normal")
  
  # Prepare data for plotting using Diagnosis_Sub
  if ("Diagnosis_Sub" %in% colnames(gene_expression_filtered)) {
    expression_data <- data.frame(
      Expression = gene_expression_filtered$target_gene_expression,
      Group = gene_expression_filtered$Diagnosis_Sub,
      Highlight = gene_expression_filtered$highlight
    )
    expression_data$Group <- factor(expression_data$Group, levels = c(Group_1, Group_2, Group_3))
  } else {
    # Fallback to Diagnosis if Diagnosis_Sub not available
    expression_data <- data.frame(
      Expression = gene_expression_filtered$target_gene_expression,
      Group = gene_expression_filtered$Diagnosis,
      Highlight = gene_expression_filtered$highlight
    )
  }
  
  # Create box plot with outliers
  plot_target <- ggplot(expression_data, aes(x = Group, y = Expression, fill = Group)) +
    geom_hline(yintercept = yinter, linetype = "dashed", color = "gray", size = 0.3) +
    geom_boxplot(width = 0.7, outlier.shape = NA, position = position_dodge(width = 0.3)) +
    geom_jitter(
      data = subset(expression_data, Highlight == "normal"),
      aes(color = Highlight),
      width = 0.2, size = 0.8, alpha = 0.4
    ) +
    geom_jitter(
      data = subset(expression_data, Highlight == "highlight"),
      aes(color = Highlight),
      width = 0.2, size = 1.0, alpha = 1.0, shape = 4, stroke = 1.0
    ) +
    scale_fill_manual(values = group_colors) + 
    scale_color_manual(values = c("highlight" = "darkred", "normal" = "gray")) + 
    labs(y = paste0(Target_gene, " Normalized TPM Expression")) +
    coord_cartesian(ylim = y_coord) + 
    theme_minimal() + 
    theme(
      axis.title =  element_text(size = title_text_size, family = font_family, face = title_face),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = axis_text_size, family = font_family),
      legend.title = element_blank(),
      legend.text = element_blank(),
      legend.position = "none",
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA), 
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # Save plot as PNG
  file_name_png <- paste0(output_dir, "/Boxplot_", Target_gene, "_with_Outliers.png")
  ggsave(file_name_png, plot = plot_target, width = 2.5, height = 4, dpi = plot_dpi, units = plot_unit)
}

