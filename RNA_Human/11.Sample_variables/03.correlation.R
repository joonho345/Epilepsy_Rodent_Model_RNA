library(dplyr)
library(ggplot2)


human_l2fc_df <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/11.Sample_variables/02.deseq/DESeq_ALL.txt",
                            sep = "\t", header = TRUE, row.names = 1, fill = TRUE)
#Human_groups <- c('HS', 'withoutHS', 'childhood', 'adolescence', 'adult', 'X20', 'X40', 'X60', 'IP', 'NO', 'IP_Sei', 'IP_Other')
Human_groups <- c('HS', 'withoutHS', 'IP_Sei', 'IP_Other', 'NO')
Human_groups_name <- '_variables'


# Choose Species
Species <- 'Mouse'
orthologs_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_mm_17128.txt"
l2fc_df <- DESeq_l2fc_df_M # import data from 00.Matrix_M.R
coldata_df_M <- DESeq_coldata_df_M # import data from 00.Matrix_M.R

human_l2fc_df[] <- lapply(human_l2fc_df, as.numeric)
l2fc_df[] <- lapply(l2fc_df, as.numeric)

##
ortholog_one2one <- read.table(orthologs_path, sep = "\t", header = TRUE, fill = TRUE)
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one$Gene.name), ]
Species_Gene_name <- paste0(Species, '.gene.name')
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one[[Species_Gene_name]]), ]
human_genes <- ortholog_one2one$Gene.name
animal_genes <- ortholog_one2one[[Species_Gene_name]]

### correlation_matrix ###
# Initialize a matrix to store correlation results
correlation_matrix <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))
pvalue_matrix <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))
rownames(correlation_matrix) <- Human_groups
colnames(correlation_matrix) <- colnames(l2fc_df)
rownames(pvalue_matrix) <- Human_groups
colnames(pvalue_matrix) <- colnames(l2fc_df)

for (Human_group in Human_groups) {
  human_data_filtered <- human_l2fc_df[rownames(human_l2fc_df) %in% human_genes, ]
  animal_data_filtered <- l2fc_df[rownames(l2fc_df) %in% animal_genes, ]
  rownames(animal_data_filtered) <- human_genes[match(rownames(animal_data_filtered), animal_genes)]
  
  common_genes <- intersect(rownames(human_data_filtered), rownames(animal_data_filtered))
  human_data_filtered <- human_data_filtered[common_genes, ]
  animal_data_filtered <- animal_data_filtered[common_genes, ]
  
  human_l2fc <- human_data_filtered[[Human_group]]  # Select the relevant human group column
  animal_l2fc <- animal_data_filtered
  
  correlation_results <- sapply(animal_l2fc, function(x) cor(human_l2fc, x, use = "complete.obs"))
  pvalue_results <- sapply(animal_l2fc, function(x) cor.test(human_l2fc, x, use = "complete.obs")$p.value)
  
  correlation_matrix[Human_group, ] <- correlation_results
  pvalue_matrix[Human_group, ] <- pvalue_results
}

correlation_matrix_M <- correlation_matrix
pvalue_matrix_M <- pvalue_matrix



#################
Species <- 'Rat'
orthologs_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_rn_16211.txt"
l2fc_df <- DESeq_l2fc_df_R # import data from 00.Matrix_M.R
coldata_df_R <- DESeq_coldata_df_R # import data from 00.Matrix_M.R

##
ortholog_one2one <- read.table(orthologs_path, sep = "\t", header = TRUE, fill = TRUE)
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one$Gene.name), ]
Species_Gene_name <- paste0(Species, '.gene.name')
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one[[Species_Gene_name]]), ]
human_genes <- ortholog_one2one$Gene.name
animal_genes <- ortholog_one2one[[Species_Gene_name]]

### correlation_matrix ###
# Initialize a matrix to store correlation results
correlation_matrix <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))
pvalue_matrix <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))
rownames(correlation_matrix) <- Human_groups
colnames(correlation_matrix) <- colnames(l2fc_df)
rownames(pvalue_matrix) <- Human_groups
colnames(pvalue_matrix) <- colnames(l2fc_df)

for (Human_group in Human_groups) {
  human_data_filtered <- human_l2fc_df[rownames(human_l2fc_df) %in% human_genes, ]
  animal_data_filtered <- l2fc_df[rownames(l2fc_df) %in% animal_genes, ]
  rownames(animal_data_filtered) <- human_genes[match(rownames(animal_data_filtered), animal_genes)]
  
  common_genes <- intersect(rownames(human_data_filtered), rownames(animal_data_filtered))
  human_data_filtered <- human_data_filtered[common_genes, ]
  animal_data_filtered <- animal_data_filtered[common_genes, ]
  
  human_l2fc <- human_data_filtered[[Human_group]]  # Select the relevant human group column
  animal_l2fc <- animal_data_filtered
  
  correlation_results <- sapply(animal_l2fc, function(x) cor(human_l2fc, x, use = "complete.obs"))
  pvalue_results <- sapply(animal_l2fc, function(x) cor.test(human_l2fc, x, use = "complete.obs")$p.value)
  
  correlation_matrix[Human_group, ] <- correlation_results
  pvalue_matrix[Human_group, ] <- pvalue_results
}

correlation_matrix_R <- correlation_matrix
pvalue_matrix_R <- pvalue_matrix


################
correlation_matrix_W <- cbind(correlation_matrix_M, correlation_matrix_R)
pvalue_matrix_W <- cbind(pvalue_matrix_M, pvalue_matrix_R)

# Load permutation p-value matrices if they exist
perm_pvalue_file_M <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/11.Sample_variables/03.correlation/01.Correlation_Models_Mouse_perm_pvalue.txt"
perm_pvalue_file_R <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/11.Sample_variables/03.correlation/01.Correlation_Models_Rat_perm_pvalue.txt"

if (file.exists(perm_pvalue_file_M) && file.exists(perm_pvalue_file_R)) {
  perm_pvalue_M <- read.table(perm_pvalue_file_M, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  perm_pvalue_R <- read.table(perm_pvalue_file_R, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  perm_pvalue_W <- cbind(perm_pvalue_M, perm_pvalue_R)
  use_perm_pvalue <- TRUE
} else {
  use_perm_pvalue <- FALSE
  warning("Permutation p-value files not found. Using correlation test p-values instead.")
}

# Define two plot configurations
plot_configs <- list(
  list(
    models = c("M_KAI_IH_IPSI_A_CR", "R_PPS_A_CR"),
    suffix = "_CR",
    model_levels = c("R_PPS_A_CR", "M_KAI_IH_IPSI_A_CR")
  ),
  list(
    models = c("M_KAI_IH_IPSI_A_IM", "R_PPS_A_DOFS"),
    suffix = "_BEST",
    model_levels = c("R_PPS_A_DOFS", "M_KAI_IH_IPSI_A_IM")
  )
)

# Generate plots for each configuration
for (config in plot_configs) {
  # Select models
  correlation_matrix_subset <- correlation_matrix_W[, config$models, drop = FALSE]
  pvalue_matrix_subset <- pvalue_matrix_W[, config$models, drop = FALSE]
  
  correlation_matrix_subset <- t(correlation_matrix_subset)
  pvalue_matrix_subset <- t(pvalue_matrix_subset)
  
  # Save correlation and p-value matrices to text files
  output_dir <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/11.Sample_variables/03.correlation/"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save correlation matrix
  corr_file_name <- paste0(output_dir, "03.Correlation_Matrix", config$suffix, ".txt")
  write.table(correlation_matrix_subset, file = corr_file_name, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  cat("Correlation matrix saved to:", corr_file_name, "\n")
  
  # Save p-value matrix
  pval_file_name <- paste0(output_dir, "04.Pvalue_Matrix", config$suffix, ".txt")
  write.table(pvalue_matrix_subset, file = pval_file_name, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  cat("P-value matrix saved to:", pval_file_name, "\n")
  
  cols <- colorRamp2(c(-max(correlation_matrix_subset), 0, max(correlation_matrix_subset)), c("dodgerblue", "white", "firebrick"))
  
  #### Heatmap without clustering ####
  heat_plot <- Heatmap(correlation_matrix_subset,
                       name = "Correlation", 
                       col = cols,  # Use the color scale
                       rect_gp = gpar(col = "white", lwd = 2),  # Style the grid
                       row_title = paste0(Species, "Groups"),
                       column_title = paste0(Species, "Groups"),
                       column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                       cluster_rows = FALSE,  # Disable row clustering
                       cluster_columns = FALSE,  # Disable column clustering
                       show_row_dend = FALSE,  # Disable row dendrogram
                       show_column_dend = FALSE,  # Disable column dendrogram
                       row_names_gp = gpar(fontsize = 12),  # Set font size for row names
                       column_names_gp = gpar(fontsize = 12))  # Set font size for column names
  
  # Save the heatmap
  file_name <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/11.Sample_variables/03.correlation/01.Correlation_Models", config$suffix, ".png")
  CairoPNG(file = file_name, width = 1600, height = 700) 
  draw(heat_plot, heatmap_legend_side = "right", annotation_legend_side = "right", padding = unit(c(10, 10, 10, 40), "mm"))  # Adjust padding for better visualization
  dev.off()
  
  ######################
  correlation_long <- reshape2::melt(correlation_matrix_subset)
  pvalue_long <- reshape2::melt(pvalue_matrix_subset)
  dotplot_data <- merge(correlation_long, pvalue_long, by = c("Var1", "Var2"))
  colnames(dotplot_data) <- c("Model", "Variable", "Correlation", "PValue")
  dotplot_data$Model <- factor(dotplot_data$Model, levels = config$model_levels)  # Set desired order
  
  # Add permutation p-values if available, use them for coloring
  if (use_perm_pvalue) {
    perm_pvalue_subset <- perm_pvalue_W[, config$models, drop = FALSE]
    perm_pvalue_subset <- t(perm_pvalue_subset)
    perm_pvalue_long <- reshape2::melt(perm_pvalue_subset)
    colnames(perm_pvalue_long) <- c("Model", "Variable", "PermPValue")
    dotplot_data <- merge(dotplot_data, perm_pvalue_long, by = c("Model", "Variable"), all.x = TRUE)
    # Use permutation p-value for logP calculation (for coloring)
    dotplot_data$logP <- -log10(ifelse(!is.na(dotplot_data$PermPValue), dotplot_data$PermPValue, dotplot_data$PValue))
  } else {
    # Use correlation test p-value
    dotplot_data$logP <- -log10(dotplot_data$PValue)
  }
  
  # Set size to 0 for non-significant points (p-value >= 0.05)
  dotplot_data <- dotplot_data %>%
    mutate(
      Correlation_size = ifelse(PValue < 0.05, abs(Correlation), 0),
      LogP_size = ifelse(PValue < 0.05, logP, 0)
    )
  
  # Calculate limits for scales
  max_logp <- max(dotplot_data$logP, na.rm = TRUE)
  max_correlation <- max(abs(dotplot_data$Correlation), na.rm = TRUE)
  min_correlation <- min(dotplot_data$Correlation, na.rm = TRUE)
  
  # Set color gradient limits based on configuration
  color_limit <- ifelse(config$suffix == "_BEST", 250, 180)
  
  # Plot
  axis_text_size <- 8
  legend_text_size <- 8
  title_text_size <- 12
  font_family <- "Arial"
  
  # Original plot: size = Correlation, color = LogP (only show significant points)
  plot_target <- ggplot(dotplot_data, aes(x = Variable, y = Model)) +
    geom_hline(aes(yintercept = as.numeric(Model)), linetype = "dashed", color = "gray70") +
    geom_vline(aes(xintercept = as.numeric(Variable)), linetype = "dashed", color = "gray70") +
    geom_point(aes(size = Correlation_size, color = logP)) +
    scale_color_gradient(low = "gray95", high = "firebrick", limits = c(0, color_limit)) +  
    scale_size_continuous(range = c(2, 14), limits = c(0, max_correlation)) +  
    theme_minimal() +
    labs(size = "Correlation",
         color = "-log10(P-value)") +
    theme(
      plot.title = element_blank(),
      axis.title =  element_blank(),
      axis.text.x = element_text(size = axis_text_size, family = font_family),
      axis.text.y = element_text(size = axis_text_size, family = font_family), 
      legend.title = element_text(size = legend_text_size, family = font_family),
      legend.text = element_text(size = legend_text_size, family = font_family), 
      legend.position = "right",
      axis.line = element_line(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA), 
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 150, r = 20, b = 150, l = 20, unit = "pt")
    )
  
  # Swapped plot: size = LogP, color = Correlation (diverging color scale, only show significant points)
  plot_target_swapped <- ggplot(dotplot_data, aes(x = Variable, y = Model)) +
    geom_hline(aes(yintercept = as.numeric(Model)), linetype = "dashed", color = "gray70") +
    geom_vline(aes(xintercept = as.numeric(Variable)), linetype = "dashed", color = "gray70") +
    geom_point(aes(size = LogP_size, color = Correlation)) +
    scale_color_gradient2(low = "dodgerblue", mid = "gray95", high = "firebrick", midpoint = 0, limits = c(min_correlation, max_correlation)) +
    scale_size_continuous(range = c(2, 14), limits = c(0, max_logp)) +  
    theme_minimal() +
    labs(size = "-log10(P-value)",
         color = "Correlation") +
    theme(
      plot.title = element_blank(),
      axis.title =  element_blank(),
      axis.text.x = element_text(size = axis_text_size, family = font_family),
      axis.text.y = element_text(size = axis_text_size, family = font_family), 
      legend.title = element_text(size = legend_text_size, family = font_family),
      legend.text = element_text(size = legend_text_size, family = font_family), 
      legend.position = "right",
      axis.line = element_line(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA), 
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 150, r = 20, b = 150, l = 20, unit = "pt")
    )
  
  plot_width <- 7
  plot_height <- 6
  plot_unit <- 'in'
  plot_dpi <- 300
  
  # Save original plot
  file_name <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/11.Sample_variables/03.correlation/02.Correlation_Models", config$suffix, ".png")
  ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)
  
  file_name <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/11.Sample_variables/03.correlation/02.Correlation_Models", config$suffix, ".svg")
  ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit, device = 'svg')
  
  # Save swapped plot
  file_name <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/11.Sample_variables/03.correlation/02.Correlation_Models", config$suffix, "_swapped.png")
  ggsave(file_name, plot = plot_target_swapped, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)
  
  file_name <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/11.Sample_variables/03.correlation/02.Correlation_Models", config$suffix, "_swapped.svg")
  ggsave(file_name, plot = plot_target_swapped, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit, device = 'svg')
}

