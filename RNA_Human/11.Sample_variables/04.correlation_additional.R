# Load required libraries
library(ggplot2)
library(dplyr)

# Load human data
human_l2fc_df <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/11.Sample_variables/02.deseq/DESeq_ALL.txt",
                            sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

# Define human groups to plot
Human_groups <- c('HS', 'withoutHS', 'IP_Sei', 'IP_Other', 'NO')

# Check if animal data exists in environment
if (!exists("DESeq_l2fc_df_M") || !exists("DESeq_l2fc_df_R")) {
  stop("Animal DESeq data not found. Please run 00.Matrix_M.R first to load DESeq_l2fc_df_M and DESeq_l2fc_df_R")
}

# Convert to numeric
human_l2fc_df[] <- lapply(human_l2fc_df, as.numeric)

############## Robust Statistics Function ##############
# Function for bootstrapping and permutation testing
get_robust_stats <- function(human_vec, animal_vec, n_perm = 1000) {
  # Remove NA pairs
  complete_idx <- complete.cases(human_vec, animal_vec)
  human_vec <- human_vec[complete_idx]
  animal_vec <- animal_vec[complete_idx]
  
  # Check if we have enough data points
  if (length(human_vec) < 3) {
    return(list(r = NA, p_perm = NA, ci_low = NA, ci_high = NA))
  }
  
  # 1. Observed Correlation
  obs_cor <- cor(human_vec, animal_vec, use = "complete.obs")
  
  # 2. Permutation p-value
  # Shuffle one vector n times and see how often random cor >= observed cor
  perm_cors <- replicate(n_perm, {
    cor(sample(human_vec), animal_vec, use = "complete.obs")
  })
  perm_p <- sum(abs(perm_cors) >= abs(obs_cor)) / n_perm
  
  # 3. Bootstrap CI (95%)
  boot_cors <- replicate(n_perm, {
    idx <- sample(1:length(human_vec), replace = TRUE)
    cor(human_vec[idx], animal_vec[idx], use = "complete.obs")
  })
  boot_ci <- quantile(boot_cors, probs = c(0.025, 0.975), na.rm = TRUE)
  
  return(list(r = obs_cor, p_perm = perm_p, ci_low = boot_ci[1], ci_high = boot_ci[2]))
}

############## Calculate correlation matrices with robust statistics ##############
# Mouse
Species <- 'Mouse'
orthologs_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_mm_17128.txt"
l2fc_df <- DESeq_l2fc_df_M
l2fc_df[] <- lapply(l2fc_df, as.numeric)

ortholog_one2one <- read.table(orthologs_path, sep = "\t", header = TRUE, fill = TRUE)
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one$Gene.name), ]
Species_Gene_name <- paste0(Species, '.gene.name')
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one[[Species_Gene_name]]), ]
human_genes <- ortholog_one2one$Gene.name
animal_genes <- ortholog_one2one[[Species_Gene_name]]

# Initialize matrices
correlation_matrix_M <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))
pvalue_matrix_M <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))
perm_pvalue_matrix_M <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))
bootstrap_ci_low_M <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))
bootstrap_ci_high_M <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))

rownames(correlation_matrix_M) <- Human_groups
colnames(correlation_matrix_M) <- colnames(l2fc_df)
rownames(pvalue_matrix_M) <- Human_groups
colnames(pvalue_matrix_M) <- colnames(l2fc_df)
rownames(perm_pvalue_matrix_M) <- Human_groups
colnames(perm_pvalue_matrix_M) <- colnames(l2fc_df)
rownames(bootstrap_ci_low_M) <- Human_groups
colnames(bootstrap_ci_low_M) <- colnames(l2fc_df)
rownames(bootstrap_ci_high_M) <- Human_groups
colnames(bootstrap_ci_high_M) <- colnames(l2fc_df)

n_perm <- 1000

cat("Calculating correlations for Mouse...\n")
for (Human_group in Human_groups) {
  cat("  Processing", Human_group, "...\n")
  human_data_filtered <- human_l2fc_df[rownames(human_l2fc_df) %in% human_genes, ]
  animal_data_filtered <- l2fc_df[rownames(l2fc_df) %in% animal_genes, ]
  rownames(animal_data_filtered) <- human_genes[match(rownames(animal_data_filtered), animal_genes)]
  
  common_genes <- intersect(rownames(human_data_filtered), rownames(animal_data_filtered))
  human_data_filtered <- human_data_filtered[common_genes, ]
  animal_data_filtered <- animal_data_filtered[common_genes, ]
  
  human_l2fc <- human_data_filtered[[Human_group]]
  animal_l2fc <- animal_data_filtered
  
  # Calculate correlation using robust statistics
  robust_results <- lapply(animal_l2fc, function(x) get_robust_stats(human_l2fc, x, n_perm = n_perm))
  
  correlation_matrix_M[Human_group, ] <- sapply(robust_results, function(x) x$r)
  pvalue_matrix_M[Human_group, ] <- sapply(animal_l2fc, function(x) cor.test(human_l2fc, x, use = "complete.obs")$p.value)
  perm_pvalue_matrix_M[Human_group, ] <- sapply(robust_results, function(x) x$p_perm)
  bootstrap_ci_low_M[Human_group, ] <- sapply(robust_results, function(x) x$ci_low)
  bootstrap_ci_high_M[Human_group, ] <- sapply(robust_results, function(x) x$ci_high)
}

# Rat
Species <- 'Rat'
orthologs_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_rn_16211.txt"
l2fc_df <- DESeq_l2fc_df_R
l2fc_df[] <- lapply(l2fc_df, as.numeric)

ortholog_one2one <- read.table(orthologs_path, sep = "\t", header = TRUE, fill = TRUE)
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one$Gene.name), ]
Species_Gene_name <- paste0(Species, '.gene.name')
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one[[Species_Gene_name]]), ]
human_genes <- ortholog_one2one$Gene.name
animal_genes <- ortholog_one2one[[Species_Gene_name]]

# Initialize matrices
correlation_matrix_R <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))
pvalue_matrix_R <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))
perm_pvalue_matrix_R <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))
bootstrap_ci_low_R <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))
bootstrap_ci_high_R <- matrix(nrow = length(Human_groups), ncol = ncol(l2fc_df))

rownames(correlation_matrix_R) <- Human_groups
colnames(correlation_matrix_R) <- colnames(l2fc_df)
rownames(pvalue_matrix_R) <- Human_groups
colnames(pvalue_matrix_R) <- colnames(l2fc_df)
rownames(perm_pvalue_matrix_R) <- Human_groups
colnames(perm_pvalue_matrix_R) <- colnames(l2fc_df)
rownames(bootstrap_ci_low_R) <- Human_groups
colnames(bootstrap_ci_low_R) <- colnames(l2fc_df)
rownames(bootstrap_ci_high_R) <- Human_groups
colnames(bootstrap_ci_high_R) <- colnames(l2fc_df)

cat("Calculating correlations for Rat...\n")
for (Human_group in Human_groups) {
  cat("  Processing", Human_group, "...\n")
  human_data_filtered <- human_l2fc_df[rownames(human_l2fc_df) %in% human_genes, ]
  animal_data_filtered <- l2fc_df[rownames(l2fc_df) %in% animal_genes, ]
  rownames(animal_data_filtered) <- human_genes[match(rownames(animal_data_filtered), animal_genes)]
  
  common_genes <- intersect(rownames(human_data_filtered), rownames(animal_data_filtered))
  human_data_filtered <- human_data_filtered[common_genes, ]
  animal_data_filtered <- animal_data_filtered[common_genes, ]
  
  human_l2fc <- human_data_filtered[[Human_group]]
  animal_l2fc <- animal_data_filtered
  
  # Calculate correlation using robust statistics
  robust_results <- lapply(animal_l2fc, function(x) get_robust_stats(human_l2fc, x, n_perm = n_perm))
  
  correlation_matrix_R[Human_group, ] <- sapply(robust_results, function(x) x$r)
  pvalue_matrix_R[Human_group, ] <- sapply(animal_l2fc, function(x) cor.test(human_l2fc, x, use = "complete.obs")$p.value)
  perm_pvalue_matrix_R[Human_group, ] <- sapply(robust_results, function(x) x$p_perm)
  bootstrap_ci_low_R[Human_group, ] <- sapply(robust_results, function(x) x$ci_low)
  bootstrap_ci_high_R[Human_group, ] <- sapply(robust_results, function(x) x$ci_high)
}

# Combine matrices
correlation_matrix_W <- cbind(correlation_matrix_M, correlation_matrix_R)
pvalue_matrix_W <- cbind(pvalue_matrix_M, pvalue_matrix_R)
perm_pvalue_matrix_W <- cbind(perm_pvalue_matrix_M, perm_pvalue_matrix_R)
bootstrap_ci_low_W <- cbind(bootstrap_ci_low_M, bootstrap_ci_low_R)
bootstrap_ci_high_W <- cbind(bootstrap_ci_high_M, bootstrap_ci_high_R)

# Normalize column names (convert underscores to dashes)
colnames(correlation_matrix_W) <- gsub("_", "-", colnames(correlation_matrix_W))
colnames(pvalue_matrix_W) <- gsub("_", "-", colnames(pvalue_matrix_W))
colnames(perm_pvalue_matrix_W) <- gsub("_", "-", colnames(perm_pvalue_matrix_W))
colnames(bootstrap_ci_low_W) <- gsub("_", "-", colnames(bootstrap_ci_low_W))
colnames(bootstrap_ci_high_W) <- gsub("_", "-", colnames(bootstrap_ci_high_W))

# Define group order (same as in 02.Dotplot_correlation_WholeGene.R)
group_order <- c(
  "M-KAI-IH-IPSI-A-HA", "M-KAI-IH-IPSI-A-AC", "M-KAI-IH-IPSI-A-IM", "M-KAI-IH-IPSI-A-CR",
  "M-KAI-IH-CON-A-AC", "M-KAI-IH-CON-A-IM", "M-KAI-IH-CON-A-CR", 
  "M-KAI-IA-IPSI-A-AC", "M-KAI-IA-IPSI-A-IM", "M-KAI-IA-IPSI-A-CR",
  "M-KAI-IP-A-HA", "R-KAI-IP-A-CR", "R-KAI-SUB-I-IM", 
  "M-PILO-IP-A-HA", "M-PILO-IP-A-AC", "M-PILO-IP-A-IM", "M-PILO-IP-A-CR", 
  "R-PILO-IP-A-CR", "R-PPS-A-AC", "R-PPS-A-DOFS", "R-PPS-A-IM", "R-PPS-A-CR", 
  "R-AMG-IPSI-A-CR", "R-TBI-IPSI-A-CR"
)

# Filter to only include groups that exist in the correlation matrix
group_order <- group_order[group_order %in% colnames(correlation_matrix_W)]

# Create filtered dataframes for lines (same grouping as reference script)
line_group_1 <- group_order[group_order %in% c("M-KAI-IH-IPSI-A-HA", "M-KAI-IH-IPSI-A-AC", "M-KAI-IH-IPSI-A-IM", "M-KAI-IH-IPSI-A-CR",
                                                "M-KAI-IH-CON-A-AC", "M-KAI-IH-CON-A-IM", "M-KAI-IH-CON-A-CR")]
line_group_2 <- group_order[group_order %in% c("M-KAI-IA-IPSI-A-AC", "M-KAI-IA-IPSI-A-IM", "M-KAI-IA-IPSI-A-CR")]
line_group_3 <- group_order[group_order %in% c("M-KAI-IP-A-HA", "R-KAI-IP-A-CR")]
line_group_4 <- group_order[group_order %in% c("M-PILO-IP-A-HA", "M-PILO-IP-A-AC", "M-PILO-IP-A-IM", "M-PILO-IP-A-CR", "R-PILO-IP-A-CR")]
line_group_5 <- group_order[group_order %in% c("R-PPS-A-AC", "R-PPS-A-DOFS", "R-PPS-A-IM", "R-PPS-A-CR")]

# Loop through each human group
for (Target_comparison in Human_groups) {
  cat("\nProcessing plots for:", Target_comparison, "\n")
  
  # Extract data for this human group
  correlation_data <- as.numeric(correlation_matrix_W[Target_comparison, ])
  bootstrap_ci_low_data <- as.numeric(bootstrap_ci_low_W[Target_comparison, ])
  bootstrap_ci_high_data <- as.numeric(bootstrap_ci_high_W[Target_comparison, ])
  perm_pvalue_data <- as.numeric(perm_pvalue_matrix_W[Target_comparison, ])
  
  # Create maps
  correlation_map <- setNames(correlation_data, colnames(correlation_matrix_W))
  bootstrap_ci_low_map <- setNames(bootstrap_ci_low_data, colnames(bootstrap_ci_low_W))
  bootstrap_ci_high_map <- setNames(bootstrap_ci_high_data, colnames(bootstrap_ci_high_W))
  perm_pvalue_map <- setNames(perm_pvalue_data, colnames(perm_pvalue_matrix_W))
  
  # Create dataframe
  transformed_df <- data.frame(
    Group = names(correlation_map),
    Correlation = unname(correlation_map),
    CI_low = unname(bootstrap_ci_low_map[match(names(correlation_map), names(bootstrap_ci_low_map))]),
    CI_high = unname(bootstrap_ci_high_map[match(names(correlation_map), names(bootstrap_ci_high_map))]),
    Perm_pvalue = unname(perm_pvalue_map[match(names(correlation_map), names(perm_pvalue_map))]),
    stringsAsFactors = FALSE
  )
  
  # Create categorical variable for significance (p < 0.05 = significant)
  transformed_df$Significance <- ifelse(transformed_df$Perm_pvalue < 0.05, "Significant", "Non-significant")
  
  # Convert Group names to match group_order format
  transformed_df <- transformed_df %>%
    mutate(Group = gsub("\\.", "-", Group))
  
  # Filter to only include groups in group_order
  transformed_df <- transformed_df %>%
    filter(Group %in% group_order) %>%
    mutate(Group = factor(Group, levels = group_order)) %>%
    arrange(Group)
  
  # Create line groups
  transformed_df_line_1 <- transformed_df %>% filter(Group %in% line_group_1)
  transformed_df_line_2 <- transformed_df %>% filter(Group %in% line_group_2)
  transformed_df_line_3 <- transformed_df %>% filter(Group %in% line_group_3)
  transformed_df_line_4 <- transformed_df %>% filter(Group %in% line_group_4)
  transformed_df_line_5 <- transformed_df %>% filter(Group %in% line_group_5)
  
  ############## plot with ticks ###########
  axis_text_size <- 10
  legend_text_size <- 8
  title_text_size <- 12
  title_face <- "bold"
  font_family <- "Arial"
  plot_width <- 13
  plot_height <- 4
  plot_unit <- 'in'
  plot_dpi <- 300
  
  plot_target <- ggplot(transformed_df, aes(x = Group, y = Correlation)) +
    geom_hline(yintercept = c(0.1, 0.2, 0.3), linetype = "dashed", color = "gray", size = 0.3, alpha = 0.3) +
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.3, size = 0.5, alpha = 0.6, color = "#d73027") +
    geom_line(data = transformed_df_line_1, aes(x = Group, y = Correlation, group = 1), size = 1, alpha = 0.5, color = "#d73027") +
    geom_line(data = transformed_df_line_2, aes(x = Group, y = Correlation, group = 1), size = 1, alpha = 0.5, color = "#d73027") +
    geom_line(data = transformed_df_line_3, aes(x = Group, y = Correlation, group = 1), size = 1, alpha = 0.5, color = "#d73027") +
    geom_line(data = transformed_df_line_4, aes(x = Group, y = Correlation, group = 1), size = 1, alpha = 0.5, color = "#d73027") +
    geom_line(data = transformed_df_line_5, aes(x = Group, y = Correlation, group = 1), size = 1, alpha = 0.5, color = "#d73027") +
    geom_point(aes(color = Significance), size = 3, alpha = 1) +
    scale_color_manual(
      values = c("Significant" = "#d73027", "Non-significant" = "#808080"),
      name = "Permutation\np-value",
      guide = guide_legend(title.position = "top", title.hjust = 0.5)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      axis.title =  element_text(size = title_text_size, family = font_family, face = title_face),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = axis_text_size, family = font_family, face = title_face,
                                 angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = axis_text_size, family = font_family), 
      legend.title = element_text(size = legend_text_size, family = font_family),
      legend.text = element_text(size = legend_text_size, family = font_family), 
      legend.position = "right",
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA), 
      plot.background = element_rect(fill = "white", color = NA), 
      plot.margin = margin(t = 10, r = 10, b = 10, l = 25, unit = "pt")
    )
  
  output_dir <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/11.Sample_variables/04.correlation_dotplot/"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  file_name <- paste0(output_dir, "01.Dotplot_correlation_", Target_comparison, ".png")
  ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)
  cat("  Plot saved to:", file_name, "\n")
  
  ############## plot without ticks ###########
  axis_text_size <- 10
  legend_text_size <- 8
  title_text_size <- 10
  title_face <- "bold"
  font_family <- "Arial"
  plot_width <- 13
  plot_height <- 1.8
  plot_unit <- 'in'
  plot_dpi <- 300
  
  plot_target <- ggplot(transformed_df, aes(x = Group, y = Correlation)) +
    geom_hline(yintercept = c(0, 0.1, 0.2, 0.3), linetype = "dashed", color = "gray", size = 0.3, alpha = 0.3) +
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.3, size = 0.5, alpha = 0.6, color = "#d73027") +
    geom_line(data = transformed_df_line_1, aes(x = Group, y = Correlation, group = 1), size = 1, alpha = 0.5, color = "#d73027") +
    geom_line(data = transformed_df_line_2, aes(x = Group, y = Correlation, group = 1), size = 1, alpha = 0.5, color = "#d73027") +
    geom_line(data = transformed_df_line_3, aes(x = Group, y = Correlation, group = 1), size = 1, alpha = 0.5, color = "#d73027") +
    geom_line(data = transformed_df_line_4, aes(x = Group, y = Correlation, group = 1), size = 1, alpha = 0.5, color = "#d73027") +
    geom_line(data = transformed_df_line_5, aes(x = Group, y = Correlation, group = 1), size = 1, alpha = 0.5, color = "#d73027") +
    geom_point(aes(color = Significance), size = 3, alpha = 1) +
    scale_color_manual(
      values = c("Significant" = "#d73027", "Non-significant" = "#808080"),
      name = "Permutation\np-value",
      guide = guide_legend(title.position = "top", title.hjust = 0.5)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      axis.title =  element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = axis_text_size, family = font_family), 
      legend.title = element_blank(),
      legend.text = element_blank(), 
      legend.position = "none",
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA), 
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  file_name <- paste0(output_dir, "01.Dotplot_correlation_", Target_comparison, "_tick_line.png")
  ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)
  cat("  Plot saved to:", file_name, "\n")
  
  file_name <- paste0(output_dir, "01.Dotplot_correlation_", Target_comparison, "_tick_line.svg")
  ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit, device = 'svg')
  cat("  Plot saved to:", file_name, "\n\n")
}

cat("\nAll plots generated successfully!\n")
