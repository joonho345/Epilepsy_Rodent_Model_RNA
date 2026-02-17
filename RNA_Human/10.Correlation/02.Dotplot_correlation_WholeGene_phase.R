# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ggrepel)

human_l2fc_df <- DESeq_l2fc_df_FILTERED
Human_groups <- c('FILTERED_1_MTLEALL_NL', 'FILTERED_2_MTLEHS_NL', 'FILTERED_3_MTLE_NL')
Target_comparison <- 'FILTERED_1_MTLEALL_NL'
Target_group_name <- '1'

# Choose Species
Species <- 'Mouse'
orthologs_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_mm_17128.txt"
l2fc_df <- DESeq_l2fc_df_M # import data from 00.Matrix_M.R
coldata_df_M <- DESeq_coldata_df_M # import data from 00.Matrix_M.R

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

correlation_matrix_W <- correlation_matrix_W[Target_comparison, ]
pvalue_matrix_W <- pvalue_matrix_W[Target_comparison, ]
  
treatment_color_vector <- c("KAI_IH" = "#fc8d62", "KAI_IA" = "#999999", "KAI_IP" = "#e5c494", "PILO_IP" = "#66c2a5",
                            "TBI" = "#e78ac3", "AMG" = "#a6d854", "KAI_SUB" = "#ffd92f", "PPS" = "#8da0cb")

coldata_df_R
coldata_df_M

coldata_df_W <- rbind(coldata_df_M, coldata_df_R)
treatments <- setNames(coldata_df_W$Treatment, coldata_df_W$Group)
phases <- setNames(coldata_df_W$Phase, coldata_df_W$Group)
ages <- setNames(coldata_df_W$Age, coldata_df_W$Group)


# Extract correlation and p-values for FILTERED_1_MTLEALL_NL
correlation_data <- data.frame(
  Group = names(correlation_matrix_W),
  Correlation = correlation_matrix_W,
  Pvalue = -log10(pvalue_matrix_W)
)

# Merge with metadata
combined_dotplot_data <- correlation_data %>%
  left_join(coldata_df_W, by = "Group") %>%
  mutate(
    Color = treatment_color_vector[Treatment],
    Species = ifelse(Species == "M", "Mouse", "Rat")
  )

combined_dotplot_data$Group <- factor(
  combined_dotplot_data$Group, 
  levels = combined_dotplot_data$Group[order(combined_dotplot_data$Correlation, decreasing = TRUE)]
)
# Create display version with hyphens instead of underscores
combined_dotplot_data$Group_display <- gsub("_", "-", combined_dotplot_data$Group)
combined_dotplot_data$Group_display <- factor(
  combined_dotplot_data$Group_display,
  levels = combined_dotplot_data$Group_display[order(combined_dotplot_data$Correlation, decreasing = TRUE)]
)
combined_dotplot_data_HA <- combined_dotplot_data[combined_dotplot_data$Phase == 'hyperacute', ]
combined_dotplot_data_HA$Group_display <- factor(combined_dotplot_data_HA$Group_display, 
                                                  levels = unique(combined_dotplot_data_HA$Group_display[order(combined_dotplot_data_HA$Correlation, decreasing = TRUE)]))
combined_dotplot_data_AC <- combined_dotplot_data[combined_dotplot_data$Phase == 'acute', ]
combined_dotplot_data_AC$Group_display <- factor(combined_dotplot_data_AC$Group_display, 
                                                 levels = unique(combined_dotplot_data_AC$Group_display[order(combined_dotplot_data_AC$Correlation, decreasing = TRUE)]))
combined_dotplot_data_IM <- combined_dotplot_data[combined_dotplot_data$Phase == 'intermediate', ]
combined_dotplot_data_IM$Group_display <- factor(combined_dotplot_data_IM$Group_display, 
                                                 levels = unique(combined_dotplot_data_IM$Group_display[order(combined_dotplot_data_IM$Correlation, decreasing = TRUE)]))
combined_dotplot_data_CR <- combined_dotplot_data[combined_dotplot_data$Phase == 'chronic', ]
combined_dotplot_data_CR$Group_display <- factor(combined_dotplot_data_CR$Group_display, 
                                                 levels = unique(combined_dotplot_data_CR$Group_display[order(combined_dotplot_data_CR$Correlation, decreasing = TRUE)]))


# Define phases and their respective data subsets
phase_list <- list(
  "hyperacute" = combined_dotplot_data_HA,
  "acute" = combined_dotplot_data_AC,
  "intermediate" = combined_dotplot_data_IM,
  "chronic" = combined_dotplot_data_CR
)


############## plot ###########
# Create the dot plot
axis_text_size <- 10
legend_text_size <- 8
title_text_size <- 12
title_face <- "bold"
font_family <- "Arial"
plot_width <- 6
plot_height <- 4
plot_unit <- 'in'
plot_dpi <- 300

for (phase in names(phase_list)) {
  
  plot_data <- phase_list[[phase]]
  plot_target <- ggplot(plot_data, aes(x = Group_display, y = Correlation)) +
    geom_hline(yintercept = c(0.1, 0.2, 0.3), linetype = "dashed", color = "gray", size = 0.3, alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray", size = 0.5, alpha = 0.3) +

    geom_line(aes(group = 1), size = 1, alpha = 0.5, color = "#d73027") +
    geom_point(size = 3, alpha = 0.8, color = "#d73027") +
    
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
  file_name <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/02.Dotplot_correlation/01.Dotplot_correlation_WholeGene_", Target_group_name, "_Phase_", phase, ".png")
  ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)
}

