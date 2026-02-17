library(ggplot2)
library(ggpubr)
library(extrafont)
library(dplyr)
library(svglite)

#### Read TPM matrix ####
tpm_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/03.Normalization/adjusted_merged_matrix_1_TPM.txt"
Target_matrix <- read.table(tpm_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

#### Read coldata ####
coldata_df <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/02.Quantification/filtered_coldata.txt",
                         sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE, 
                         fill = TRUE, quote = "", comment.char = "")

#### Filter for TLE (MTLEALL) and NL only - HIPPOCAMPUS ONLY ####
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

#### Plot settings ####
axis_text_size <- 8
legend_text_size <- 8
title_text_size <- 12
title_face <- "bold"
font_family <- "Arial"
plot_width <- 0.9
plot_height <- 1.7
plot_unit <- 'in'
plot_dpi <- 300

# Initialize results storage
results_summary <- data.frame(
  Gene = character(),
  Ttest_MTLEHS_vs_NL_pvalue = numeric(),
  Ttest_MTLE_vs_NL_pvalue = numeric(),
  stringsAsFactors = FALSE
)

for (Target_gene in names(target_genes_config)) {
  # Get gene-specific settings
  yinter <- target_genes_config[[Target_gene]]$yinter
  y_coord <- target_genes_config[[Target_gene]]$y_coord
  
  # Extract the expression data for the target gene
  gene_expression <- Target_matrix[Target_gene, ]
  group1 <- gene_expression[Target_coldata$Diagnosis_Sub == Group_1]
  group2 <- gene_expression[Target_coldata$Diagnosis_Sub == Group_2]
  group3 <- gene_expression[Target_coldata$Diagnosis_Sub == Group_3]
  group1 <- as.numeric(group1)
  group2 <- as.numeric(group2)
  group3 <- as.numeric(group3)
  
  expression_data <- data.frame(Expression = c(group1, group2, group3), 
                                Group = factor(rep(c(Group_1, Group_2, Group_3), 
                                                   c(length(group1), length(group2), length(group3))),
                                               levels = c(Group_1, Group_2, Group_3)))
  
  # Perform t-tests
  t_test_group1_vs_group3 <- t.test(group1, group3, alternative = "two.sided", var.equal = FALSE)
  t_test_group2_vs_group3 <- t.test(group2, group3, alternative = "two.sided", var.equal = FALSE)
  p_value_group1_vs_group3 <- t_test_group1_vs_group3$p.value
  p_value_group2_vs_group3 <- t_test_group2_vs_group3$p.value

  # Store results
  results_summary <- rbind(results_summary, data.frame(
    Gene = Target_gene,
    Ttest_MTLEHS_vs_NL_pvalue = p_value_group1_vs_group3,
    Ttest_MTLE_vs_NL_pvalue = p_value_group2_vs_group3,
    stringsAsFactors = FALSE
  ))
  
  # Create plot WITH axis title and text 
  plot_target_with_axis <- ggplot(expression_data, aes(x = Group, y = Expression, fill = Group)) +
    geom_hline(yintercept = yinter, linetype = "dashed", color = "gray") +
    geom_jitter(color = "gray", size = 0.8, alpha = 0.4, 
                position = position_jitterdodge(jitter.width = 0.40, dodge.width = 0.3)) +
    geom_boxplot(width = 0.7, outlier.shape = NA, position = position_dodge(width = 0.3)) +
    scale_fill_manual(values = group_colors) + 
    scale_color_manual(values = group_colors) + 
    labs(y = paste0(Target_gene, " Normalized TPM Expression")) +
    coord_cartesian(ylim = y_coord) + 
    theme_minimal() + 
    theme(
      axis.title = element_text(size = title_text_size, family = font_family, face = title_face),
      axis.text.x = element_text(size = axis_text_size, family = font_family),
      axis.text.y = element_text(size = axis_text_size, family = font_family),
      legend.title = element_blank(),
      legend.text = element_blank(),
      legend.position = "none",
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA), 
      plot.background = element_rect(fill = "white", color = NA)
    ) 
  
  # Save plot WITH axis as PNG 
  file_name_png_with_axis <- paste0(output_dir, "/04.TPM_Comparison_", Target_gene, "_3_with_axis.png")
  ggsave(file_name_png_with_axis, plot = plot_target_with_axis, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)
}
