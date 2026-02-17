# Load required libraries
library(ggplot2)
library(ggrepel)
library(dplyr)

# Define gene sets to analyze
gene_sets <- list(
  DEG_NT = list(
    genes = c("KCNE5", "KCNQ1", "KCNH8", "KCNK5", "PTGS2"),
    colors = c("KCNE5" = "#e5c494", "KCNQ1" = "#66c2a5", "KCNH8" = "#8da0cb", "KCNK5" = "#fc8d62", "PTGS2" = "#4daf4a"),
    y_inter = c(-2, 2, 4, 6),
    y_min = -2.0,
    y_max = 6.0
  ),
  DEG_NI = list(
    genes = c("IL1B", "OSM", "SERPINE1"),
    colors = c("IL1B" = "#e5c494", "OSM" = "#66c2a5", "SERPINE1" = "#8da0cb"),
    y_inter = c(-3, -1.5, 1.5, 3),
    y_min = -5.0,
    y_max = 7.0
  ),
  DEG_NG = list(
    genes = c("SHANK3", "EGR2"),
    colors = c("SHANK3" = "#e5c494", "EGR2" = "#66c2a5"),
    y_inter = c(-3, -1.5, 1.5, 3),
    y_min = -3.0,
    y_max = 5.0
  )
)

# Define treatment configurations
treatment_configs <- list(
  KAI_IH_IPSI = list(
    file_paths = list(
      FILTERED_1_MTLEALL_NL = "/home/joonho345/1_Epilepsy_RNA/RNA_Human/07.DESeq/01.DESeq_files/DESeq_FILTERED_1_MTLEALL_NL.txt",
      M_KAI_IH_IPSI_A_HA = "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/05.DESeq/02.DESeq_files/DESeq_M_KAI_IH_IPSI_A_HA.txt",
      M_KAI_IH_IPSI_A_AC = "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/05.DESeq/02.DESeq_files/DESeq_M_KAI_IH_IPSI_A_AC.txt",
      M_KAI_IH_IPSI_A_IM = "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/05.DESeq/02.DESeq_files/DESeq_M_KAI_IH_IPSI_A_IM.txt",
      M_KAI_IH_IPSI_A_CR = "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/05.DESeq/02.DESeq_files/DESeq_M_KAI_IH_IPSI_A_CR.txt"
    ),
    phase_level = c("HA", "AC", "IM", "CR", "MTLEALL")
  ),
  PPS = list(
    file_paths = list(
      FILTERED_1_MTLEALL_NL = "/home/joonho345/1_Epilepsy_RNA/RNA_Human/07.DESeq/01.DESeq_files/DESeq_FILTERED_1_MTLEALL_NL.txt",
      R_PPS_A_AC = "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/05.DESeq/02.DESeq_files/DESeq_R_PPS_A_AC.txt",
      R_PPS_A_IM = "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/05.DESeq/02.DESeq_files/DESeq_R_PPS_A_IM.txt",
      R_PPS_A_DOFS = "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/05.DESeq/02.DESeq_files/DESeq_R_PPS_A_DOFS.txt",
      R_PPS_A_CR = "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/05.DESeq/02.DESeq_files/DESeq_R_PPS_A_CR.txt"
    ),
    phase_level = c("AC", "DOFS", "IM", "CR", "MTLEALL")
  )
)

# Load orthologs
ortholog_hs_mm_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_mm_17128.txt"
ortholog_hs_rn_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_rn_16211.txt"
ortholog_hs_mm <- read.table(ortholog_hs_mm_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ortholog_hs_rn <- read.table(ortholog_hs_rn_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mouse_to_human <- setNames(ortholog_hs_mm$Gene.name, ortholog_hs_mm$Mouse.gene.name)
rat_to_human <- setNames(ortholog_hs_rn$Gene.name, ortholog_hs_rn$Rat.gene.name)

# Output directory
output_dir <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/10.Gene_Dotplot/03.Dotplot_DEGs/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Loop through each treatment configuration
for (treatment_name in names(treatment_configs)) {
  treatment_config <- treatment_configs[[treatment_name]]
  file_paths <- treatment_config$file_paths
  Phase_level <- treatment_config$phase_level
  
  # Loop through each gene set
  for (gene_set_name in names(gene_sets)) {
    gene_set <- gene_sets[[gene_set_name]]
    genes_to_analyze <- gene_set$genes
    gene_colors <- gene_set$colors
    y_inter <- gene_set$y_inter
    y_min <- gene_set$y_min
    y_max <- gene_set$y_max
    
    cat("Processing:", treatment_name, "-", gene_set_name, "\n")
    
    # Process genes
    all_filtered_dfs <- list()
    for (Gene_name in genes_to_analyze) {
      filtered_dfs <- list()
      for (id in names(file_paths)) {
        df <- read.table(file_paths[[id]], header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
        df$original_gene <- rownames(df)
        if (id == "FILTERED_1_MTLEALL_NL") {
          df$mapped_gene <- df$original_gene  # Human data, no mapping needed
        } else if (grepl("^M", id)) {  # Mouse
          df$mapped_gene <- ifelse(!is.na(mouse_to_human[df$original_gene]), mouse_to_human[df$original_gene], df$original_gene)
        } else if (grepl("^R", id)) {  # Rat
          df$mapped_gene <- ifelse(!is.na(rat_to_human[df$original_gene]), rat_to_human[df$original_gene], df$original_gene)
        } else {
          df$mapped_gene <- df$original_gene
        }
        df$Phase <- ifelse(grepl("HA", id), "HA",
                           ifelse(grepl("AC", id), "AC",
                                  ifelse(grepl("IM", id), "IM",
                                         ifelse(grepl("CR", id), "CR",
                                                ifelse(grepl("DOFS", id), "DOFS",
                                                       ifelse(id == "FILTERED_1_MTLEALL_NL", "MTLEALL", "None"))))))
        filtered_df <- df[df$mapped_gene == Gene_name, , drop = FALSE]
        if (nrow(filtered_df) > 0) {
          filtered_df$ID <- id
          filtered_df$Gene <- Gene_name
          filtered_dfs[[id]] <- filtered_df
        }
      }
      result <- do.call(rbind, filtered_dfs)
      all_filtered_dfs[[Gene_name]] <- result
    }
    
    combined_result <- do.call(rbind, all_filtered_dfs)
    if (nrow(combined_result) == 0) {
      cat("  No data found for", gene_set_name, "- skipping\n")
      next
    }
    
    plot_data <- combined_result[combined_result$Phase != "None", ]  # Exclude rows where Phase is "None"
    plot_data$Phase <- factor(plot_data$Phase, levels = Phase_level)
    
    # Create summary table with log2FC values
    log2fc_table <- plot_data %>%
      select(Gene, ID, Phase, original_gene, mapped_gene, log2FoldChange, baseMean, pvalue, padj) %>%
      arrange(Gene, Phase)
    
    # Reorder Phase column to match Phase_level order
    log2fc_table$Phase <- factor(log2fc_table$Phase, levels = Phase_level)
    log2fc_table <- log2fc_table[order(log2fc_table$Gene, log2fc_table$Phase), ]
    
    # Save log2FC values to text file
    txt_file_name <- paste0(output_dir, "Log2FC_", treatment_name, "_", gene_set_name, ".txt")
    write.table(log2fc_table, file = txt_file_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    cat("  Log2FC values saved to:", txt_file_name, "\n")
    
    # Plot the data
    axis_text_size <- 8
    legend_text_size <- 8
    title_text_size <- 12
    title_face <- "bold"
    font_family <- "Arial"
    
    plot_target <- ggplot(plot_data, aes(x = Phase, y = log2FoldChange, group = Gene, color = Gene)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", size = 0.6) + 
      geom_hline(yintercept = y_inter, linetype = "dashed", color = "gray", size = 0.2) +
      geom_vline(xintercept = which(levels(plot_data$Phase) == "MTLEALL") - 0.5, linetype = "dotted", color = "red", size = 0.8) +  
      geom_line(size = 1, alpha = 0.5) +  # Connect dots with lines
      geom_point(size = 3.0, alpha = 1.0) +  # Add dots
      scale_color_manual(values = gene_colors) +  # Use the defined gene colors
      scale_x_discrete(drop = FALSE) +  # Keep all phases on x-axis even if empty
      labs(
        y = "Log2 Fold Change (Case vs. NL)",
        color = "Gene"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_blank(),
        axis.title =  element_text(size = title_text_size, family = font_family, face = title_face),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, family = font_family, face = title_face, angle = 45, vjust = 1.0, hjust = 1.0),
        axis.text.y = element_text(size = axis_text_size, family = font_family), 
        legend.title = element_text(size = legend_text_size, family = font_family),
        legend.text = element_text(size = legend_text_size, family = font_family), 
        legend.position = "right",
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA), 
        plot.margin = margin(t = 10, r = 20, b = 20, l = 20, unit = "pt")
      ) +
      ylim(y_min, y_max)
    
    # Save the plot
    plot_width <- 6
    plot_height <- 5
    plot_unit <- 'in'
    plot_dpi <- 300
    file_name <- paste0(output_dir, "Dotplot_", treatment_name, "_", gene_set_name, ".png")
    ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)
    cat("  Plot saved to:", file_name, "\n\n")
  }
}

