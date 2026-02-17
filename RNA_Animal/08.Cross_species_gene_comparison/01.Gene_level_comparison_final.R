library(ggplot2)
library(reshape2)

### generation of dfs
# Human (same for both versions)
file_path <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/07.DESeq/01.DESeq_files/DESeq_FILTERED_1_MTLEALL_NL.txt"
DESeq_FILTERED_1_MTLEALL_NL <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
DESeq_FILTERED_1_MTLEALL_NL <- DESeq_FILTERED_1_MTLEALL_NL[, c("log2FoldChange", "padj"), drop = FALSE]

configs <- list(
  list(
    mouse_file = "DESeq_M_KAI_IH_IPSI_A_IM.txt",
    rat_file = "DESeq_R_PPS_A_DOFS.txt",
    prefix = "02.",
    suffix = "_BEST"
  )
)

# Loop through each configuration
for (config in configs) {
  # Mouse
  file_path <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/05.DESeq/02.DESeq_files/", config$mouse_file)
  DESeq_M <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  Species <- 'Mouse'
  orthologs_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_mm_17128.txt"
  ortholog_one2one <- read.table(orthologs_path, sep = "\t", header = TRUE, fill = TRUE)
  ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one$Gene.name), ]
  Species_Gene_name <- paste0(Species, '.gene.name')
  ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one[[Species_Gene_name]]), ]
  mouse_to_human_map <- setNames(ortholog_one2one$Gene.name, ortholog_one2one[[Species_Gene_name]])
  
  DESeq_M$Human_Gene <- rownames(DESeq_M)
  DESeq_M$Human_Gene <- mouse_to_human_map[DESeq_M$Human_Gene]
  DESeq_M <- DESeq_M[!is.na(DESeq_M$Human_Gene), ]
  rownames(DESeq_M) <- DESeq_M$Human_Gene
  DESeq_M <- DESeq_M[, c("log2FoldChange", "padj"), drop = FALSE]
  
  # Rat
  file_path <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/05.DESeq/02.DESeq_files/", config$rat_file)
  DESeq_R <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  Species <- 'Rat'
  orthologs_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_rn_16211.txt"
  ortholog_one2one <- read.table(orthologs_path, sep = "\t", header = TRUE, fill = TRUE)
  ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one$Gene.name), ]
  Species_Gene_name <- paste0(Species, '.gene.name')
  ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one[[Species_Gene_name]]), ]
  rat_to_human_map <- setNames(ortholog_one2one$Gene.name, ortholog_one2one[[Species_Gene_name]])
  
  DESeq_R$Human_Gene <- rownames(DESeq_R)
  DESeq_R$Human_Gene <- rat_to_human_map[DESeq_R$Human_Gene]
  DESeq_R <- DESeq_R[!is.na(DESeq_R$Human_Gene), ]
  rownames(DESeq_R) <- DESeq_R$Human_Gene
  DESeq_R <- DESeq_R[, c("log2FoldChange", "padj"), drop = FALSE]
  
  ### cbind all dfs
  # Count total Human DEGs (before filtering by shared orthologs)
  n_human_total <- nrow(DESeq_FILTERED_1_MTLEALL_NL)
  human_degs_total <- DESeq_FILTERED_1_MTLEALL_NL[DESeq_FILTERED_1_MTLEALL_NL$padj < 0.05 & 
                                                    (DESeq_FILTERED_1_MTLEALL_NL$log2FoldChange > 1 | DESeq_FILTERED_1_MTLEALL_NL$log2FoldChange < -1), ]
  n_human_degs_total <- nrow(human_degs_total)
  
  # Count Human DEGs in each gene set (before filtering by shared orthologs)
  n_NT_degs_total <- sum(rownames(human_degs_total) %in% union(H_gene_set_NT, H_gene_set_IC))
  n_NI_degs_total <- sum(rownames(human_degs_total) %in% union(H_gene_set_NI_GO0050728_173, H_gene_set_NI_GO0050729_153))
  n_NG_degs_total <- sum(rownames(human_degs_total) %in% union(union(H_gene_set_NG_GO0050769_331, H_gene_set_NG_GO0050768_163),
                                                               union(H_gene_set_GG_GO0014014_59, H_gene_set_GG_GO0014015_94)))
  
  n_mouse_mapped <- nrow(DESeq_M)
  n_rat_mapped <- nrow(DESeq_R)
  
  shared_genes <- Reduce(intersect, list(
    rownames(DESeq_FILTERED_1_MTLEALL_NL),
    rownames(DESeq_M),
    rownames(DESeq_R)))
  n_shared_orthologs <- length(shared_genes)
  
  DESeq_FILTERED_1_MTLEALL_NL_shared <- DESeq_FILTERED_1_MTLEALL_NL[shared_genes, , drop = FALSE]
  DESeq_M_shared <- DESeq_M[shared_genes, , drop = FALSE]
  DESeq_R_shared <- DESeq_R[shared_genes, , drop = FALSE]
  
  colnames(DESeq_FILTERED_1_MTLEALL_NL_shared) <- c("Human", "padj_Human")
  colnames(DESeq_M_shared) <- c("Mouse", "padj_Mouse")
  colnames(DESeq_R_shared) <- c("Rat", "padj_Rat")
  
  combined_df <- cbind(
    DESeq_FILTERED_1_MTLEALL_NL_shared,
    DESeq_M_shared,
    DESeq_R_shared
  )
  
  combined_df$Gene <- rownames(combined_df)
  
  # Count human DEGs (padj < 0.05 and |log2FC| > 1) in shared orthologs
  human_degs_shared <- combined_df[combined_df$padj_Human < 0.05 & (combined_df$Human > 1 | combined_df$Human < -1), ]
  n_human_degs_shared <- nrow(human_degs_shared)
  
  # Count how many human DEGs are in each gene set (from shared orthologs)
  n_NT_degs_shared <- sum(human_degs_shared$Gene %in% union(H_gene_set_NT, H_gene_set_IC))
  n_NI_degs_shared <- sum(human_degs_shared$Gene %in% union(H_gene_set_NI_GO0050728_173, H_gene_set_NI_GO0050729_153))
  n_NG_degs_shared <- sum(human_degs_shared$Gene %in% union(union(H_gene_set_NG_GO0050769_331, H_gene_set_NG_GO0050768_163),
                                                             union(H_gene_set_GG_GO0014014_59, H_gene_set_GG_GO0014015_94)))
  
  # Filter by gene sets and DEG criteria
  filtered_df_NT <- combined_df[combined_df$Gene %in% union(H_gene_set_NT, H_gene_set_IC), , drop = FALSE]
  filtered_df_NT <- filtered_df_NT[filtered_df_NT$padj_Human < 0.05 & (filtered_df_NT$Human > 1 | filtered_df_NT$Human < -1), , drop = FALSE]
  n_NT_final <- nrow(filtered_df_NT)
  
  filtered_df_NI <- combined_df[combined_df$Gene %in% union(H_gene_set_NI_GO0050728_173, H_gene_set_NI_GO0050729_153), , drop = FALSE]
  filtered_df_NI <- filtered_df_NI[filtered_df_NI$padj_Human < 0.05 & (filtered_df_NI$Human > 1 | filtered_df_NI$Human < -1), , drop = FALSE]
  n_NI_final <- nrow(filtered_df_NI)
  
  filtered_df_NG <- combined_df[combined_df$Gene %in% union(union(H_gene_set_NG_GO0050769_331, H_gene_set_NG_GO0050768_163),
                                                            union(H_gene_set_GG_GO0014014_59, H_gene_set_GG_GO0014015_94)), , drop = FALSE]
  filtered_df_NG <- filtered_df_NG[filtered_df_NG$padj_Human < 0.05 & (filtered_df_NG$Human > 1 | filtered_df_NG$Human < -1), , drop = FALSE]
  n_NG_final <- nrow(filtered_df_NG)
  
  # Print summary for this configuration
  cat("\n========== Configuration:", config$suffix, "==========\n")
  cat("Total Human DEGs (padj < 0.05 & |log2FC| > 1):", n_human_degs_total, "\n")
  cat("Human DEGs in shared orthologs:", n_human_degs_shared, "\n")
  cat("Genes filtered out (Total DEGs - Shared orthologs DEGs):", n_human_degs_total - n_human_degs_shared, "\n")
  cat("\n--- Filtering by Gene Set Category ---\n")
  cat("NT: Total DEGs =", n_NT_degs_total, "→ In shared orthologs =", n_NT_degs_shared, "→ Final =", n_NT_final, "(filtered out:", n_NT_degs_total - n_NT_degs_shared, "from total,", n_NT_degs_shared - n_NT_final, "from shared)\n")
  cat("NI: Total DEGs =", n_NI_degs_total, "→ In shared orthologs =", n_NI_degs_shared, "→ Final =", n_NI_final, "(filtered out:", n_NI_degs_total - n_NI_degs_shared, "from total,", n_NI_degs_shared - n_NI_final, "from shared)\n")
  cat("NG: Total DEGs =", n_NG_degs_total, "→ In shared orthologs =", n_NG_degs_shared, "→ Final =", n_NG_final, "(filtered out:", n_NG_degs_total - n_NG_degs_shared, "from total,", n_NG_degs_shared - n_NG_final, "from shared)\n")
  cat("========================================\n\n")
  
  # Save summary to file
  summary_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/09.Cross-species_gene_comparison/05.Gene_level_final/", config$prefix, "filtering_summary", config$suffix, ".txt")
  summary_lines <- c(
    paste0("========== Configuration: ", config$suffix, " =========="),
    "",
    "Overall Statistics:",
    paste0("Total Human DEGs (padj < 0.05 & |log2FC| > 1): ", n_human_degs_total),
    paste0("Human DEGs in shared orthologs: ", n_human_degs_shared),
    paste0("Genes filtered out (Total DEGs - Shared orthologs DEGs): ", n_human_degs_total - n_human_degs_shared),
    "",
    "Filtering by Gene Set Category:",
    paste0("NT: Total DEGs = ", n_NT_degs_total, " → In shared orthologs = ", n_NT_degs_shared, " → Final = ", n_NT_final, " (filtered out: ", n_NT_degs_total - n_NT_degs_shared, " from total, ", n_NT_degs_shared - n_NT_final, " from shared)"),
    paste0("NI: Total DEGs = ", n_NI_degs_total, " → In shared orthologs = ", n_NI_degs_shared, " → Final = ", n_NI_final, " (filtered out: ", n_NI_degs_total - n_NI_degs_shared, " from total, ", n_NI_degs_shared - n_NI_final, " from shared)"),
    paste0("NG: Total DEGs = ", n_NG_degs_total, " → In shared orthologs = ", n_NG_degs_shared, " → Final = ", n_NG_final, " (filtered out: ", n_NG_degs_total - n_NG_degs_shared, " from total, ", n_NG_degs_shared - n_NG_final, " from shared)"),
    "",
    "========================================"
  )
  writeLines(summary_lines, summary_file)
  
  # Define target dataframes and names
  target_configs <- list(
    list(df = filtered_df_NT, name = 'NT'),
    list(df = filtered_df_NI, name = 'NI'),
    list(df = filtered_df_NG, name = 'NG')
  )
  
  # Loop through each target dataframe
  for (target_config in target_configs) {
    target_df <- target_config$df
    name <- target_config$name
    
    # Transform filtered_df into long format for ggplot
    df_long <- melt(target_df, id.vars = "Gene", 
                    measure.vars = c("Human", "Mouse", "Rat"),
                    variable.name = "Species", value.name = "log2FoldChange")
    
    # Create padj long format for significance annotation
    padj_long <- melt(target_df, id.vars = "Gene",
                      measure.vars = c("padj_Human", "padj_Mouse", "padj_Rat"),
                      variable.name = "padj_Species", value.name = "padj")
    padj_long$Species <- ifelse(padj_long$padj_Species == "padj_Human", "Human",
                               ifelse(padj_long$padj_Species == "padj_Mouse", "Mouse", "Rat"))
    
    # Merge padj information into df_long
    df_long <- merge(df_long, padj_long[, c("Gene", "Species", "padj")], by = c("Gene", "Species"), all.x = TRUE)
    
    # Create significance flag for Mouse and Rat (Human is already filtered to be significant)
    df_long$is_significant <- ifelse(df_long$Species == "Human", TRUE,
                                    df_long$padj < 0.05 & !is.na(df_long$padj))
    
    # Order genes by 'Human' log2FoldChange values (descending)
    ordered_genes <- target_df$Gene[order(target_df$Human, decreasing = TRUE)]
    
    # Save gene list to text file
    gene_list_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/09.Cross-species_gene_comparison/05.Gene_level_final/", config$prefix, "gene_list_", name, config$suffix, ".txt")
    writeLines(ordered_genes, gene_list_file)
    
    # Ensure df_long$Gene is a factor with levels explicitly set in the correct order
    df_long$Gene <- factor(df_long$Gene, levels = rev(ordered_genes))
    
    # Generate numerical y-axis mapping based on the new order
    df_long$Gene_Num <- as.numeric(factor(df_long$Gene, levels = rev(ordered_genes)))
    
    
    # Define colors for species
    species_color <- c("Human" = "#d73027", "Mouse" = "#313695", "Rat" = "#1a9850")
    
    # Define plot aesthetics
    axis_text_size <- 8
    legend_text_size <- 8
    title_text_size <- 12
    title_face <- "bold"
    font_family <- "Arial"
    plot_width <- 6
    row_height <- 0.27
    plot_unit <- 'in'
    plot_dpi <- 300
    num_genes <- length(unique(df_long$Gene))
    plot_height <- num_genes * row_height
    
    plot_target <- ggplot(df_long, aes(x = log2FoldChange, y = Gene_Num, color = Species)) +
      geom_rect(aes(xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf), fill = "grey96", inherit.aes = FALSE) +
      geom_vline(xintercept = c(-2.5, -1, 1, 2.5), linetype = "dashed", color = "gray", size = 0.5) +
      geom_segment(aes(x = 0, xend = log2FoldChange, y = Gene_Num, yend = Gene_Num), 
                   color = "gray70", size = 1.0) +
      geom_vline(xintercept = 0, color = "black", size = 0.8) +
      geom_point(size = 4, position = position_dodge(width = 0.6)) + # Dots with separation for species
      # Add gray border for significant Mouse and Rat dots (padj < 0.05)
      geom_point(data = subset(df_long, is_significant & Species != "Human"),
                 aes(x = log2FoldChange, y = Gene_Num),
                 shape = 21, size = 4, stroke = 1.5, color = "gray50", alpha = 0.8, fill = NA,
                 position = position_dodge(width = 0.6)) +
      
      scale_color_manual(values = species_color) +
      scale_y_continuous(
        breaks = df_long$Gene_Num, 
        labels = df_long$Gene) +   # **Ensures correct y-axis labels**
      
      labs(x = "Log2 Fold Change (Case vs. Control)") +
      theme_minimal() +
      theme(
        plot.title = element_blank(),
        axis.title =  element_text(size = title_text_size, family = font_family, face = title_face),
        axis.text.x = element_text(size = axis_text_size, family = font_family),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = title_text_size, family = font_family, face = title_face), 
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    
    # Save the plot
    file_name <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/09.Cross-species_gene_comparison/05.Gene_level_final/", config$prefix, "dotplot_", name, "_genes", config$suffix, ".png")
    ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = 300, units = 'in')

    
    ########## with line
    plot_target <- ggplot(df_long, aes(x = log2FoldChange, y = Gene_Num, color = Species)) +
      geom_rect(aes(xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf), 
                fill = "grey96", inherit.aes = FALSE) +
      geom_vline(xintercept = c(-2.5, -1, 1, 2.5), 
                 linetype = "dashed", color = "gray", size = 0.5) +
      geom_segment(aes(x = 0, xend = log2FoldChange, 
                       y = Gene_Num, yend = Gene_Num), 
                   color = "gray70", size = 1.0) +
      geom_vline(xintercept = 0, color = "black", size = 0.8) +
      geom_point(size = 4, position = position_dodge(width = 0.6)) +
      # Add gray border for significant Mouse and Rat dots (padj < 0.05)
      geom_point(data = subset(df_long, is_significant & Species != "Human"),
                 aes(x = log2FoldChange, y = Gene_Num),
                 shape = 21, size = 4, stroke = 1.5, color = "gray50", alpha = 0.8, fill = NA,
                 position = position_dodge(width = 0.6)) +
      # Add a line that connects only the Human dots (red)
      geom_line(data = subset(df_long, Species == "Human"),
                aes(x = log2FoldChange, y = Gene_Num),
                color = "#d73027", size = 1) +
      
      scale_color_manual(values = species_color) +
      scale_y_continuous(
        breaks = df_long$Gene_Num, 
        labels = df_long$Gene) +
      labs(x = "Log2 Fold Change (Case vs. Control)") +
      theme_minimal() +
      theme(
        plot.title = element_blank(),
        axis.title = element_text(size = title_text_size, family = font_family, face = title_face),
        axis.text.x = element_text(size = axis_text_size, family = font_family),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = title_text_size, family = font_family, face = title_face), 
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    # Save the plot
    file_name <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/09.Cross-species_gene_comparison/05.Gene_level_final/", config$prefix, "dotplot_", name, "_genes_withline", config$suffix, ".png")
    ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = 300, units = 'in')
    
    file_name <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/09.Cross-species_gene_comparison/05.Gene_level_final/", config$prefix, "dotplot_", name, "_genes_withline", config$suffix, ".svg")
    ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = 300, units = 'in', device = 'svg')
    
  }  # End of target_config loop
  }  # End of config loop
