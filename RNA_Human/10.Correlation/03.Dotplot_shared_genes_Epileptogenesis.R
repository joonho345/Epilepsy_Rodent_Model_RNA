# Load required libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)  # To use rownames_to_column


########## Human gene count ##########
# Define Target Comparison
Target_comparison <- 'FILTERED_1_MTLEALL_NL'
Target_group_name <- '1'

# File path for the DESeq results
file_path <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/07.DESeq/01.DESeq_files/DESeq_", Target_comparison, ".txt")
res_df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
res_df <- rownames_to_column(res_df, var = "gene")

log2fc_threshold_up <- 1
log2fc_threshold_down <- -1
padj_threshold <- 0.05

MTLE_up_genes <- res_df %>%
  filter(log2FoldChange > log2fc_threshold_up & padj < padj_threshold) %>%
  pull(gene)
MTLE_down_genes <- res_df %>%
  filter(log2FoldChange < log2fc_threshold_down & padj < padj_threshold) %>%
  pull(gene)

# Filter to only H_gene_set_W genes (epileptogenesis GO gene sets)
MTLE_up_genes <- MTLE_up_genes[MTLE_up_genes %in% H_gene_set_W]
MTLE_down_genes <- MTLE_down_genes[MTLE_down_genes %in% H_gene_set_W]


########## Mouse gene count ###############
# Import Mouse data
Species <- 'Mouse'
orthologs_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_mm_17128.txt"
ortholog_one2one <- read.table(orthologs_path, sep = "\t", header = TRUE, fill = TRUE)
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one$Gene.name), ]
Species_Gene_name <- paste0(Species, '.gene.name')
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one[[Species_Gene_name]]), ]
mouse_to_human_map <- setNames(ortholog_one2one$Gene.name, ortholog_one2one[[Species_Gene_name]])

# Import 
csv_file <- "/data/project/1_Epilepsy_RNA/RNA_Animal/Raw_Data/Animal_sheet_FTP_3.csv"
file_data <- read.csv(csv_file, stringsAsFactors = FALSE)
file_data <- file_data %>% filter(Species == "M")
files_and_labels <- setNames(file_data$Filepath, file_data$Group)

for (group_name in names(files_and_labels)) {
  file_path <- files_and_labels[group_name]
  gene_data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  gene_data <- rownames_to_column(gene_data, var = "gene")  # Ensure 'gene' is a column
  up_genes <- gene_data %>%
    filter(log2FoldChange > log2fc_threshold_up & padj < padj_threshold) %>%
    pull(gene)
  down_genes <- gene_data %>%
    filter(log2FoldChange < log2fc_threshold_down & padj < padj_threshold) %>%
    pull(gene)
  up_genes_human <- mouse_to_human_map[up_genes]
  down_genes_human <- mouse_to_human_map[down_genes]
  up_genes_human <- as.character(na.omit(up_genes_human))
  down_genes_human <- as.character(na.omit(down_genes_human))
  assign(paste0(group_name, "_up_genes"), up_genes_human, envir = .GlobalEnv)
  assign(paste0(group_name, "_down_genes"), down_genes_human, envir = .GlobalEnv)
}

intersections <- list()
for (group_name in names(files_and_labels)) {
  if (exists(paste0(group_name, "_up_genes")) && exists(paste0(group_name, "_down_genes"))) {
    up_genes <- get(paste0(group_name, "_up_genes"))
    down_genes <- get(paste0(group_name, "_down_genes"))
    # Filter to only H_gene_set_W genes
    up_genes <- up_genes[up_genes %in% H_gene_set_W]
    down_genes <- down_genes[down_genes %in% H_gene_set_W]
    intersections[[paste0("MTLE_", group_name, "_up")]] <- intersect(MTLE_up_genes, up_genes)
    intersections[[paste0("MTLE_", group_name, "_down")]] <- intersect(MTLE_down_genes, down_genes)
  }
}


########## Rat gene count ###############
# Import Rat data
Species <- 'Rat'
orthologs_path <- "/home/joonho345/resources/GeneSet/Orthologs/One2one_hs_rn_16211.txt"
ortholog_one2one <- read.table(orthologs_path, sep = "\t", header = TRUE, fill = TRUE)
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one$Gene.name), ]
Species_Gene_name <- paste0(Species, '.gene.name')
ortholog_one2one <- ortholog_one2one[!duplicated(ortholog_one2one[[Species_Gene_name]]), ]
rat_to_human_map <- setNames(ortholog_one2one$Gene.name, ortholog_one2one[[Species_Gene_name]])

# Import 
csv_file <- "/data/project/1_Epilepsy_RNA/RNA_Animal/Raw_Data/Animal_sheet_FTP_3.csv"
file_data <- read.csv(csv_file, stringsAsFactors = FALSE)
file_data <- file_data %>% filter(Species == "R")
files_and_labels <- setNames(file_data$Filepath, file_data$Group)

for (group_name in names(files_and_labels)) {
  file_path <- files_and_labels[group_name]
  gene_data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  gene_data <- rownames_to_column(gene_data, var = "gene")  # Ensure 'gene' is a column
  up_genes <- gene_data %>%
    filter(log2FoldChange > log2fc_threshold_up & padj < padj_threshold) %>%
    pull(gene)
  down_genes <- gene_data %>%
    filter(log2FoldChange < log2fc_threshold_down & padj < padj_threshold) %>%
    pull(gene)
  up_genes_human <- rat_to_human_map[up_genes]
  down_genes_human <- rat_to_human_map[down_genes]
  up_genes_human <- as.character(na.omit(up_genes_human))
  down_genes_human <- as.character(na.omit(down_genes_human))
  assign(paste0(group_name, "_up_genes"), up_genes_human, envir = .GlobalEnv)
  assign(paste0(group_name, "_down_genes"), down_genes_human, envir = .GlobalEnv)
}

for (group_name in names(files_and_labels)) {
  if (exists(paste0(group_name, "_up_genes")) && exists(paste0(group_name, "_down_genes"))) {
    up_genes <- get(paste0(group_name, "_up_genes"))
    down_genes <- get(paste0(group_name, "_down_genes"))
    # Filter to only H_gene_set_W genes
    up_genes <- up_genes[up_genes %in% H_gene_set_W]
    down_genes <- down_genes[down_genes %in% H_gene_set_W]
    intersections[[paste0("MTLE_", group_name, "_up")]] <- intersect(MTLE_up_genes, up_genes)
    intersections[[paste0("MTLE_", group_name, "_down")]] <- intersect(MTLE_down_genes, down_genes)
  }
}


############ Shared gene counts #############
intersection_df <- data.frame(
  Vector_Name = character(),
  Gene_Count = integer(),
  stringsAsFactors = FALSE
)

for (name in names(intersections)) {
  gene_count <- length(intersections[[name]])
  intersection_df <- rbind(intersection_df, data.frame(Vector_Name = name, Gene_Count = gene_count))
}
intersection_df <- intersection_df %>%
  mutate(Type = ifelse(grepl("_up$", Vector_Name), "up", "down"))

transformed_df <- intersection_df %>%
  mutate(Group = gsub("_", ".", gsub("MTLE_", "", gsub("_(up|down)$", "", Vector_Name)))) %>%
  group_by(Group) %>%
  summarise(
    Gene_Count_up = sum(Gene_Count[Type == "up"]),
    Gene_Count_down = sum(Gene_Count[Type == "down"])
  ) %>%
  arrange(Group)


############ Shared gene counts + Correlation #############
input_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/01.Correlation/02.Correlation_Epileptogenesis_Mouse_FILTERED.txt")
correlation_matrix_SAMPLE_M <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
input_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/01.Correlation/02.Correlation_Epileptogenesis_Rat_FILTERED.txt")
correlation_matrix_SAMPLE_R <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

colnames(correlation_matrix_SAMPLE_M) <- gsub("_", ".", colnames(correlation_matrix_SAMPLE_M))
colnames(correlation_matrix_SAMPLE_R) <- gsub("_", ".", colnames(correlation_matrix_SAMPLE_R))
correlation_data_M <- as.numeric(correlation_matrix_SAMPLE_M[Target_comparison, ])
correlation_data_R <- as.numeric(correlation_matrix_SAMPLE_R[Target_comparison, ])
correlation_map_M <- setNames(correlation_data_M, colnames(correlation_matrix_SAMPLE_M))
correlation_map_R <- setNames(correlation_data_R, colnames(correlation_matrix_SAMPLE_R))

transformed_df <- transformed_df %>%
  mutate(
    Correlation = ifelse(
      startsWith(Group, "M"),
      correlation_map_M[Group],  # Use Mouse correlations for groups starting with M
      correlation_map_R[Group]   # Use Rat correlations for groups starting with R
    )
  )
transformed_df <- transformed_df %>%
  mutate(Group = gsub("\\.", "-", Group))

group_order <- c(
  "M-KAI-IH-IPSI-A-HA", "M-KAI-IH-IPSI-A-AC", "M-KAI-IH-IPSI-A-IM", "M-KAI-IH-IPSI-A-CR",
  "M-KAI-IH-CON-A-AC", "M-KAI-IH-CON-A-IM", "M-KAI-IH-CON-A-CR", 
  "M-KAI-IA-IPSI-A-AC", "M-KAI-IA-IPSI-A-IM", "M-KAI-IA-IPSI-A-CR",
  "M-KAI-IP-A-HA", "R-KAI-IP-A-CR", "R-KAI-SUB-I-IM", 
  "M-PILO-IP-A-HA", "M-PILO-IP-A-AC", "M-PILO-IP-A-IM", "M-PILO-IP-A-CR", 
  "R-PILO-IP-A-CR", "R-PPS-A-AC", "R-PPS-A-DOFS", "R-PPS-A-IM", "R-PPS-A-CR", 
  "R-AMG-IPSI-A-CR", "R-TBI-IPSI-A-CR"
)

# Filter to only include groups in group_order and set factor levels
transformed_df <- transformed_df %>%
  filter(Group %in% group_order) %>%
  mutate(Group = factor(Group, levels = group_order)) %>%
  arrange(Group)
print(transformed_df, n = Inf)

###########
transformed_df_long <- transformed_df %>%
  dplyr::select(Group, Gene_Count_up, Gene_Count_down) %>%
  tidyr::pivot_longer(cols = c(Gene_Count_up, Gene_Count_down), 
                      names_to = "Gene_Type", 
                      values_to = "Gene_Count") %>%
  mutate(Group = factor(Group, levels = group_order))

gene_colors <- c("Gene_Count_up" = "#313695", "Gene_Count_down" = "#1a9850")

# Create filtered dataframes for lines (multiple groups)
# Group 1: models 1-7
line_group_1 <- group_order[1:7]
transformed_df_line_1_up <- transformed_df_long %>%
  filter(Group %in% line_group_1 & Gene_Type == "Gene_Count_up") %>%
  mutate(Group = factor(Group, levels = group_order))
transformed_df_line_1_down <- transformed_df_long %>%
  filter(Group %in% line_group_1 & Gene_Type == "Gene_Count_down") %>%
  mutate(Group = factor(Group, levels = group_order))

# Group 2: models 8-10
line_group_2 <- group_order[8:10]
transformed_df_line_2_up <- transformed_df_long %>%
  filter(Group %in% line_group_2 & Gene_Type == "Gene_Count_up") %>%
  mutate(Group = factor(Group, levels = group_order))
transformed_df_line_2_down <- transformed_df_long %>%
  filter(Group %in% line_group_2 & Gene_Type == "Gene_Count_down") %>%
  mutate(Group = factor(Group, levels = group_order))

# Group 3: models 11-12
line_group_3 <- group_order[11:12]
transformed_df_line_3_up <- transformed_df_long %>%
  filter(Group %in% line_group_3 & Gene_Type == "Gene_Count_up") %>%
  mutate(Group = factor(Group, levels = group_order))
transformed_df_line_3_down <- transformed_df_long %>%
  filter(Group %in% line_group_3 & Gene_Type == "Gene_Count_down") %>%
  mutate(Group = factor(Group, levels = group_order))

# Group 4: models 14-18 (skipping 13)
line_group_4 <- group_order[14:18]
transformed_df_line_4_up <- transformed_df_long %>%
  filter(Group %in% line_group_4 & Gene_Type == "Gene_Count_up") %>%
  mutate(Group = factor(Group, levels = group_order))
transformed_df_line_4_down <- transformed_df_long %>%
  filter(Group %in% line_group_4 & Gene_Type == "Gene_Count_down") %>%
  mutate(Group = factor(Group, levels = group_order))

# Group 5: models 19-22
line_group_5 <- group_order[19:22]
transformed_df_line_5_up <- transformed_df_long %>%
  filter(Group %in% line_group_5 & Gene_Type == "Gene_Count_up") %>%
  mutate(Group = factor(Group, levels = group_order))
transformed_df_line_5_down <- transformed_df_long %>%
  filter(Group %in% line_group_5 & Gene_Type == "Gene_Count_down") %>%
  mutate(Group = factor(Group, levels = group_order))

############## plot ###########
# Create the dot plot
axis_text_size <- 10
legend_text_size <- 8
title_text_size <- 12
title_face <- "bold"
font_family <- "Arial"
plot_width <- 13
plot_height <- 4
plot_unit <- 'in'
plot_dpi <- 300

plot_target <- ggplot(transformed_df_long, aes(x = Group, y = Gene_Count, color = Gene_Type, group = Gene_Type)) +
  geom_hline(yintercept = c(25, 50, 75, 100), linetype = "dashed", color = "gray", size = 0.3, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray", size = 0.5, alpha = 0.3) +
  geom_vline(xintercept = 20.5, linetype = "dashed", color = "black", size = 0.5, alpha = 0.5) + 
  geom_vline(xintercept = 13.5, linetype = "dashed", color = "gray", size = 0.5, alpha = 0.5) + 
  geom_vline(xintercept = 25.5, linetype = "dashed", color = "gray", size = 0.5, alpha = 0.5) + 
  
  # Up-regulated genes lines
  geom_line(data = transformed_df_line_1_up, aes(x = Group, y = Gene_Count, group = 1), size = 1, alpha = 0.5, color = "#313695") +
  geom_line(data = transformed_df_line_2_up, aes(x = Group, y = Gene_Count, group = 1), size = 1, alpha = 0.5, color = "#313695") +
  geom_line(data = transformed_df_line_3_up, aes(x = Group, y = Gene_Count, group = 1), size = 1, alpha = 0.5, color = "#313695") +
  geom_line(data = transformed_df_line_4_up, aes(x = Group, y = Gene_Count, group = 1), size = 1, alpha = 0.5, color = "#313695") +
  geom_line(data = transformed_df_line_5_up, aes(x = Group, y = Gene_Count, group = 1), size = 1, alpha = 0.5, color = "#313695") +
  # Down-regulated genes lines
  geom_line(data = transformed_df_line_1_down, aes(x = Group, y = Gene_Count, group = 1), size = 1, alpha = 0.5, color = "#1a9850") +
  geom_line(data = transformed_df_line_2_down, aes(x = Group, y = Gene_Count, group = 1), size = 1, alpha = 0.5, color = "#1a9850") +
  geom_line(data = transformed_df_line_3_down, aes(x = Group, y = Gene_Count, group = 1), size = 1, alpha = 0.5, color = "#1a9850") +
  geom_line(data = transformed_df_line_4_down, aes(x = Group, y = Gene_Count, group = 1), size = 1, alpha = 0.5, color = "#1a9850") +
  geom_line(data = transformed_df_line_5_down, aes(x = Group, y = Gene_Count, group = 1), size = 1, alpha = 0.5, color = "#1a9850") +
  geom_point(size = 3, alpha = 1.0) +
  scale_x_discrete(limits = group_order) +
  scale_color_manual(values = gene_colors, labels = c("Up-regulated Genes", "Down-regulated Genes")) +
  
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
  ) +
  labs(x = "Group", y = "Gene Count", color = "Gene Type")

file_name <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/03.Dotplot_shared_genes/02.Dotplot_shared_genes_Epileptogenesis_", Target_group_name, ".png")
ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)
