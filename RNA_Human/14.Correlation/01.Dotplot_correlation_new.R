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
file_path <- paste0("/home/joonho345/3_RNA/RNA_Human/06.DESeq/01.DESeq_files/DESeq_", Target_comparison, ".txt")
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
csv_file <- "/home/joonho345/3_RNA/RNA_Animal/Raw_Data/Animal_sheet_FTP_2.csv"
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
    intersections[[paste0("MTLE_", group_name, "_up")]] <- intersect(MTLE_up_genes, up_genes)
    intersections[[paste0("MTLE_", group_name, "_down")]] <- intersect(MTLE_down_genes, down_genes)
  }
}
print(intersections)


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
csv_file <- "/home/joonho345/3_RNA/RNA_Animal/Raw_Data/Animal_sheet_FTP_2.csv"
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
    intersections[[paste0("MTLE_", group_name, "_up")]] <- intersect(MTLE_up_genes, up_genes)
    intersections[[paste0("MTLE_", group_name, "_down")]] <- intersect(MTLE_down_genes, down_genes)
  }
}
print(intersections)


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

print(transformed_df)


############ Shared gene counts + Correlation #############
input_file <- paste0("/home/joonho345/3_RNA/RNA_Human/10.Correlation/00.Correlation_Matrix_Mouse_FILTERED.txt")
correlation_matrix_SAMPLE_M <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
input_file <- paste0("/home/joonho345/3_RNA/RNA_Human/10.Correlation/00.Correlation_Matrix_Rat_FILTERED.txt")
correlation_matrix_SAMPLE_R <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

colnames(correlation_matrix_SAMPLE_M) <- gsub("_", ".", colnames(correlation_matrix_SAMPLE_M))
colnames(correlation_matrix_SAMPLE_R) <- gsub("_", ".", colnames(correlation_matrix_SAMPLE_R))
correlation_data_M <- as.numeric(correlation_matrix_SAMPLE_M["MTLEALL vs. NL", ])
correlation_data_R <- as.numeric(correlation_matrix_SAMPLE_R["MTLEALL vs. NL", ])
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
  "M-KAI-ST-BO-O-HA", "M-KAI-ST-BO-O-AC", "M-KAI-ST-BO-O-IM", "M-KAI-ST-BO-O-CR",
  "M-KAI-ST-IPSI-Y-AC", "M-KAI-ST-IPSI-Y-IM", "M-KAI-ST-IPSI-Y-CR", 
  "M-KAI-ST-CON-Y-AC", "M-KAI-ST-CON-Y-IM", "M-KAI-ST-CON-Y-CR", 
  "M-KAI-IP-O-HA", "R-KAI-IP-Y-CR", "R-KAI-SUB-I-IM", 
  "M-PILO-IP-Y-HA", "M-PILO-IP-Y-AC", "M-PILO-IP-Y-IM", 
  "M-PILO-IP-O-HA", "M-PILO-IP-O-IM", "M-PILO-IP-O-CR", 
  "R-PILO-IP-Y-CR", "R-PPS-Y-AC", "R-PPS-Y-DOFS", "R-PPS-Y-IM", "R-PPS-Y-CR", 
  "R-AMG-IPSI-Y-CR", "R-TBI-IPSI-Y-CR"
)

transformed_df <- transformed_df %>%
  mutate(Group = factor(Group, levels = group_order)) %>%
  arrange(Group)
print(transformed_df)




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

plot_target <- ggplot(transformed_df, aes(x = Group, y = Correlation)) +
  geom_hline(yintercept = c(0.1, 0.2, 0.3), linetype = "dashed", color = "gray", size = 0.3, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray", size = 0.5, alpha = 0.3) +
  geom_vline(xintercept = 20.5, linetype = "dashed", color = "black", size = 0.5, alpha = 0.5) + 
  geom_vline(xintercept = 13.5, linetype = "dashed", color = "gray", size = 0.5, alpha = 0.5) + 
  geom_vline(xintercept = 25.5, linetype = "dashed", color = "gray", size = 0.5, alpha = 0.5) + 
  
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
file_name <- "/home/joonho345/3_RNA/RNA_Human/14.DEG_shared_count/01.Dotplot_correlation_new.png"
ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)


### without tick ###
axis_text_size <- 10
legend_text_size <- 8
title_text_size <- 10
title_face <- "bold"
font_family <- "Arial"
plot_width <- 13
plot_height <- 1.5
plot_unit <- 'in'
plot_dpi <- 300
plot_target <- ggplot(transformed_df, aes(x = Group, y = Correlation)) +
  geom_hline(yintercept = c(0.1, 0.2, 0.3), linetype = "dashed", color = "gray", size = 0.3, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray", size = 0.5, alpha = 0.3) +
  geom_vline(xintercept = 20.5, linetype = "dashed", color = "black", size = 0.5, alpha = 0.5) + 
  geom_vline(xintercept = 13.5, linetype = "dashed", color = "gray", size = 0.5, alpha = 0.5) + 
  geom_vline(xintercept = 25.5, linetype = "dashed", color = "gray", size = 0.5, alpha = 0.5) + 

  geom_line(aes(group = 1), size = 1, alpha = 0.5, color = "#d73027") +
  geom_point(size = 3, alpha = 0.8, color = "#d73027") +
  
  theme_minimal() +
  theme(
    plot.title = element_blank(),
    axis.title =  element_text(size = title_text_size, family = font_family),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = axis_text_size, family = font_family), 
    legend.title = element_blank(),
    legend.text = element_blank(), 
    legend.position = "none",
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA), 
  )

file_name <- "/home/joonho345/3_RNA/RNA_Human/14.DEG_shared_count/01.Dotplot_correlation_new_tick.png"
ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)

# without tick and line
plot_target <- ggplot(transformed_df, aes(x = Group, y = Correlation)) +
  geom_hline(yintercept = c(0, 0.1, 0.2, 0.3), linetype = "dashed", color = "gray", size = 0.3, alpha = 0.3) +
  #geom_hline(yintercept = 0, linetype = "dashed", color = "gray", size = 0.5, alpha = 0.3) +
  #geom_vline(xintercept = 20.5, linetype = "dashed", color = "black", size = 0.5, alpha = 0.5) + 
  #geom_vline(xintercept = 13.5, linetype = "dashed", color = "gray", size = 0.5, alpha = 0.5) + 
  #geom_vline(xintercept = 25.5, linetype = "dashed", color = "gray", size = 0.5, alpha = 0.5) + 

  geom_line(aes(group = 1), size = 1, alpha = 0.5, color = "#d73027") +
  geom_point(size = 3, alpha = 0.8, color = "#d73027") +

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
    plot.background = element_rect(fill = "white", color = NA), 
  )

file_name <- "/home/joonho345/3_RNA/RNA_Human/14.DEG_shared_count/01.Dotplot_correlation_new_tick_line.png"
ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)

file_name <- "/home/joonho345/3_RNA/RNA_Human/14.DEG_shared_count/01.Dotplot_correlation_new_tick_line.svg"
ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit, device = 'svg')

