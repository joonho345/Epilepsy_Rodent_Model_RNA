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

log2fc_threshold_up <- 0.5
log2fc_threshold_down <- -0.5
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

############### txt file
genes_to_save <- union(intersections$MTLE_M_KAI_ST_BO_O_IM_up, intersections$MTLE_M_KAI_ST_BO_O_CR_up)
output_file <- "/home/joonho345/3_RNA/RNA_Human/14.DEG_shared_count/02.Shared_Genes_M_KAI_ST_BO_05.txt"
writeLines(genes_to_save, output_file)

genes_to_save <- union(intersections$MTLE_M_PILO_IP_Y_AC_up, intersections$MTLE_M_PILO_IP_Y_IM_up)
output_file <- "/home/joonho345/3_RNA/RNA_Human/14.DEG_shared_count/02.Shared_Genes_M_PILO_Y_05.txt"
writeLines(genes_to_save, output_file)

genes_to_save <- union(intersections$MTLE_R_PPS_Y_DOFS_up, intersections$MTLE_R_PPS_Y_CR_up)
output_file <- "/home/joonho345/3_RNA/RNA_Human/14.DEG_shared_count/02.Shared_Genes_R_PPS_O_05.txt"
writeLines(genes_to_save, output_file)


######## Epileptogenesis ########
# Function to filter genes by H_gene_set_W, print them, and save to file
filter_and_save_genes <- function(genes, gene_set, output_file) {
  filtered_genes <- genes[genes %in% gene_set]  # Filter genes present in H_gene_set_W
  cat("Filtered Genes to be Saved (", output_file, "):\n", sep = "")
  print(filtered_genes)  # Print the filtered genes
  writeLines(filtered_genes, output_file)  # Save filtered genes to the specified file
}

genes_to_save <- union(intersections$MTLE_M_KAI_ST_BO_O_IM_up, intersections$MTLE_M_KAI_ST_BO_O_CR_up)
output_file <- "/home/joonho345/3_RNA/RNA_Human/14.DEG_shared_count/02.Shared_Genes_M_KAI_ST_BO_05_filtered.txt"
filter_and_save_genes(genes_to_save, H_gene_set_W, output_file)

genes_to_save <- union(intersections$MTLE_M_PILO_IP_Y_AC_up, intersections$MTLE_M_PILO_IP_Y_IM_up)
output_file <- "/home/joonho345/3_RNA/RNA_Human/14.DEG_shared_count/02.Shared_Genes_M_PILO_Y_05_filtered.txt"
filter_and_save_genes(genes_to_save, H_gene_set_W, output_file)

genes_to_save <- union(intersections$MTLE_R_PPS_Y_DOFS_up, intersections$MTLE_R_PPS_Y_CR_up)
output_file <- "/home/joonho345/3_RNA/RNA_Human/14.DEG_shared_count/02.Shared_Genes_R_PPS_O_05_filtered.txt"
filter_and_save_genes(genes_to_save, H_gene_set_W, output_file)

