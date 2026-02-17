# Load required libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)  # To use rownames_to_column

log2fc_threshold_up <- 1
log2fc_threshold_down <- -1
padj_threshold <- 0.05

########## Human MTLE gene count ###############
# Import Human MTLE data
file_path <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/07.DESeq/01.DESeq_files/DESeq_FILTERED_1_MTLEALL_NL.txt"
res_df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
res_df <- rownames_to_column(res_df, var = "gene")

MTLE_up <- res_df %>%
  filter(log2FoldChange > log2fc_threshold_up & padj < padj_threshold) %>%
  pull(gene)
MTLE_down <- res_df %>%
  filter(log2FoldChange < log2fc_threshold_down & padj < padj_threshold) %>%
  pull(gene)

print(paste0("the number of elevated genes in epileptogenesis gene set is ", length(intersect(MTLE_up, H_gene_set_W))))
print(paste0("the number of decreased genes in epileptogenesis gene set is ", length(intersect(MTLE_down, H_gene_set_W))))

MTLE_up <- intersect(MTLE_up, H_gene_set_W)
MTLE_up <- as.character(MTLE_up)
output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/01.MTLE_up.txt"
writeLines(MTLE_up, output_file)
MTLE_down <- intersect(MTLE_down, H_gene_set_W)
MTLE_down <- as.character(MTLE_down)
output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/01.MTLE_down.txt"
writeLines(MTLE_down, output_file)

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

shared_genes_list_R_M <- list()  # Initialize a list to store the shared genes
for (group_name in names(files_and_labels)) {
  up_genes_var <- paste0(group_name, "_up_genes")
  down_genes_var <- paste0(group_name, "_down_genes")
  if (exists(up_genes_var) && exists(down_genes_var)) {
    up_genes_human <- get(up_genes_var)
    down_genes_human <- get(down_genes_var)
    shared_up_genes <- intersect(up_genes_human, H_gene_set_W)
    shared_down_genes <- intersect(down_genes_human, H_gene_set_W)
    shared_genes_list_R_M[[paste0(group_name, "_shared_up_genes")]] <- shared_up_genes
    shared_genes_list_R_M[[paste0(group_name, "_shared_down_genes")]] <- shared_down_genes
  }
}

for (name in names(shared_genes_list_R_M)) {
  cat("\nShared genes for:", name, "\n")
  print(shared_genes_list_R_M[[name]])
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

shared_genes_list_R <- list()  # Initialize a list to store the shared genes
for (group_name in names(files_and_labels)) {
  up_genes_var <- paste0(group_name, "_up_genes")
  down_genes_var <- paste0(group_name, "_down_genes")
  if (exists(up_genes_var) && exists(down_genes_var)) {
    up_genes_human <- get(up_genes_var)
    down_genes_human <- get(down_genes_var)
    shared_up_genes <- intersect(up_genes_human, H_gene_set_W)
    shared_down_genes <- intersect(down_genes_human, H_gene_set_W)
    shared_genes_list_R[[paste0(group_name, "_shared_up_genes")]] <- shared_up_genes
    shared_genes_list_R[[paste0(group_name, "_shared_down_genes")]] <- shared_down_genes
  }
}

for (name in names(shared_genes_list_R)) {
  cat("\nShared genes for:", name, "\n")
  print(shared_genes_list_R[[name]])
}


########################
# M_KAI_IH_IPSI_A
M_KAI_IH_IPSI_A_up <- union(union(shared_genes_list_R_M$M_KAI_IH_IPSI_A_HA_shared_up_genes, shared_genes_list_R_M$M_KAI_IH_IPSI_A_AC_shared_up_genes),
                       union(shared_genes_list_R_M$M_KAI_IH_IPSI_A_IM_shared_up_genes, shared_genes_list_R_M$M_KAI_IH_IPSI_A_CR_shared_up_genes))
M_KAI_IH_IPSI_A_up <- as.character(M_KAI_IH_IPSI_A_up)
output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/01.M_KAI_IH_IPSI_A_up.txt"
writeLines(M_KAI_IH_IPSI_A_up, output_file)
M_KAI_IH_IPSI_A_down <- union(union(shared_genes_list_R_M$M_KAI_IH_IPSI_A_HA_shared_down_genes, shared_genes_list_R_M$M_KAI_IH_IPSI_A_AC_shared_down_genes),
                       union(shared_genes_list_R_M$M_KAI_IH_IPSI_A_IM_shared_down_genes, shared_genes_list_R_M$M_KAI_IH_IPSI_A_CR_shared_down_genes))
M_KAI_IH_IPSI_A_down <- as.character(M_KAI_IH_IPSI_A_down)
output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/01.M_KAI_IH_IPSI_A_down.txt"
writeLines(M_KAI_IH_IPSI_A_down, output_file)

# M_PILO_IP_A
M_PILO_IP_A_up <- union(union(shared_genes_list_R_M$M_PILO_IP_A_HA_shared_up_genes, shared_genes_list_R_M$M_PILO_IP_A_AC_shared_up_genes),
                       shared_genes_list_R_M$M_PILO_IP_A_IM_shared_up_genes)
M_PILO_IP_A_up <- as.character(M_PILO_IP_A_up)
output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/01.M_PILO_IP_A_up.txt"
writeLines(M_PILO_IP_A_up, output_file)
M_PILO_IP_A_down <- union(union(shared_genes_list_R_M$M_PILO_IP_A_HA_shared_down_genes, shared_genes_list_R_M$M_PILO_IP_A_AC_shared_down_genes),
                       shared_genes_list_R_M$M_PILO_IP_A_IM_shared_down_genes)
M_PILO_IP_A_down <- as.character(M_PILO_IP_A_down)
output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/01.M_PILO_IP_A_down.txt"
writeLines(M_PILO_IP_A_down, output_file)

# R_PPS_A
R_PPS_A_up <- union(union(shared_genes_list_R$R_PPS_A_AC_shared_up_genes, shared_genes_list_R$R_PPS_A_DOFS_shared_up_genes),
                       union(shared_genes_list_R$R_PPS_A_IM_shared_up_genes, shared_genes_list_R$R_PPS_A_CR_shared_up_genes))
R_PPS_A_up <- as.character(R_PPS_A_up)
output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/01.R_PPS_A_up.txt"
writeLines(R_PPS_A_up, output_file)
R_PPS_A_down <- union(union(shared_genes_list_R$R_PPS_A_AC_shared_down_genes, shared_genes_list_R$R_PPS_A_DOFS_shared_down_genes),
                       union(shared_genes_list_R$R_PPS_A_IM_shared_down_genes, shared_genes_list_R$R_PPS_A_CR_shared_down_genes))
R_PPS_A_down <- as.character(R_PPS_A_down)
output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/01.R_PPS_A_down.txt"
writeLines(R_PPS_A_down, output_file)


## up (intersection & unique)
intersection_M_KAI_IH_IPSI_R_PPS_A_up <- intersect(M_KAI_IH_IPSI_A_up, R_PPS_A_up)
intersection_M_KAI_IH_IPSI_R_PPS_A_up <- as.character(intersection_M_KAI_IH_IPSI_R_PPS_A_up)
unique_M_KAI_IH_IPSI_A_up <- setdiff(M_KAI_IH_IPSI_A_up, R_PPS_A_up)
unique_M_KAI_IH_IPSI_A_up <- as.character(unique_M_KAI_IH_IPSI_A_up)
unique_R_PPS_A_up <- setdiff(R_PPS_A_up, M_KAI_IH_IPSI_A_up)
unique_R_PPS_A_up <- as.character(unique_R_PPS_A_up)

output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/02.intersection_M_KAI_IH_IPSI_R_PPS_A_up.txt"
writeLines(intersection_M_KAI_IH_IPSI_R_PPS_A_up, output_file)
output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/02.unique_M_KAI_IH_IPSI_A_up.txt"
writeLines(unique_M_KAI_IH_IPSI_A_up, output_file)
output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/02.unique_R_PPS_A_up.txt"
writeLines(unique_R_PPS_A_up, output_file)

## down (intersection & unique)
intersection_M_KAI_IH_IPSI_R_PPS_A_down <- intersect(M_KAI_IH_IPSI_A_down, R_PPS_A_down)
intersection_M_KAI_IH_IPSI_R_PPS_A_down <- as.character(intersection_M_KAI_IH_IPSI_R_PPS_A_down)
unique_M_KAI_IH_IPSI_A_down <- setdiff(M_KAI_IH_IPSI_A_down, R_PPS_A_down)
unique_M_KAI_IH_IPSI_A_down <- as.character(unique_M_KAI_IH_IPSI_A_down)
unique_R_PPS_A_down <- setdiff(R_PPS_A_down, M_KAI_IH_IPSI_A_down)
unique_R_PPS_A_down <- as.character(unique_R_PPS_A_down)

output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/02.intersection_M_KAI_IH_IPSI_R_PPS_A_down.txt"
writeLines(intersection_M_KAI_IH_IPSI_R_PPS_A_down, output_file)
output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/02.unique_M_KAI_IH_IPSI_A_down.txt"
writeLines(unique_M_KAI_IH_IPSI_A_down, output_file)
output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/02.unique_R_PPS_A_down.txt"
writeLines(unique_R_PPS_A_down, output_file)


########################
## M_KAI_IH_IPSI_A and MTLE (up - intersection & unique)
if (length(MTLE_up) > 0) {
  intersection_M_KAI_IH_IPSI_MTLE_up <- intersect(M_KAI_IH_IPSI_A_up, MTLE_up)
  intersection_M_KAI_IH_IPSI_MTLE_up <- as.character(intersection_M_KAI_IH_IPSI_MTLE_up)
  unique_M_KAI_IH_IPSI_A_MTLE_up <- setdiff(M_KAI_IH_IPSI_A_up, MTLE_up)
  unique_M_KAI_IH_IPSI_A_MTLE_up <- as.character(unique_M_KAI_IH_IPSI_A_MTLE_up)
  unique_MTLE_M_KAI_IH_IPSI_up <- setdiff(MTLE_up, M_KAI_IH_IPSI_A_up)
  unique_MTLE_M_KAI_IH_IPSI_up <- as.character(unique_MTLE_M_KAI_IH_IPSI_up)
  
  output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/03.intersection_M_KAI_IH_IPSI_MTLE_up.txt"
  writeLines(intersection_M_KAI_IH_IPSI_MTLE_up, output_file)
  output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/03.unique_M_KAI_IH_IPSI_A_MTLE_up.txt"
  writeLines(unique_M_KAI_IH_IPSI_A_MTLE_up, output_file)
  output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/03.unique_MTLE_M_KAI_IH_IPSI_up.txt"
  writeLines(unique_MTLE_M_KAI_IH_IPSI_up, output_file)
}

## M_KAI_IH_IPSI_A and MTLE (down - intersection & unique)
if (length(MTLE_down) > 0) {
  intersection_M_KAI_IH_IPSI_MTLE_down <- intersect(M_KAI_IH_IPSI_A_down, MTLE_down)
  intersection_M_KAI_IH_IPSI_MTLE_down <- as.character(intersection_M_KAI_IH_IPSI_MTLE_down)
  unique_M_KAI_IH_IPSI_A_MTLE_down <- setdiff(M_KAI_IH_IPSI_A_down, MTLE_down)
  unique_M_KAI_IH_IPSI_A_MTLE_down <- as.character(unique_M_KAI_IH_IPSI_A_MTLE_down)
  unique_MTLE_M_KAI_IH_IPSI_down <- setdiff(MTLE_down, M_KAI_IH_IPSI_A_down)
  unique_MTLE_M_KAI_IH_IPSI_down <- as.character(unique_MTLE_M_KAI_IH_IPSI_down)
  
  output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/03.intersection_M_KAI_IH_IPSI_MTLE_down.txt"
  writeLines(intersection_M_KAI_IH_IPSI_MTLE_down, output_file)
  output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/03.unique_M_KAI_IH_IPSI_A_MTLE_down.txt"
  writeLines(unique_M_KAI_IH_IPSI_A_MTLE_down, output_file)
  output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/03.unique_MTLE_M_KAI_IH_IPSI_down.txt"
  writeLines(unique_MTLE_M_KAI_IH_IPSI_down, output_file)
}

## R_PPS_A and MTLE (up - intersection & unique)
if (length(MTLE_up) > 0) {
  intersection_R_PPS_MTLE_up <- intersect(R_PPS_A_up, MTLE_up)
  intersection_R_PPS_MTLE_up <- as.character(intersection_R_PPS_MTLE_up)
  unique_R_PPS_A_MTLE_up <- setdiff(R_PPS_A_up, MTLE_up)
  unique_R_PPS_A_MTLE_up <- as.character(unique_R_PPS_A_MTLE_up)
  unique_MTLE_R_PPS_up <- setdiff(MTLE_up, R_PPS_A_up)
  unique_MTLE_R_PPS_up <- as.character(unique_MTLE_R_PPS_up)
  
  output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/04.intersection_R_PPS_MTLE_up.txt"
  writeLines(intersection_R_PPS_MTLE_up, output_file)
  output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/04.unique_R_PPS_A_MTLE_up.txt"
  writeLines(unique_R_PPS_A_MTLE_up, output_file)
  output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/04.unique_MTLE_R_PPS_up.txt"
  writeLines(unique_MTLE_R_PPS_up, output_file)
}

## R_PPS_A and MTLE (down - intersection & unique)
if (length(MTLE_down) > 0) {
  intersection_R_PPS_MTLE_down <- intersect(R_PPS_A_down, MTLE_down)
  intersection_R_PPS_MTLE_down <- as.character(intersection_R_PPS_MTLE_down)
  unique_R_PPS_A_MTLE_down <- setdiff(R_PPS_A_down, MTLE_down)
  unique_R_PPS_A_MTLE_down <- as.character(unique_R_PPS_A_MTLE_down)
  unique_MTLE_R_PPS_down <- setdiff(MTLE_down, R_PPS_A_down)
  unique_MTLE_R_PPS_down <- as.character(unique_MTLE_R_PPS_down)
  
  output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/04.intersection_R_PPS_MTLE_down.txt"
  writeLines(intersection_R_PPS_MTLE_down, output_file)
  output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/04.unique_R_PPS_A_MTLE_down.txt"
  writeLines(unique_R_PPS_A_MTLE_down, output_file)
  output_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_01.Matrix/04.unique_MTLE_R_PPS_down.txt"
  writeLines(unique_MTLE_R_PPS_down, output_file)
}

