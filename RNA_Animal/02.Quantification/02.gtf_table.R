library(dplyr)
library(tidyr)

gtf_path_M <- "/data/resource/reference/mouse/Mus_musculus.GRCm39.112.gtf"
gtf_table_M <- read.table(gtf_path_M, header = FALSE, skip = 5, sep = "\t", quote = "", comment.char = "")
gtf_table_M <- gtf_table_M %>% dplyr::select(-c(V1, V4, V5, V6, V7, V8))
gtf_table_M <- gtf_table_M[gtf_table_M$V3 == "gene", ]

gtf_path_R <- "/data/resource/reference/mouse/Rattus_norvegicus.mRatBN7.2.112.gtf"
gtf_table_R <- read.table(gtf_path_R, header = FALSE, skip = 5, sep = "\t", quote = "", comment.char = "")
gtf_table_R <- gtf_table_R %>% dplyr::select(-c(V1, V4, V5, V6, V7, V8))
gtf_table_R <- gtf_table_R[gtf_table_R$V3 == "gene", ]


### Create gene_id column
extract_gene_id <- function(v9_string) {
  # Use a regular expression to find and extract the gene_id value
  matches <- regmatches(v9_string, regexpr('gene_id "\\S+"', v9_string))
  gene_id <- ifelse(length(matches) > 0, sub('gene_id "(\\S+)"', '\\1', matches), NA)
  return(gene_id)
}
gene_ids_M <- sapply(gtf_table_M$V9, extract_gene_id)
gtf_table_M$gene_ids_M <- gene_ids_M
gene_ids_R <- sapply(gtf_table_R$V9, extract_gene_id)
gtf_table_R$gene_ids_R <- gene_ids_R


### Create gene_name column
extract_gene_name <- function(v9_string) {
  # Use a regular expression to find and extract the gene_id value
  matches <- regmatches(v9_string, regexpr('gene_name "\\S+"', v9_string))
  gene_name <- ifelse(length(matches) > 0, sub('gene_name "(\\S+)"', '\\1', matches), NA)
  return(gene_name)
}
gene_names_M <- sapply(gtf_table_M$V9, extract_gene_name)
gtf_table_M$gene_names_M <- gene_names_M
gene_names_R <- sapply(gtf_table_R$V9, extract_gene_name)
gtf_table_R$gene_names_R <- gene_names_R

# Remove rows where 'gene_name' is NA
gtf_table_M <- gtf_table_M %>% filter(!is.na(gene_names_M))
gtf_table_R <- gtf_table_R %>% filter(!is.na(gene_names_R))
