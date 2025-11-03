library(dplyr)
library(tidyr)

gtf_path <- "/home/joonho345/resources/Reference/Homo_sapiens.GRCh38.112.gtf"
gtf_table <- read.table(gtf_path, header = FALSE, skip = 5, sep = "\t", quote = "", comment.char = "")
gtf_table <- gtf_table %>% dplyr::select(-c(V1, V4, V5, V6, V7, V8))
gtf_table <- gtf_table[gtf_table$V3 == "gene", ]


### Create gene_id column
extract_gene_id <- function(v9_string) {
  # Use a regular expression to find and extract the gene_id value
  matches <- regmatches(v9_string, regexpr('gene_id "\\S+"', v9_string))
  gene_id <- ifelse(length(matches) > 0, sub('gene_id "(\\S+)"', '\\1', matches), NA)
  return(gene_id)
}
gene_id <- sapply(gtf_table$V9, extract_gene_id)
gtf_table$gene_id <- gene_id


### Create gene_name column
extract_gene_name <- function(v9_string) {
  # Use a regular expression to find and extract the gene_id value
  matches <- regmatches(v9_string, regexpr('gene_name "\\S+"', v9_string))
  gene_name <- ifelse(length(matches) > 0, sub('gene_name "(\\S+)"', '\\1', matches), NA)
  return(gene_name)
}
gene_name <- sapply(gtf_table$V9, extract_gene_name)
gtf_table$gene_name <- gene_name


### capitalize gene name
gtf_table$gene_name <- toupper(gtf_table$gene_name)

# Remove rows where 'gene_name' is NA
gtf_table <- gtf_table %>% filter(!is.na(gene_name))
