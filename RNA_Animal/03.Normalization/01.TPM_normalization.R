library(dplyr)
library(tidyr)
library(readr)

############# Mouse ############
#### Create gene coverage matrix
gtf_path_M <- "/data/resource/reference/mouse/Mus_musculus.GRCm39.112.gtf"
gtf_table_M <- read.table(gtf_path_M, header = FALSE, skip = 5, sep = "\t", quote = "", comment.char = "")

colnames(gtf_table_M) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
exons <- gtf_table_M %>% filter(feature == "exon")
exons$gene_name <- str_extract(exons$attribute, "gene_name \"[^\"]+\"")
exons$gene_name <- str_replace_all(exons$gene_name, "gene_name |\"", "")
exons <- exons %>% mutate(length = end - start + 1)

gene_cov_M <- exons %>%
  group_by(gene_name) %>%
  arrange(start, end) %>%  # Arrange by start position
  mutate(
    adjusted_start = ifelse(row_number() == 1, start, pmax(start, lag(end) + 1)),
    non_overlap_length = ifelse(adjusted_start <= end, end - adjusted_start + 1, 0)
  ) %>%
  summarize(gene_length = sum(non_overlap_length, na.rm = TRUE)) %>%
  ungroup()

#### Calculate TPM - merged_matrix_M
raw_counts_path_M <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/02.Quantification/merged_matrix_M.txt"
output_path_M <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/03.Normalization/merged_matrix_TPM_M.txt"

counts_M <- read.table(raw_counts_path_M, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
counts_M <- as.matrix(counts_M)

# 1. Convert gene length from base pairs to kilobases
matched_cov_M <- gene_cov_M[match(rownames(counts_M), gene_cov_M$gene_name), ]
gene_lengths_kb_M <- matched_cov_M$gene_length / 1000
# 2. Calculate RPK (Reads Per Kilobase) RPK = raw counts / gene length in kilobases
rpk_M <- counts_M / gene_lengths_kb_M
# 3. Calculate scaling factor (sum of RPKs for each sample)
scaling_factors_M <- colSums(rpk_M)
# 4. Calculate TPM, TPM = RPK / scaling factor * 1e6
tpm_matrix_M <- t(t(rpk_M) / scaling_factors_M) * 1e6
head(tpm_matrix_M)

write.table(tpm_matrix_M, file = output_path_M, sep = "\t", quote = FALSE)

#### Calculate TPM - adjusted_merged_matrix_M
raw_counts_path_M <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/02.Quantification/adjusted_merged_matrix_M.txt"
output_path_M <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/03.Normalization/adjusted_merged_matrix_TPM_M.txt"

counts_M <- read.table(raw_counts_path_M, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
counts_M <- as.matrix(counts_M)

# 1. Convert gene length from base pairs to kilobases
matched_cov_M <- gene_cov_M[match(rownames(counts_M), gene_cov_M$gene_name), ]
gene_lengths_kb_M <- matched_cov_M$gene_length / 1000
# 2. Calculate RPK (Reads Per Kilobase) RPK = raw counts / gene length in kilobases
rpk_M <- counts_M / gene_lengths_kb_M
# 3. Calculate scaling factor (sum of RPKs for each sample)
scaling_factors_M <- colSums(rpk_M)
# 4. Calculate TPM, TPM = RPK / scaling factor * 1e6
tpm_matrix_M <- t(t(rpk_M) / scaling_factors_M) * 1e6
head(tpm_matrix_M)

write.table(tpm_matrix_M, file = output_path_M, sep = "\t", quote = FALSE)


############# Rat ############
#### Create gene coverage matrix
gtf_path_R <- "/data/resource/reference/mouse/Rattus_norvegicus.mRatBN7.2.112.gtf"
gtf_table_R <- read.table(gtf_path_R, header = FALSE, skip = 5, sep = "\t", quote = "", comment.char = "")

colnames(gtf_table_R) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
exons <- gtf_table_R %>% filter(feature == "exon")
exons$gene_name <- str_extract(exons$attribute, "gene_name \"[^\"]+\"")
exons$gene_name <- str_replace_all(exons$gene_name, "gene_name |\"", "")
exons <- exons %>% mutate(length = end - start + 1)

gene_cov_R <- exons %>%
  group_by(gene_name) %>%
  arrange(start, end) %>%  # Arrange by start position
  mutate(
    adjusted_start = ifelse(row_number() == 1, start, pmax(start, lag(end) + 1)),
    non_overlap_length = ifelse(adjusted_start <= end, end - adjusted_start + 1, 0)
  ) %>%
  summarize(gene_length = sum(non_overlap_length, na.rm = TRUE)) %>%
  ungroup()

#### Calculate TPM - merged_matrix_R
raw_counts_path_R <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/02.Quantification/merged_matrix_R.txt"
output_path_R <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/03.Normalization/merged_matrix_TPM_R.txt"

counts_R <- read.table(raw_counts_path_R, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
counts_R <- as.matrix(counts_R)

# 1. Convert gene length from base pairs to kilobases
matched_cov_R <- gene_cov_R[match(rownames(counts_R), gene_cov_R$gene_name), ]
gene_lengths_kb_R <- matched_cov_R$gene_length / 1000
# 2. Calculate RPK (Reads Per Kilobase) RPK = raw counts / gene length in kilobases
rpk_R <- counts_R / gene_lengths_kb_R
# 3. Calculate scaling factor (sum of RPKs for each sample)
scaling_factors_R <- colSums(rpk_R)
# 4. Calculate TPM, TPM = RPK / scaling factor * 1e6
tpm_matrix_R <- t(t(rpk_R) / scaling_factors_R) * 1e6
head(tpm_matrix_R)

write.table(tpm_matrix_R, file = output_path_R, sep = "\t", quote = FALSE)

#### Calculate TPM - adjusted_merged_matrix_R
raw_counts_path_R <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/02.Quantification/adjusted_merged_matrix_R.txt"
output_path_R <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/03.Normalization/adjusted_merged_matrix_TPM_R.txt"

counts_R <- read.table(raw_counts_path_R, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
counts_R <- as.matrix(counts_R)

# 1. Convert gene length from base pairs to kilobases
matched_cov_R <- gene_cov_R[match(rownames(counts_R), gene_cov_R$gene_name), ]
gene_lengths_kb_R <- matched_cov_R$gene_length / 1000
# 2. Calculate RPK (Reads Per Kilobase) RPK = raw counts / gene length in kilobases
rpk_R <- counts_R / gene_lengths_kb_R
# 3. Calculate scaling factor (sum of RPKs for each sample)
scaling_factors_R <- colSums(rpk_R)
# 4. Calculate TPM, TPM = RPK / scaling factor * 1e6
tpm_matrix_R <- t(t(rpk_R) / scaling_factors_R) * 1e6
head(tpm_matrix_R)

write.table(tpm_matrix_R, file = output_path_R, sep = "\t", quote = FALSE)

