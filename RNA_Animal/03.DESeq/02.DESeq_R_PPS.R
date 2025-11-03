library(DESeq2)
library(ggplot2)
library(pasilla)


### Choose one matrix to analyze ###
coldata_df_R_PPS_O_AC <- coldata_df_R_G %>% filter(TYPE_1 == "R_PPS_Y_AC" | TYPE_1 == "PPS_CTL")
filtered_samples_R_PPS_O_AC <- rownames(coldata_df_R_PPS_O_AC)
adjusted_df_R_PPS_O_AC <- adjusted_df_R_G[, colnames(adjusted_df_R_G) %in% filtered_samples_R_PPS_O_AC]
dds_data <- adjusted_df_R_PPS_O_AC
dds_coldata <- coldata_df_R_PPS_O_AC
Target_comparison <- 'R_PPS_Y_AC'

coldata_df_R_PPS_O_IM <- coldata_df_R_G %>% filter(TYPE_1 == "R_PPS_Y_IM" | TYPE_1 == "PPS_CTL")
filtered_samples_R_PPS_O_IM <- rownames(coldata_df_R_PPS_O_IM)
adjusted_df_R_PPS_O_IM <- adjusted_df_R_G[, colnames(adjusted_df_R_G) %in% filtered_samples_R_PPS_O_IM]
dds_data <- adjusted_df_R_PPS_O_IM
dds_coldata <- coldata_df_R_PPS_O_IM
Target_comparison <- 'R_PPS_Y_IM'

coldata_df_R_PPS_O_CR <- coldata_df_R_G %>% filter(TYPE_1 == "R_PPS_Y_CR" | TYPE_1 == "PPS_CTL")
filtered_samples_R_PPS_O_CR <- rownames(coldata_df_R_PPS_O_CR)
adjusted_df_R_PPS_O_CR <- adjusted_df_R_G[, colnames(adjusted_df_R_G) %in% filtered_samples_R_PPS_O_CR]
dds_data <- adjusted_df_R_PPS_O_CR
dds_coldata <- coldata_df_R_PPS_O_CR
Target_comparison <- 'R_PPS_Y_CR'

coldata_df_R_PPS_O_DOFS <- coldata_df_R_G %>% filter(TYPE_1 == "R_PPS_Y_DOFS" | TYPE_1 == "PPS_CTL")
filtered_samples_R_PPS_O_DOFS <- rownames(coldata_df_R_PPS_O_DOFS)
adjusted_df_R_PPS_O_DOFS <- adjusted_df_R_G[, colnames(adjusted_df_R_G) %in% filtered_samples_R_PPS_O_DOFS]
dds_data <- adjusted_df_R_PPS_O_DOFS
dds_coldata <- coldata_df_R_PPS_O_DOFS
Target_comparison <- 'R_PPS_Y_DOFS'

################################################
contrast_vector <- c("Treatment_Lateral","PPS","CTL_PPS")

# Check if coldata and count matrix have the same order with respect to samples
if (!all(rownames(dds_coldata) == colnames(dds_data))) {
  stop("Rownames of coldata do not match colnames of count data")
}

## dds
dds <- DESeqDataSetFromMatrix(countData = dds_data,
                              colData = dds_coldata,
                              design = ~ Treatment_Lateral)
## Prefiltering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
## Differential expression analysis
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, contrast=contrast_vector)

## filter out the rows where the adjusted p-value (padj) is not available (NA)
res <- res[!is.na(res$padj),]
res_df <- as.data.frame(res@listData)
rownames(res_df) <- res@rownames
res_df$'-log10(padj)' <- -log10(res_df$padj)

#### Significance: log2FoldChange >= 2 & padj < 0.05 ####
res_df$sig <- ifelse(res_df$log2FoldChange > 1 & res$padj < 0.05, "Up Regulated",
                     ifelse(res_df$log2FoldChange < -1 & res$padj < 0.05, "Down Regulated", "No change"))
# filter only significant genes
res_sig <- res[abs(res$log2FoldChange) >= 1 & res$padj < 0.05,]
res_sig_df <- as.data.frame(res_sig@listData)

#### Export files
file_name <- paste("/home/joonho345/3_RNA/RNA_Animal/03.DESeq_2/02.DESeq_files/DESeq_",
                   Target_comparison, ".txt", sep="")
write.table(res_df, file_name, sep = "\t", quote = FALSE)


