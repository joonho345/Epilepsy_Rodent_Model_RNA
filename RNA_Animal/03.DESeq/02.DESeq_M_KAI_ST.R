library(DESeq2)
library(ggplot2)
library(pasilla)


### Choose one matrix to analyze ###
coldata_df_M_KAI_ST_BO_O_HA <- coldata_df_M_E %>% filter(TYPE_1 == "M_KAI_ST_BO_O_HA")
filtered_samples_M_KAI_ST_BO_O_HA <- rownames(coldata_df_M_KAI_ST_BO_O_HA)
adjusted_df_M_KAI_ST_BO_O_HA <- adjusted_df_M_E[, colnames(adjusted_df_M_E) %in% filtered_samples_M_KAI_ST_BO_O_HA]
dds_data <- adjusted_df_M_KAI_ST_BO_O_HA
dds_coldata <- coldata_df_M_KAI_ST_BO_O_HA
Target_comparison <- 'M_KAI_ST_BO_O_HA'

coldata_df_M_KAI_ST_BO_O_AC <- coldata_df_M_E %>% filter(TYPE_1 == "M_KAI_ST_BO_O_AC")
filtered_samples_M_KAI_ST_BO_O_AC <- rownames(coldata_df_M_KAI_ST_BO_O_AC)
adjusted_df_M_KAI_ST_BO_O_AC <- adjusted_df_M_E[, colnames(adjusted_df_M_E) %in% filtered_samples_M_KAI_ST_BO_O_AC]
dds_data <- adjusted_df_M_KAI_ST_BO_O_AC
dds_coldata <- coldata_df_M_KAI_ST_BO_O_AC
Target_comparison <- 'M_KAI_ST_BO_O_AC'

coldata_df_M_KAI_ST_BO_O_IM <- coldata_df_M_E %>% filter(TYPE_1 == "M_KAI_ST_BO_O_IM")
filtered_samples_M_KAI_ST_BO_O_IM <- rownames(coldata_df_M_KAI_ST_BO_O_IM)
adjusted_df_M_KAI_ST_BO_O_IM <- adjusted_df_M_E[, colnames(adjusted_df_M_E) %in% filtered_samples_M_KAI_ST_BO_O_IM]
dds_data <- adjusted_df_M_KAI_ST_BO_O_IM
dds_coldata <- coldata_df_M_KAI_ST_BO_O_IM
Target_comparison <- 'M_KAI_ST_BO_O_IM'

coldata_df_M_KAI_ST_BO_O_CR <- coldata_df_M_E %>% filter(TYPE_1 == "M_KAI_ST_BO_O_CR")
filtered_samples_M_KAI_ST_BO_O_CR <- rownames(coldata_df_M_KAI_ST_BO_O_CR)
adjusted_df_M_KAI_ST_BO_O_CR <- adjusted_df_M_E[, colnames(adjusted_df_M_E) %in% filtered_samples_M_KAI_ST_BO_O_CR]
dds_data <- adjusted_df_M_KAI_ST_BO_O_CR
dds_coldata <- coldata_df_M_KAI_ST_BO_O_CR
Target_comparison <- 'M_KAI_ST_BO_O_CR'

contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_BOTH","CTL_SAL_HIPPO_BOTH")

################################################
coldata_df_M_KAI_ST_IPSI_Y_AC <- coldata_df_M_O %>% filter(TYPE_1 == "M_KAI_ST_IPSI_Y_AC")
filtered_samples_M_KAI_ST_IPSI_Y_AC <- rownames(coldata_df_M_KAI_ST_IPSI_Y_AC)
adjusted_df_M_KAI_ST_IPSI_Y_AC <- adjusted_df_M_O[, colnames(adjusted_df_M_O) %in% filtered_samples_M_KAI_ST_IPSI_Y_AC]
dds_data <- adjusted_df_M_KAI_ST_IPSI_Y_AC
dds_coldata <- coldata_df_M_KAI_ST_IPSI_Y_AC
Target_comparison <- 'M_KAI_ST_IPSI_Y_AC'

coldata_df_M_KAI_ST_IPSI_Y_IM <- coldata_df_M_O %>% filter(TYPE_1 == "M_KAI_ST_IPSI_Y_IM")
filtered_samples_M_KAI_ST_IPSI_Y_IM <- rownames(coldata_df_M_KAI_ST_IPSI_Y_IM)
adjusted_df_M_KAI_ST_IPSI_Y_IM <- adjusted_df_M_O[, colnames(adjusted_df_M_O) %in% filtered_samples_M_KAI_ST_IPSI_Y_IM]
dds_data <- adjusted_df_M_KAI_ST_IPSI_Y_IM
dds_coldata <- coldata_df_M_KAI_ST_IPSI_Y_IM
Target_comparison <- 'M_KAI_ST_IPSI_Y_IM'

coldata_df_M_KAI_ST_IPSI_Y_CR <- coldata_df_M_O %>% filter(TYPE_1 == "M_KAI_ST_IPSI_Y_CR")
filtered_samples_M_KAI_ST_IPSI_Y_CR <- rownames(coldata_df_M_KAI_ST_IPSI_Y_CR)
adjusted_df_M_KAI_ST_IPSI_Y_CR <- adjusted_df_M_O[, colnames(adjusted_df_M_O) %in% filtered_samples_M_KAI_ST_IPSI_Y_CR]
dds_data <- adjusted_df_M_KAI_ST_IPSI_Y_CR
dds_coldata <- coldata_df_M_KAI_ST_IPSI_Y_CR
Target_comparison <- 'M_KAI_ST_IPSI_Y_CR'

contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_IPSI","CTL_SAL_HIPPO_IPSI")

################################################
coldata_df_M_KAI_ST_CON_Y_AC <- coldata_df_M_O %>% filter(TYPE_1 == "M_KAI_ST_CON_Y_AC")
filtered_samples_M_KAI_ST_CON_Y_AC <- rownames(coldata_df_M_KAI_ST_CON_Y_AC)
adjusted_df_M_KAI_ST_CON_Y_AC <- adjusted_df_M_O[, colnames(adjusted_df_M_O) %in% filtered_samples_M_KAI_ST_CON_Y_AC]
dds_data <- adjusted_df_M_KAI_ST_CON_Y_AC
dds_coldata <- coldata_df_M_KAI_ST_CON_Y_AC
Target_comparison <- 'M_KAI_ST_CON_Y_AC'

coldata_df_M_KAI_ST_CON_Y_IM <- coldata_df_M_O %>% filter(TYPE_1 == "M_KAI_ST_CON_Y_IM")
filtered_samples_M_KAI_ST_CON_Y_IM <- rownames(coldata_df_M_KAI_ST_CON_Y_IM)
adjusted_df_M_KAI_ST_CON_Y_IM <- adjusted_df_M_O[, colnames(adjusted_df_M_O) %in% filtered_samples_M_KAI_ST_CON_Y_IM]
dds_data <- adjusted_df_M_KAI_ST_CON_Y_IM
dds_coldata <- coldata_df_M_KAI_ST_CON_Y_IM
Target_comparison <- 'M_KAI_ST_CON_Y_IM'

coldata_df_M_KAI_ST_CON_Y_CR <- coldata_df_M_O %>% filter(TYPE_1 == "M_KAI_ST_CON_Y_CR")
filtered_samples_M_KAI_ST_CON_Y_CR <- rownames(coldata_df_M_KAI_ST_CON_Y_CR)
adjusted_df_M_KAI_ST_CON_Y_CR <- adjusted_df_M_O[, colnames(adjusted_df_M_O) %in% filtered_samples_M_KAI_ST_CON_Y_CR]
dds_data <- adjusted_df_M_KAI_ST_CON_Y_CR
dds_coldata <- coldata_df_M_KAI_ST_CON_Y_CR
Target_comparison <- 'M_KAI_ST_CON_Y_CR'

contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_CONTRA","CTL_SAL_HIPPO_CONTRA")

################################################
# Check if coldata and count matrix have the same order with respect to samples
if (!all(rownames(dds_coldata) == colnames(dds_data))) {
  stop("Rownames of coldata do not match colnames of count data")
}

## dds
dds <- DESeqDataSetFromMatrix(countData = dds_data,
                              colData = dds_coldata,
                              design = ~ Treatment_Lateral)
## Prefiltering
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]
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


