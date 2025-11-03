#######################################
library(DESeq2)
library(ggplot2)
library(pasilla)

#### Wihtout Outliers -> New DESeq ####
### Choose one matrix to analyze ###
### raw data
dds_data <- adjusted_df_1
dds_coldata <- coldata_df_1
Target_comparison <- 'FILTERED_1_MTLEALL_NL'
contrast_vector <- c("Diagnosis","MTLEALL","NL")

dds_data <- adjusted_df_2
dds_coldata <- coldata_df_2
Target_comparison <- 'FILTERED_2_MTLEHS_NL'
contrast_vector <- c("Diagnosis","MTLEHS","NL")

dds_data <- adjusted_df_3
dds_coldata <- coldata_df_3
Target_comparison <- 'FILTERED_3_MTLE_NL'
contrast_vector <- c("Diagnosis","MTLE","NL")


# Check if coldata and count matrix have the same order with respect to samples
if (!all(rownames(dds_coldata) == colnames(dds_data))) {
  stop("Rownames of coldata do not match colnames of count data")
}


## dds
dds <- DESeqDataSetFromMatrix(countData = dds_data,
                              colData = dds_coldata,
                              design = ~ Diagnosis)
## Prefiltering
keep <- rowSums(counts(dds)) >= 200
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
res_df$sig <- ifelse(res_df$log2FoldChange > 2 & res$padj < 0.05, "Up Regulated",
                     ifelse(res_df$log2FoldChange < -2 & res$padj < 0.05, "Down Regulated", "No change"))
# filter only significant genes
res_sig <- res[abs(res$log2FoldChange) >= 2 & res$padj < 0.05,]
res_sig_df <- as.data.frame(res_sig@listData)


#### Export files
file_name <- paste("/home/joonho345/3_RNA/RNA_Human/06.DESeq/01.DESeq_files/DESeq_",
                   Target_comparison, ".txt", sep="")
write.table(res_df, file_name, sep = "\t", quote = FALSE)





