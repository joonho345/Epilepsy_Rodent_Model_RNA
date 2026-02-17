library(DESeq2)
library(ggplot2)
library(pasilla)

# Define all comparisons
comparisons <- list(
  list(name = 'HS', data = 'adjusted_df_HS', coldata = 'coldata_df_HS'),
  list(name = 'withoutHS', data = 'adjusted_df_withoutHS', coldata = 'coldata_df_withoutHS'),
  list(name = 'IP', data = 'adjusted_df_IP', coldata = 'coldata_df_IP'),
  list(name = 'NO', data = 'adjusted_df_NO', coldata = 'coldata_df_NO'),
  list(name = 'IP_Sei', data = 'adjusted_df_IP_Sei', coldata = 'coldata_df_IP_Sei'),
  list(name = 'IP_Other', data = 'adjusted_df_IP_Other', coldata = 'coldata_df_IP_Other')
)

contrast_vector <- c("Diagnosis","MTLEALL","NL")

# Loop through each comparison
for (i in 1:length(comparisons)) {
  Target_comparison <- comparisons[[i]]$name
  dds_data <- get(comparisons[[i]]$data)
  dds_coldata <- get(comparisons[[i]]$coldata)
  
  cat("\nProcessing:", Target_comparison, "\n")

# Check if coldata and count matrix have the same order with respect to samples
if (!all(rownames(dds_coldata) == colnames(dds_data))) {
    stop(paste("Rownames of coldata do not match colnames of count data for", Target_comparison))
}

## dds
dds <- DESeqDataSetFromMatrix(countData = dds_data,
                              colData = dds_coldata,
                              design = ~ Diagnosis)
## Prefiltering
keep <- rowSums(counts(dds)) >= 400
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
  file_name <- paste("/home/joonho345/1_Epilepsy_RNA/RNA_Human/11.Sample_variables/02.deseq/DESeq_",
                   Target_comparison, ".txt", sep="")
write.table(res_df, file_name, sep = "\t", quote = FALSE)
  
  cat("Completed:", Target_comparison, "-", nrow(res_df), "genes\n")
}


