library(DESeq2)
library(ggplot2)
library(pasilla)


### Choose one matrix to analyze ###
### raw data
dds_data <- Expr_matrix_1
dds_coldata <- coldata_df_1
Target_comparison <- 'FILTERED_1_MTLEALL_NL'
contrast_vector <- c("Diagnosis","MTLEALL","NL")

## dds
dds_data <- round(dds_data)
dds <- DESeqDataSetFromMatrix(countData = dds_data,
                              colData = dds_coldata,
                              design = ~ Diagnosis)
## Prefiltering
keep <- rowSums(counts(dds)) >= 100
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

#### Export files
file_name <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/06.DEG_02DESeq/",
                   celltype, "_",  Dataset_H, "_", Target_comparison, ".txt")
write.table(res_df, file_name, sep = "\t", quote = FALSE)



