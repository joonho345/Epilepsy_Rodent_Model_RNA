library(DESeq2)
library(ggplot2)
library(pasilla)


### Choose one matrix to analyze ###
### raw data
dds_data <- adjusted_df_M_KAI_IA_IPSI_A_AC
dds_coldata <- coldata_df_M_KAI_IA_IPSI_A_AC
Target_comparison <- 'M_KAI_IA_IPSI_A_AC'
contrast_vector <- c("Treatment_Lateral","KAI_AMG_IPSI","CTL_SAL_AMG_IPSI")

dds_data <- adjusted_df_M_KAI_IA_IPSI_A_IM
dds_coldata <- coldata_df_M_KAI_IA_IPSI_A_IM
Target_comparison <- 'M_KAI_IA_IPSI_A_IM'
contrast_vector <- c("Treatment_Lateral","KAI_AMG_IPSI","CTL_SAL_AMG_IPSI")

dds_data <- adjusted_df_M_KAI_IA_IPSI_A_CR
dds_coldata <- coldata_df_M_KAI_IA_IPSI_A_CR
Target_comparison <- 'M_KAI_IA_IPSI_A_CR'
contrast_vector <- c("Treatment_Lateral","KAI_AMG_IPSI","CTL_SAL_AMG_IPSI")

dds_data <- adjusted_df_M_KAI_IH_IPSI_A_HA
dds_coldata <- coldata_df_M_KAI_IH_IPSI_A_HA
Target_comparison <- 'M_KAI_IH_IPSI_A_HA'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_IPSI","CTL_SAL_HIPPO_IPSI")

dds_data <- adjusted_df_M_KAI_IH_IPSI_A_AC
dds_coldata <- coldata_df_M_KAI_IH_IPSI_A_AC
Target_comparison <- 'M_KAI_IH_IPSI_A_AC'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_IPSI","CTL_SAL_HIPPO_IPSI")

dds_data <- adjusted_df_M_KAI_IH_IPSI_A_IM
dds_coldata <- coldata_df_M_KAI_IH_IPSI_A_IM
Target_comparison <- 'M_KAI_IH_IPSI_A_IM'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_IPSI","CTL_SAL_HIPPO_IPSI")

dds_data <- adjusted_df_M_KAI_IH_IPSI_A_CR
dds_coldata <- coldata_df_M_KAI_IH_IPSI_A_CR
Target_comparison <- 'M_KAI_IH_IPSI_A_CR'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_IPSI","CTL_SAL_HIPPO_IPSI")

dds_data <- adjusted_df_M_KAI_IH_CON_A_AC
dds_coldata <- coldata_df_M_KAI_IH_CON_A_AC
Target_comparison <- 'M_KAI_IH_CON_A_AC'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_CONTRA","CTL_SAL_HIPPO_CONTRA")

dds_data <- adjusted_df_M_KAI_IH_CON_A_IM
dds_coldata <- coldata_df_M_KAI_IH_CON_A_IM
Target_comparison <- 'M_KAI_IH_CON_A_IM'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_CONTRA","CTL_SAL_HIPPO_CONTRA")

dds_data <- adjusted_df_M_KAI_IH_CON_A_CR
dds_coldata <- coldata_df_M_KAI_IH_CON_A_CR
Target_comparison <- 'M_KAI_IH_CON_A_CR'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_CONTRA","CTL_SAL_HIPPO_CONTRA")

dds_data <- adjusted_df_M_KAI_IP_A_HA
dds_coldata <- coldata_df_M_KAI_IP_A_HA
Target_comparison <- 'M_KAI_IP_A_HA'
contrast_vector <- c("Treatment_Lateral","KAI_IP","CTL_SAL_IP")

dds_data <- adjusted_df_M_PILO_IP_A_HA
dds_coldata <- coldata_df_M_PILO_IP_A_HA
Target_comparison <- 'M_PILO_IP_A_HA'
contrast_vector <- c("Treatment_Lateral","PILO","CTL_SAL_IP")

dds_data <- adjusted_df_M_PILO_IP_A_AC
dds_coldata <- coldata_df_M_PILO_IP_A_AC
Target_comparison <- 'M_PILO_IP_A_AC'
contrast_vector <- c("Treatment_Lateral","PILO","CTL_SAL_IP")

dds_data <- adjusted_df_M_PILO_IP_A_IM
dds_coldata <- coldata_df_M_PILO_IP_A_IM
Target_comparison <- 'M_PILO_IP_A_IM'
contrast_vector <- c("Treatment_Lateral","PILO","CTL_SAL_IP")

dds_data <- adjusted_df_M_PILO_IP_A_CR
dds_coldata <- coldata_df_M_PILO_IP_A_CR
Target_comparison <- 'M_PILO_IP_A_CR'
contrast_vector <- c("Treatment_Lateral","PILO","CTL_SAL_IP")



# Check if coldata and count matrix have the same order with respect to samples
if (!all(rownames(dds_coldata) == colnames(dds_data))) {
  stop("Rownames of coldata do not match colnames of count data")
}

## dds
dds_data <- round(dds_data)
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

#### Export files
file_name <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/06.DEG_02DESeq/",
                    celltype, "_",  Dataset_M, "_", Target_comparison, ".txt")
write.table(res_df, file_name, sep = "\t", quote = FALSE)


