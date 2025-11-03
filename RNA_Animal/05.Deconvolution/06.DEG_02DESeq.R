library(DESeq2)
library(ggplot2)
library(pasilla)


### Choose one matrix to analyze ###
### raw data
#1
dds_data <- adjusted_df_M_PILO_O_HA
dds_coldata <- coldata_df_M_PILO_O_HA
Target_comparison <- 'M_PILO_IP_O_HA'
contrast_vector <- c("Treatment_Lateral","PILO","CTL_SAL_IP")
#2
dds_data <- adjusted_df_M_PILO_O_IM
dds_coldata <- coldata_df_M_PILO_O_IM
Target_comparison <- 'M_PILO_IP_O_IM'
contrast_vector <- c("Treatment_Lateral","PILO","CTL_SAL_IP")
#3
dds_data <- adjusted_df_M_PILO_O_CR
dds_coldata <- coldata_df_M_PILO_O_CR
Target_comparison <- 'M_PILO_IP_O_CR'
contrast_vector <- c("Treatment_Lateral","PILO","CTL_SAL_IP")
#4
dds_data <- adjusted_df_M_PILO_Y_HA
dds_coldata <- coldata_df_M_PILO_Y_HA
Target_comparison <- 'M_PILO_IP_Y_HA'
contrast_vector <- c("Treatment_Lateral","PILO","CTL_SAL_IP")
#5
dds_data <- adjusted_df_M_PILO_Y_AC
dds_coldata <- coldata_df_M_PILO_Y_AC
Target_comparison <- 'M_PILO_IP_Y_AC'
contrast_vector <- c("Treatment_Lateral","PILO","CTL_SAL_IP")
#6
dds_data <- adjusted_df_M_PILO_Y_IM
dds_coldata <- coldata_df_M_PILO_Y_IM
Target_comparison <- 'M_PILO_IP_Y_IM'
contrast_vector <- c("Treatment_Lateral","PILO","CTL_SAL_IP")
#7
dds_data <- adjusted_df_M_KAI_IP_O_HA
dds_coldata <- coldata_df_M_KAI_IP_O_HA
Target_comparison <- 'M_KAI_IP_O_HA'
contrast_vector <- c("Treatment_Lateral","KAI_IP","CTL_SAL_IP")
#8
dds_data <- adjusted_df_M_KAI_ST_BO_O_HA
dds_coldata <- coldata_df_M_KAI_ST_BO_O_HA
Target_comparison <- 'M_KAI_ST_BO_O_HA'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_BOTH","CTL_SAL_HIPPO_BOTH")
#9
dds_data <- adjusted_df_M_KAI_ST_BO_O_AC
dds_coldata <- coldata_df_M_KAI_ST_BO_O_AC
Target_comparison <- 'M_KAI_ST_BO_O_AC'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_BOTH","CTL_SAL_HIPPO_BOTH")
#10
dds_data <- adjusted_df_M_KAI_ST_BO_O_IM
dds_coldata <- coldata_df_M_KAI_ST_BO_O_IM
Target_comparison <- 'M_KAI_ST_BO_O_IM'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_BOTH","CTL_SAL_HIPPO_BOTH")
#11
dds_data <- adjusted_df_M_KAI_ST_BO_O_CR
dds_coldata <- coldata_df_M_KAI_ST_BO_O_CR
Target_comparison <- 'M_KAI_ST_BO_O_CR'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_BOTH","CTL_SAL_HIPPO_BOTH")
#12
dds_data <- adjusted_df_M_KAI_ST_IPSI_Y_AC
dds_coldata <- coldata_df_M_KAI_ST_IPSI_Y_AC
Target_comparison <- 'M_KAI_ST_IPSI_Y_AC'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_IPSI","CTL_SAL_HIPPO_IPSI")
#13
dds_data <- adjusted_df_M_KAI_ST_IPSI_Y_IM
dds_coldata <- coldata_df_M_KAI_ST_IPSI_Y_IM
Target_comparison <- 'M_KAI_ST_IPSI_Y_IM'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_IPSI","CTL_SAL_HIPPO_IPSI")
#14
dds_data <- adjusted_df_M_KAI_ST_IPSI_Y_CR
dds_coldata <- coldata_df_M_KAI_ST_IPSI_Y_CR
Target_comparison <- 'M_KAI_ST_IPSI_Y_CR'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_IPSI","CTL_SAL_HIPPO_IPSI")
#15
dds_data <- adjusted_df_M_KAI_ST_CON_Y_AC
dds_coldata <- coldata_df_M_KAI_ST_CON_Y_AC
Target_comparison <- 'M_KAI_ST_CON_Y_AC'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_CONTRA","CTL_SAL_HIPPO_CONTRA")
#16
dds_data <- adjusted_df_M_KAI_ST_CON_Y_IM
dds_coldata <- coldata_df_M_KAI_ST_CON_Y_IM
Target_comparison <- 'M_KAI_ST_CON_Y_IM'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_CONTRA","CTL_SAL_HIPPO_CONTRA")
#17
dds_data <- adjusted_df_M_KAI_ST_CON_Y_CR
dds_coldata <- coldata_df_M_KAI_ST_CON_Y_CR
Target_comparison <- 'M_KAI_ST_CON_Y_CR'
contrast_vector <- c("Treatment_Lateral","KAI_HIPPO_CONTRA","CTL_SAL_HIPPO_CONTRA")


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
file_name <- paste0("/home/joonho345/3_RNA/RNA_Animal/06.Deconvolution/06.DEG_02DESeq/",
                    celltype, "_",  Dataset_M, "_", Target_comparison, ".txt")
write.table(res_df, file_name, sep = "\t", quote = FALSE)


