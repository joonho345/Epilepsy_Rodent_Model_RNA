library(dplyr)

#### import RAW Count matrix adjusted VERSION ####
adjusted_df_M <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/02.Quantification/adjusted_merged_matrix_M.txt",
                          sep = "\t", header = TRUE, row.names = 1)
adjusted_df_R <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/02.Quantification/adjusted_merged_matrix_R.txt",
                            sep = "\t", header = TRUE, row.names = 1)

#### import TPM matrix adjusted VERSION ####
adjusted_df_M <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/03.Normalization/adjusted_merged_matrix_TPM_M.txt",
                            sep = "\t", header = TRUE, row.names = 1)
adjusted_df_R <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/03.Normalization/adjusted_merged_matrix_TPM_R.txt",
                            sep = "\t", header = TRUE, row.names = 1)

#### import coldata ####
coldata_df_M <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/02.Quantification/filtered_coldata_M.txt",
                         sep = "\t", header = TRUE, row.names = 1, fill = TRUE)
coldata_df_R <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/02.Quantification/filtered_coldata_R.txt",
                           sep = "\t", header = TRUE, row.names = 1, fill = TRUE)


##################################
#### import DESeq log2fc data _3 ####
DESeq_l2fc_df_M <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/05.DESeq/DESeq_M_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)
DESeq_l2fc_df_R <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/05.DESeq/DESeq_R_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

#### IMPORT DESeq coldata _3 ####
DESeq_coldata_df_M <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/05.DESeq/DESeq_coldata_3_M.txt",
                                 sep = "\t", header = TRUE, row.names = 1, fill = TRUE)
DESeq_coldata_df_R <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/05.DESeq/DESeq_coldata_3_R.txt",
                                 sep = "\t", header = TRUE, row.names = 1, fill = TRUE)


##################################
#### import DESeq log2fc data _2 ####
DESeq_l2fc_Astro_M <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/06.DEG_02DESeq/DESeq_Astro-M_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_ExN_M <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/06.DEG_02DESeq/DESeq_ExN-M_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_Micro_M <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/06.DEG_02DESeq/DESeq_Micro-M_ALL.txt",                           
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_InN_M <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/06.DEG_02DESeq/DESeq_InN-M_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_Oligo_M <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/06.DEG_02DESeq/DESeq_Oligo-M_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_OPC_M <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/06.DEG_02DESeq/DESeq_OPC-M_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_Endo_M <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/06.DEG_02DESeq/DESeq_Endo-M_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)
