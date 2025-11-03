library(dplyr)

#### import RAW Count matrix adjusted VERSION ####
adjusted_df_M <- read.table("/home/joonho345/3_RNA/RNA_Animal/02.Quantification/adjusted_merged_matrix_M.txt",
                          sep = "\t", header = TRUE, row.names = 1)
adjusted_df_R <- read.table("/home/joonho345/3_RNA/RNA_Animal/02.Quantification/adjusted_merged_matrix_R.txt",
                            sep = "\t", header = TRUE, row.names = 1)

#### import TPM matrix adjusted VERSION ####
adjusted_df_M <- read.table("/home/joonho345/3_RNA/RNA_Animal/04.Normalization/adjusted_merged_matrix_TPM_M.txt",
                            sep = "\t", header = TRUE, row.names = 1)
adjusted_df_R <- read.table("/home/joonho345/3_RNA/RNA_Animal/04.Normalization/adjusted_merged_matrix_TPM_R.txt",
                            sep = "\t", header = TRUE, row.names = 1)

#### import coldata ####
coldata_df_M <- read.table("/home/joonho345/3_RNA/RNA_Animal/02.Quantification/filtered_coldata_M.txt",
                         sep = "\t", header = TRUE, row.names = 1, fill = TRUE)
coldata_df_R <- read.table("/home/joonho345/3_RNA/RNA_Animal/02.Quantification/filtered_coldata_R.txt",
                           sep = "\t", header = TRUE, row.names = 1, fill = TRUE)


##################################
#### import DESeq log2fc data _1 ####
#DESeq_l2fc_df_M <- read.table("/home/joonho345/3_RNA/RNA_Animal/03.DESeq/DESeq_M_ALL.txt",
#                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)
#DESeq_l2fc_df_R <- read.table("/home/joonho345/3_RNA/RNA_Animal/03.DESeq/DESeq_R_ALL.txt",
#                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

#### IMPORT DESeq coldata _1 ####
#DESeq_coldata_df_M <- read.table("/home/joonho345/3_RNA/RNA_Animal/03.DESeq/DESeq_coldata_1_M.txt",
#                           sep = "\t", header = TRUE, row.names = 1, fill = TRUE)
#DESeq_coldata_df_R <- read.table("/home/joonho345/3_RNA/RNA_Animal/03.DESeq/DESeq_coldata_1_R.txt",
#                           sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

#### import DESeq log2fc data _2 ####
DESeq_l2fc_df_M <- read.table("/home/joonho345/3_RNA/RNA_Animal/03.DESeq_2/DESeq_M_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)
DESeq_l2fc_df_R <- read.table("/home/joonho345/3_RNA/RNA_Animal/03.DESeq_2/DESeq_R_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

#### IMPORT DESeq coldata _2 ####
DESeq_coldata_df_M <- read.table("/home/joonho345/3_RNA/RNA_Animal/03.DESeq_2/DESeq_coldata_2_M.txt",
                                 sep = "\t", header = TRUE, row.names = 1, fill = TRUE)
DESeq_coldata_df_R <- read.table("/home/joonho345/3_RNA/RNA_Animal/03.DESeq_2/DESeq_coldata_2_R.txt",
                                 sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

##################################
#merged_df_M <- read.table("/home/joonho345/3_RNA/RNA_Animal/02.Quantification/merged_matrix_M.txt",
#                         sep = "\t", header = TRUE, row.names = 1)
#merged_matrix_M <- as.matrix(merged_df_M)
#merged_df_R <- read.table("/home/joonho345/3_RNA/RNA_Animal/02.Quantification/merged_matrix_R.txt",
#                         sep = "\t", header = TRUE, row.names = 1)
#merged_matrix_R <- as.matrix(merged_df_R)

#adjusted_merged_df_M <- read.table("/home/joonho345/3_RNA/RNA_Animal/02.Quantification/adjusted_merged_matrix_M.txt",
#                                 sep = "\t", header = TRUE, row.names = 1)
#adjusted_merged_matrix_M <- as.matrix(adjusted_merged_df_M)
#adjusted_merged_df_R <- read.table("/home/joonho345/3_RNA/RNA_Animal/02.Quantification/adjusted_merged_matrix_R.txt",
#                                 sep = "\t", header = TRUE, row.names = 1)
#adjusted_merged_matrix_R <- as.matrix(adjusted_merged_df_R)

#adjusted_merged_df_M_1 <- read.table("/home/joonho345/3_RNA/RNA_Animal/02.Quantification/adjusted_merged_matrix_M_1.txt",
#                                   sep = "\t", header = TRUE, row.names = 1)
#adjusted_merged_matrix_M_1 <- as.matrix(adjusted_merged_df_M_1)
#adjusted_merged_df_R_1 <- read.table("/home/joonho345/3_RNA/RNA_Animal/02.Quantification/adjusted_merged_matrix_R_1.txt",
#                                   sep = "\t", header = TRUE, row.names = 1)
#adjusted_merged_matrix_R_1 <- as.matrix(adjusted_merged_df_R_1)

##################################
#### import DESeq log2fc data _2 ####
DESeq_l2fc_Astro_M <- read.table("/home/joonho345/3_RNA/RNA_Animal/05.Deconvolution/06.DEG_02DESeq/DESeq_Astro_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_ExN_M <- read.table("/home/joonho345/3_RNA/RNA_Animal/05.Deconvolution/06.DEG_02DESeq/DESeq_ExN_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_Micro_M <- read.table("/home/joonho345/3_RNA/RNA_Animal/05.Deconvolution/06.DEG_02DESeq/DESeq_Micro_ALL.txt",                           
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_InN_M <- read.table("/home/joonho345/3_RNA/RNA_Animal/05.Deconvolution/06.DEG_02DESeq/DESeq_InN_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_ExN_DG_M <- read.table("/home/joonho345/3_RNA/RNA_Animal/05.Deconvolution/06.DEG_02DESeq/DESeq_ExN_DG_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_ExN_CA_M <- read.table("/home/joonho345/3_RNA/RNA_Animal/05.Deconvolution/06.DEG_02DESeq/DESeq_ExN_CA_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

