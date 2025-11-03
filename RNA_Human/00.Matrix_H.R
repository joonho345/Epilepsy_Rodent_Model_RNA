library(dplyr)

#### import RAW Count matrix adjusted_1 VERSION ####
adjusted_df <- read.table("/home/joonho345/3_RNA/RNA_Human/02.Quantification/adjusted_merged_matrix_1.txt",
                          sep = "\t", header = TRUE, row.names = 1)

#### import TPM matrix adjusted_1 VERSION ####
adjusted_df <- read.table("/home/joonho345/3_RNA/RNA_Human/04.Normalization/adjusted_merged_matrix_1_TPM.txt",
                          sep = "\t", header = TRUE, row.names = 1)

#### import coldata ####
coldata_df <- read.table("/home/joonho345/3_RNA/RNA_Human/02.Quantification/filtered_coldata.txt",
                         sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

##################################
#### import DESeq log2fc data ####
#DESeq_l2fc_df <- read.table("/home/joonho345/3_RNA/RNA_Human/03.DESeq/DESeq_ALL.txt",
#                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

#### import FILTERED DESeq log2fc data ####
DESeq_l2fc_df_FILTERED <- read.table("/home/joonho345/3_RNA/RNA_Human/06.DESeq/DESeq_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)


##################################
#### import DESeq log2fc data ####
#DESeq_l2fc_Astro <- read.table("/home/joonho345/3_RNA/RNA_Human/06.Deconvolution/06.DEG_02DESeq/DESeq_Astro_ALL.txt",
#                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

#DESeq_l2fc_ExN <- read.table("/home/joonho345/3_RNA/RNA_Human/06.Deconvolution/06.DEG_02DESeq/DESeq_ExN_ALL.txt",
#                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

#DESeq_l2fc_Micro <- read.table("/home/joonho345/3_RNA/RNA_Human/06.Deconvolution/06.DEG_02DESeq/DESeq_Micro_ALL.txt",                           
#                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

#DESeq_l2fc_Oligo <- read.table("/home/joonho345/3_RNA/RNA_Human/06.Deconvolution/06.DEG_02DESeq/DESeq_Oligo_ALL.txt",                        
#                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

#DESeq_l2fc_InN <- read.table("/home/joonho345/3_RNA/RNA_Human/06.Deconvolution/06.DEG_02DESeq/DESeq_InN_ALL.txt",
#                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

#DESeq_l2fc_ExN_DG <- read.table("/home/joonho345/3_RNA/RNA_Human/06.Deconvolution/06.DEG_02DESeq/DESeq_ExN_DG_ALL.txt",
#                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

#DESeq_l2fc_ExN_CA <- read.table("/home/joonho345/3_RNA/RNA_Human/06.Deconvolution/06.DEG_02DESeq/DESeq_ExN_CA_ALL.txt",
#                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

#DESeq_l2fc_ExN_SUB <- read.table("/home/joonho345/3_RNA/RNA_Human/06.Deconvolution/06.DEG_02DESeq/DESeq_ExN_SUB_ALL.txt",
#                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

#### import FILTERED DESeq log2fc data ####
DESeq_l2fc_Astro_FILTERED <- read.table("/home/joonho345/3_RNA/RNA_Human/08.Deconvolution/06.DEG_02DESeq/DESeq_Astro_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_ExN_FILTERED <- read.table("/home/joonho345/3_RNA/RNA_Human/08.Deconvolution/06.DEG_02DESeq/DESeq_ExN_ALL.txt",
                                sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_Micro_FILTERED <- read.table("/home/joonho345/3_RNA/RNA_Human/08.Deconvolution/06.DEG_02DESeq/DESeq_Micro_ALL.txt",
                                sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_Oligo_FILTERED <- read.table("/home/joonho345/3_RNA/RNA_Human/08.Deconvolution/06.DEG_02DESeq/DESeq_Oligo_ALL.txt",
                                sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_InN_FILTERED <- read.table("/home/joonho345/3_RNA/RNA_Human/08.Deconvolution/06.DEG_02DESeq/DESeq_InN_ALL.txt",
                                sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_ExN_DG_FILTERED <- read.table("/home/joonho345/3_RNA/RNA_Human/08.Deconvolution/06.DEG_02DESeq/DESeq_ExN_DG_ALL.txt",
                                sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_ExN_CA_FILTERED <- read.table("/home/joonho345/3_RNA/RNA_Human/08.Deconvolution/06.DEG_02DESeq/DESeq_ExN_CA_ALL.txt",
                                sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_ExN_SUB_FILTERED <- read.table("/home/joonho345/3_RNA/RNA_Human/08.Deconvolution/06.DEG_02DESeq/DESeq_ExN_SUB_ALL.txt",
                                sep = "\t", header = TRUE, row.names = 1, fill = TRUE)


##################################
#merged_df <- read.table("/home/joonho345/3_RNA/RNA_Human/02.Quantification/merged_matrix.txt",
#                          sep = "\t", header = TRUE, row.names = 1)
#merged_matrix <- as.matrix(merged_df)
#
#adjusted_merged_df <- read.table("/home/joonho345/3_RNA/RNA_Human/02.Quantification/adjusted_merged_matrix.txt",
#                        sep = "\t", header = TRUE, row.names = 1)
#adjusted_merged_matrix <- as.matrix(adjusted_merged_df)
#
#adjusted_merged_df_1 <- read.table("/home/joonho345/3_RNA/RNA_Human/02.Quantification/adjusted_merged_matrix_1.txt",
#                                 sep = "\t", header = TRUE, row.names = 1)
#adjusted_merged_matrix_1 <- as.matrix(adjusted_merged_df_1)
#
#adjusted_merged_df_2 <- read.table("/home/joonho345/3_RNA/RNA_Human/02.Quantification/adjusted_merged_matrix_2.txt",
#                                   sep = "\t", header = TRUE, row.names = 1)
#adjusted_merged_matrix_2 <- as.matrix(adjusted_merged_df_2)

