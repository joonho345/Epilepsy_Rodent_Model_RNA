library(dplyr)

#### import RAW Count matrix adjusted_1 VERSION ####
adjusted_df <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/02.Quantification/adjusted_merged_matrix_1.txt",
                          sep = "\t", header = TRUE, row.names = 1)

#### import TPM matrix adjusted_1 VERSION ####
adjusted_df <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/03.Normalization/adjusted_merged_matrix_1_TPM.txt",
                          sep = "\t", header = TRUE, row.names = 1)

#### import coldata ####
coldata_df <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/02.Quantification/filtered_coldata.txt",
                         sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE, 
                         fill = TRUE, quote = "", comment.char = "")

##################################
#### import FILTERED DESeq log2fc data ####
DESeq_l2fc_df_FILTERED <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/07.DESeq/DESeq_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)


#### import FILTERED DESeq log2fc data ####
DESeq_l2fc_Astro_FILTERED <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/06.DEG_02DESeq/DESeq_Astro_ALL.txt",
                              sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_ExN_FILTERED <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/06.DEG_02DESeq/DESeq_ExN_ALL.txt",
                                sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_Micro_FILTERED <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/06.DEG_02DESeq/DESeq_Micro_ALL.txt",
                                sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_Oligo_FILTERED <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/06.DEG_02DESeq/DESeq_Oligo_ALL.txt",
                                sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_InN_FILTERED <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/06.DEG_02DESeq/DESeq_InN_ALL.txt",
                                sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_ExN_DG_FILTERED <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/06.DEG_02DESeq/DESeq_ExN_DG_ALL.txt",
                                sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_ExN_CA_FILTERED <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/06.DEG_02DESeq/DESeq_ExN_CA_ALL.txt",
                                sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

DESeq_l2fc_ExN_SUB_FILTERED <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/06.DEG_02DESeq/DESeq_ExN_SUB_ALL.txt",
                                sep = "\t", header = TRUE, row.names = 1, fill = TRUE)

