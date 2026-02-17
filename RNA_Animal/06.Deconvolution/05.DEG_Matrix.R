library(dplyr)


#### import coldata_df from 00.Matrix_H.R
coldata_df_M

#### import HiRES matrix ####
celltype <- 'ExN-M'
celltype <- 'InN-M'
celltype <- 'Astro-M'
celltype <- 'Micro-M'
celltype <- 'Oligo-M'
celltype <- 'OPC-M'
celltype <- 'Endo-M'



Dataset_M <- 'GSE185862'
file_path <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_", Dataset_M,
                      "_04HiRes/", Dataset_M, "_10X_scRNA_matrix_DEGs_CellType1/Merged_", celltype, ".txt")

# print version
version <- paste0(celltype, Dataset_M)
print(version)
Expr_matrix <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#### choose groups to compare #### 
coldata_df_M_KAI_IA_IPSI_A_AC <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IA_IPSI_A_AC")
filtered_samples_M_KAI_IA_IPSI_A_AC <- rownames(coldata_df_M_KAI_IA_IPSI_A_AC)
adjusted_df_M_KAI_IA_IPSI_A_AC <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_IA_IPSI_A_AC]

coldata_df_M_KAI_IA_IPSI_A_IM <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IA_IPSI_A_IM")
filtered_samples_M_KAI_IA_IPSI_A_IM <- rownames(coldata_df_M_KAI_IA_IPSI_A_IM)
adjusted_df_M_KAI_IA_IPSI_A_IM <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_IA_IPSI_A_IM]

coldata_df_M_KAI_IA_IPSI_A_CR <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IA_IPSI_A_CR")
filtered_samples_M_KAI_IA_IPSI_A_CR <- rownames(coldata_df_M_KAI_IA_IPSI_A_CR)
adjusted_df_M_KAI_IA_IPSI_A_CR <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_IA_IPSI_A_CR]

coldata_df_M_KAI_IH_IPSI_A_HA <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IH_IPSI_A_HA")
filtered_samples_M_KAI_IH_IPSI_A_HA <- rownames(coldata_df_M_KAI_IH_IPSI_A_HA)
adjusted_df_M_KAI_IH_IPSI_A_HA <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_IH_IPSI_A_HA]

coldata_df_M_KAI_IH_IPSI_A_AC <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IH_IPSI_A_AC")
filtered_samples_M_KAI_IH_IPSI_A_AC <- rownames(coldata_df_M_KAI_IH_IPSI_A_AC)
adjusted_df_M_KAI_IH_IPSI_A_AC <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_IH_IPSI_A_AC]

coldata_df_M_KAI_IH_IPSI_A_IM <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IH_IPSI_A_IM")
filtered_samples_M_KAI_IH_IPSI_A_IM <- rownames(coldata_df_M_KAI_IH_IPSI_A_IM)
adjusted_df_M_KAI_IH_IPSI_A_IM <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_IH_IPSI_A_IM]

coldata_df_M_KAI_IH_IPSI_A_CR <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IH_IPSI_A_CR")
filtered_samples_M_KAI_IH_IPSI_A_CR <- rownames(coldata_df_M_KAI_IH_IPSI_A_CR)
adjusted_df_M_KAI_IH_IPSI_A_CR <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_IH_IPSI_A_CR]

coldata_df_M_KAI_IH_CON_A_AC <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IH_CON_A_AC")
filtered_samples_M_KAI_IH_CON_A_AC <- rownames(coldata_df_M_KAI_IH_CON_A_AC)
adjusted_df_M_KAI_IH_CON_A_AC <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_IH_CON_A_AC]

coldata_df_M_KAI_IH_CON_A_IM <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IH_CON_A_IM")
filtered_samples_M_KAI_IH_CON_A_IM <- rownames(coldata_df_M_KAI_IH_CON_A_IM)
adjusted_df_M_KAI_IH_CON_A_IM <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_IH_CON_A_IM]

coldata_df_M_KAI_IH_CON_A_CR <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IH_CON_A_CR")
filtered_samples_M_KAI_IH_CON_A_CR <- rownames(coldata_df_M_KAI_IH_CON_A_CR)
adjusted_df_M_KAI_IH_CON_A_CR <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_IH_CON_A_CR]

coldata_df_M_KAI_IP_A_HA <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IP_A_HA")
filtered_samples_M_KAI_IP_A_HA <- rownames(coldata_df_M_KAI_IP_A_HA)
adjusted_df_M_KAI_IP_A_HA <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_IP_A_HA]

coldata_df_M_PILO_IP_A_HA <- coldata_df_M %>% filter(TYPE_2 == "M_PILO_IP_A_HA" | TYPE_2 == "CTL")
filtered_samples_M_PILO_IP_A_HA <- rownames(coldata_df_M_PILO_IP_A_HA)
adjusted_df_M_PILO_IP_A_HA <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_PILO_IP_A_HA]

coldata_df_M_PILO_IP_A_AC <- coldata_df_M %>% filter(TYPE_2 == "M_PILO_IP_A_AC")
filtered_samples_M_PILO_IP_A_AC <- rownames(coldata_df_M_PILO_IP_A_AC)
adjusted_df_M_PILO_IP_A_AC <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_PILO_IP_A_AC]

coldata_df_M_PILO_IP_A_IM <- coldata_df_M %>% filter(TYPE_2 == "M_PILO_IP_A_IM" | TYPE_2 == "CTL")
filtered_samples_M_PILO_IP_A_IM <- rownames(coldata_df_M_PILO_IP_A_IM)
adjusted_df_M_PILO_IP_A_IM <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_PILO_IP_A_IM]

coldata_df_M_PILO_IP_A_CR <- coldata_df_M %>% filter(TYPE_2 == "M_PILO_IP_A_CR" | TYPE_2 == "CTL")
filtered_samples_M_PILO_IP_A_CR <- rownames(coldata_df_M_PILO_IP_A_CR)
adjusted_df_M_PILO_IP_A_CR <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_PILO_IP_A_CR]

