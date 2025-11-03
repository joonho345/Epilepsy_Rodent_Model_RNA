library(dplyr)


#### import coldata_df from 00.Matrix_H.R
coldata_df

#### import HiRES matrix ####
## cluster_1 ##
celltype <- 'ExN'
celltype <- 'InN'
celltype <- 'Astro'
celltype <- 'Micro'
# GSE185862 all cluster_1
Dataset_M <- 'GSE185862'
type <- 'all'
clustertype <- 'cluster_1'
file_path <- paste0("/home/joonho345/3_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_", Dataset_M,
                      "_04HiRes/", Dataset_M, "_scRNA_matrix_", type, "_", clustertype, "/Merged_", celltype, "_1.txt")

## cluster_2 ##
celltype <- 'ExN_DG'
celltype <- 'ExN_CA'
celltype <- 'ExN_SUB'
# GSE185862 all cluster_2
Dataset_M <- 'GSE185862'
type <- 'all'
clustertype <- 'cluster_2'
file_path <- paste0("/home/joonho345/3_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_", Dataset_M,
                      "_04HiRes/", Dataset_M, "_scRNA_matrix_", type, "_", clustertype, "/Merged_", celltype, "_1.txt")

# print version
version <- paste0(celltype, Dataset_M, type, clustertype)
print(version)
Expr_matrix <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#### choose groups to compare #### 
coldata_df_M_PILO_O_HA <- coldata_df_M_B %>% filter(TYPE_1 == "M_PILO_IP_O_HA" | TYPE_1 == "CTL") 
filtered_samples_M_PILO_O_HA <- rownames(coldata_df_M_PILO_O_HA)
adjusted_df_M_PILO_O_HA <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_PILO_O_HA]

coldata_df_M_PILO_O_IM <- coldata_df_M_B %>% filter(TYPE_1 == "M_PILO_IP_O_O_IM" | TYPE_1 == "CTL") 
filtered_samples_M_PILO_O_IM <- rownames(coldata_df_M_PILO_O_IM)
adjusted_df_M_PILO_O_IM <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_PILO_O_IM]

coldata_df_M_PILO_O_CR <- coldata_df_M_B %>% filter(TYPE_1 == "M_PILO_IP_O_CR" | TYPE_1 == "CTL")
filtered_samples_M_PILO_O_CR <- rownames(coldata_df_M_PILO_O_CR)
adjusted_df_M_PILO_O_CR <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_PILO_O_CR]

coldata_df_M_PILO_Y_HA <- coldata_df_M %>% filter(TYPE_1 == "M_PILO_IP_Y_HA")
filtered_samples_M_PILO_Y_HA <- rownames(coldata_df_M_PILO_Y_HA)
adjusted_df_M_PILO_Y_HA <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_PILO_Y_HA]

coldata_df_M_PILO_Y_AC <- coldata_df_M %>% filter(TYPE_1 == "M_PILO_IP_Y_AC")
filtered_samples_M_PILO_Y_AC <- rownames(coldata_df_M_PILO_Y_AC)
adjusted_df_M_PILO_Y_AC <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_PILO_Y_AC]

coldata_df_M_PILO_Y_IM <- coldata_df_M_L %>% filter(TYPE_1 == "M_PILO_IP_Y_IM" | TYPE == "HIP_W_CTL")
filtered_samples_M_PILO_Y_IM <- rownames(coldata_df_M_PILO_Y_IM)
adjusted_df_M_PILO_Y_IM <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_PILO_Y_IM]

coldata_df_M_KAI_IP_O_HA <- coldata_df_M_J %>% filter(TYPE_1 == "M_KAI_IP_O_HA")
filtered_samples_M_KAI_IP_O_HA <- rownames(coldata_df_M_KAI_IP_O_HA)
adjusted_df_M_KAI_IP_O_HA <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_IP_O_HA]

coldata_df_M_KAI_ST_BO_O_HA <- coldata_df_M_E %>% filter(TYPE_1 == "M_KAI_ST_BO_O_HA")
filtered_samples_M_KAI_ST_BO_O_HA <- rownames(coldata_df_M_KAI_ST_BO_O_HA)
adjusted_df_M_KAI_ST_BO_O_HA <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_ST_BO_O_HA]

coldata_df_M_KAI_ST_BO_O_AC <- coldata_df_M_E %>% filter(TYPE_1 == "M_KAI_ST_BO_O_AC")
filtered_samples_M_KAI_ST_BO_O_AC <- rownames(coldata_df_M_KAI_ST_BO_O_AC)
adjusted_df_M_KAI_ST_BO_O_AC <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_ST_BO_O_AC]

coldata_df_M_KAI_ST_BO_O_IM <- coldata_df_M_E %>% filter(TYPE_1 == "M_KAI_ST_BO_O_IM")
filtered_samples_M_KAI_ST_BO_O_IM <- rownames(coldata_df_M_KAI_ST_BO_O_IM)
adjusted_df_M_KAI_ST_BO_O_IM <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_ST_BO_O_IM]

coldata_df_M_KAI_ST_BO_O_CR <- coldata_df_M_E %>% filter(TYPE_1 == "M_KAI_ST_BO_O_CR")
filtered_samples_M_KAI_ST_BO_O_CR <- rownames(coldata_df_M_KAI_ST_BO_O_CR)
adjusted_df_M_KAI_ST_BO_O_CR <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_ST_BO_O_CR]

coldata_df_M_KAI_ST_IPSI_Y_AC <- coldata_df_M_O %>% filter(TYPE_1 == "M_KAI_ST_IPSI_Y_AC")
filtered_samples_M_KAI_ST_IPSI_Y_AC <- rownames(coldata_df_M_KAI_ST_IPSI_Y_AC)
adjusted_df_M_KAI_ST_IPSI_Y_AC <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_ST_IPSI_Y_AC]

coldata_df_M_KAI_ST_IPSI_Y_IM <- coldata_df_M_O %>% filter(TYPE_1 == "M_KAI_ST_IPSI_Y_IM")
filtered_samples_M_KAI_ST_IPSI_Y_IM <- rownames(coldata_df_M_KAI_ST_IPSI_Y_IM)
adjusted_df_M_KAI_ST_IPSI_Y_IM <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_ST_IPSI_Y_IM]

coldata_df_M_KAI_ST_IPSI_Y_CR <- coldata_df_M_O %>% filter(TYPE_1 == "M_KAI_ST_IPSI_Y_CR")
filtered_samples_M_KAI_ST_IPSI_Y_CR <- rownames(coldata_df_M_KAI_ST_IPSI_Y_CR)
adjusted_df_M_KAI_ST_IPSI_Y_CR <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_ST_IPSI_Y_CR]

coldata_df_M_KAI_ST_CON_Y_AC <- coldata_df_M_O %>% filter(TYPE_1 == "M_KAI_ST_CON_Y_AC")
filtered_samples_M_KAI_ST_CON_Y_AC <- rownames(coldata_df_M_KAI_ST_CON_Y_AC)
adjusted_df_M_KAI_ST_CON_Y_AC <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_ST_CON_Y_AC]

coldata_df_M_KAI_ST_CON_Y_IM <- coldata_df_M_O %>% filter(TYPE_1 == "M_KAI_ST_CON_Y_IM")
filtered_samples_M_KAI_ST_CON_Y_IM <- rownames(coldata_df_M_KAI_ST_CON_Y_IM)
adjusted_df_M_KAI_ST_CON_Y_IM <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_ST_CON_Y_IM]

coldata_df_M_KAI_ST_CON_Y_CR <- coldata_df_M_O %>% filter(TYPE_1 == "M_KAI_ST_CON_Y_CR")
filtered_samples_M_KAI_ST_CON_Y_CR <- rownames(coldata_df_M_KAI_ST_CON_Y_CR)
adjusted_df_M_KAI_ST_CON_Y_CR <- Expr_matrix[, colnames(Expr_matrix) %in% filtered_samples_M_KAI_ST_CON_Y_CR]




