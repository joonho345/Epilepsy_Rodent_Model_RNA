library(dplyr)

#### choose groups to compare #### 
####### COLDATA_3 ####### 
coldata_df_M_KAI_IH_IPSI_A_HA <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IH_IPSI_A_HA")
filtered_samples_M_KAI_IH_IPSI_A_HA <- rownames(coldata_df_M_KAI_IH_IPSI_A_HA)
adjusted_df_M_KAI_IH_IPSI_A_HA <- adjusted_df_M[, colnames(adjusted_df_M) %in% filtered_samples_M_KAI_IH_IPSI_A_HA]

coldata_df_M_KAI_IH_IPSI_A_AC <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IH_IPSI_A_AC")
filtered_samples_M_KAI_IH_IPSI_A_AC <- rownames(coldata_df_M_KAI_IH_IPSI_A_AC)
adjusted_df_M_KAI_IH_IPSI_A_AC <- adjusted_df_M[, colnames(adjusted_df_M) %in% filtered_samples_M_KAI_IH_IPSI_A_AC]

coldata_df_M_KAI_IH_IPSI_A_IM <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IH_IPSI_A_IM")
filtered_samples_M_KAI_IH_IPSI_A_IM <- rownames(coldata_df_M_KAI_IH_IPSI_A_IM)
adjusted_df_M_KAI_IH_IPSI_A_IM <- adjusted_df_M[, colnames(adjusted_df_M) %in% filtered_samples_M_KAI_IH_IPSI_A_IM]

coldata_df_M_KAI_IH_IPSI_A_CR <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IH_IPSI_A_CR")
filtered_samples_M_KAI_IH_IPSI_A_CR <- rownames(coldata_df_M_KAI_IH_IPSI_A_CR)
adjusted_df_M_KAI_IH_IPSI_A_CR <- adjusted_df_M[, colnames(adjusted_df_M) %in% filtered_samples_M_KAI_IH_IPSI_A_CR]

coldata_df_M_KAI_IH_CON_A_AC <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IH_CON_A_AC")
filtered_samples_M_KAI_IH_CON_A_AC <- rownames(coldata_df_M_KAI_IH_CON_A_AC)
adjusted_df_M_KAI_IH_CON_A_AC <- adjusted_df_M[, colnames(adjusted_df_M) %in% filtered_samples_M_KAI_IH_CON_A_AC]

coldata_df_M_KAI_IH_CON_A_IM <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IH_CON_A_IM")
filtered_samples_M_KAI_IH_CON_A_IM <- rownames(coldata_df_M_KAI_IH_CON_A_IM)
adjusted_df_M_KAI_IH_CON_A_IM <- adjusted_df_M[, colnames(adjusted_df_M) %in% filtered_samples_M_KAI_IH_CON_A_IM]

coldata_df_M_KAI_IH_CON_A_CR <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IH_CON_A_CR")
filtered_samples_M_KAI_IH_CON_A_CR <- rownames(coldata_df_M_KAI_IH_CON_A_CR)
adjusted_df_M_KAI_IH_CON_A_CR <- adjusted_df_M[, colnames(adjusted_df_M) %in% filtered_samples_M_KAI_IH_CON_A_CR]

coldata_df_M_KAI_IA_IPSI_A_AC <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IA_IPSI_A_AC")
filtered_samples_M_KAI_IA_IPSI_A_AC <- rownames(coldata_df_M_KAI_IA_IPSI_A_AC)
adjusted_df_M_KAI_IA_IPSI_A_AC <- adjusted_df_M[, colnames(adjusted_df_M) %in% filtered_samples_M_KAI_IA_IPSI_A_AC]

coldata_df_M_KAI_IA_IPSI_A_IM <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IA_IPSI_A_IM")
filtered_samples_M_KAI_IA_IPSI_A_IM <- rownames(coldata_df_M_KAI_IA_IPSI_A_IM)
adjusted_df_M_KAI_IA_IPSI_A_IM <- adjusted_df_M[, colnames(adjusted_df_M) %in% filtered_samples_M_KAI_IA_IPSI_A_IM]

coldata_df_M_KAI_IA_IPSI_A_CR <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IA_IPSI_A_CR")
filtered_samples_M_KAI_IA_IPSI_A_CR <- rownames(coldata_df_M_KAI_IA_IPSI_A_CR)
adjusted_df_M_KAI_IA_IPSI_A_CR <- adjusted_df_M[, colnames(adjusted_df_M) %in% filtered_samples_M_KAI_IA_IPSI_A_CR]

coldata_df_M_KAI_IP_A_HA <- coldata_df_M %>% filter(TYPE_2 == "M_KAI_IP_A_HA")
filtered_samples_M_KAI_IP_A_HA <- rownames(coldata_df_M_KAI_IP_A_HA)
adjusted_df_M_KAI_IP_A_HA <- adjusted_df_M[, colnames(adjusted_df_M) %in% filtered_samples_M_KAI_IP_A_HA]

coldata_df_M_PILO_IP_A_HA <- coldata_df_M %>% filter(TYPE_2 == "M_PILO_IP_A_HA" | TYPE_2 == "CTL")
filtered_samples_M_PILO_IP_A_HA <- rownames(coldata_df_M_PILO_IP_A_HA)
adjusted_df_M_PILO_IP_A_HA <- adjusted_df_M[, colnames(adjusted_df_M) %in% filtered_samples_M_PILO_IP_A_HA]

coldata_df_M_PILO_IP_A_AC <- coldata_df_M %>% filter(TYPE_2 == "M_PILO_IP_A_AC")
filtered_samples_M_PILO_IP_A_AC <- rownames(coldata_df_M_PILO_IP_A_AC)
adjusted_df_M_PILO_IP_A_AC <- adjusted_df_M[, colnames(adjusted_df_M) %in% filtered_samples_M_PILO_IP_A_AC]

coldata_df_M_PILO_IP_A_IM <- coldata_df_M %>% filter(TYPE_2 == "M_PILO_IP_A_IM" | TYPE_2 == "CTL")
filtered_samples_M_PILO_IP_A_IM <- rownames(coldata_df_M_PILO_IP_A_IM)
adjusted_df_M_PILO_IP_A_IM <- adjusted_df_M[, colnames(adjusted_df_M) %in% filtered_samples_M_PILO_IP_A_IM]

coldata_df_M_PILO_IP_A_CR <- coldata_df_M %>% filter(TYPE_2 == "M_PILO_IP_A_CR" | TYPE_2 == "CTL")
filtered_samples_M_PILO_IP_A_CR <- rownames(coldata_df_M_PILO_IP_A_CR)
adjusted_df_M_PILO_IP_A_CR <- adjusted_df_M[, colnames(adjusted_df_M) %in% filtered_samples_M_PILO_IP_A_CR]

coldata_df_R_KAI_IP_A_CR <- coldata_df_R %>% filter(TYPE_2 == "R_KAI_IP_A_CR")
filtered_samples_R_KAI_IP_A_CR <- rownames(coldata_df_R_KAI_IP_A_CR)
adjusted_df_R_KAI_IP_A_CR <- adjusted_df_R[, colnames(adjusted_df_R) %in% filtered_samples_R_KAI_IP_A_CR]

coldata_df_R_KAI_SUB_I_IM <- coldata_df_R %>% filter(TYPE_2 == "R_KAI_SUB_I_IM")
filtered_samples_R_KAI_SUB_I_IM <- rownames(coldata_df_R_KAI_SUB_I_IM)
adjusted_df_R_KAI_SUB_I_IM <- adjusted_df_R[, colnames(adjusted_df_R) %in% filtered_samples_R_KAI_SUB_I_IM]

coldata_df_R_PILO_IP_A_CR <- coldata_df_R %>% filter(TYPE_2 == "R_PILO_IP_A_CR")
filtered_samples_R_PILO_IP_A_CR <- rownames(coldata_df_R_PILO_IP_A_CR)
adjusted_df_R_PILO_IP_A_CR <- adjusted_df_R[, colnames(adjusted_df_R) %in% filtered_samples_R_PILO_IP_A_CR]

coldata_df_R_PPS_A_AC <- coldata_df_R %>% filter(TYPE_2 == "R_PPS_A_AC" | TYPE_2 == "PPS_CTL")
filtered_samples_R_PPS_A_AC <- rownames(coldata_df_R_PPS_A_AC)
adjusted_df_R_PPS_A_AC <- adjusted_df_R[, colnames(adjusted_df_R) %in% filtered_samples_R_PPS_A_AC]

coldata_df_R_PPS_A_IM <- coldata_df_R %>% filter(TYPE_2 == "R_PPS_A_IM" | TYPE_2 == "PPS_CTL")
filtered_samples_R_PPS_A_IM <- rownames(coldata_df_R_PPS_A_IM)
adjusted_df_R_PPS_A_IM <- adjusted_df_R[, colnames(adjusted_df_R) %in% filtered_samples_R_PPS_A_IM]

coldata_df_R_PPS_A_CR <- coldata_df_R %>% filter(TYPE_2 == "R_PPS_A_CR" | TYPE_2 == "PPS_CTL")
filtered_samples_R_PPS_A_CR <- rownames(coldata_df_R_PPS_A_CR)
adjusted_df_R_PPS_A_CR <- adjusted_df_R[, colnames(adjusted_df_R) %in% filtered_samples_R_PPS_A_CR]

coldata_df_R_PPS_A_DOFS <- coldata_df_R %>% filter(TYPE_2 == "R_PPS_A_DOFS" | TYPE_2 == "PPS_CTL")
filtered_samples_R_PPS_A_DOFS <- rownames(coldata_df_R_PPS_A_DOFS)
adjusted_df_R_PPS_A_DOFS <- adjusted_df_R[, colnames(adjusted_df_R) %in% filtered_samples_R_PPS_A_DOFS]

coldata_df_R_TBI_IPSI_A_CR <- coldata_df_R %>% filter(TYPE_2 == "R_TBI_IPSI_A_CR")
filtered_samples_R_TBI_IPSI_A_CR <- rownames(coldata_df_R_TBI_IPSI_A_CR)
adjusted_df_R_TBI_IPSI_A_CR <- adjusted_df_R[, colnames(adjusted_df_R) %in% filtered_samples_R_TBI_IPSI_A_CR]

coldata_df_R_AMG_IPSI_A_CR <- coldata_df_R %>% filter(TYPE_2 == "R_AMG_IPSI_A_CR")
filtered_samples_R_AMG_IPSI_A_CR <- rownames(coldata_df_R_AMG_IPSI_A_CR)
adjusted_df_R_AMG_IPSI_A_CR <- adjusted_df_R[, colnames(adjusted_df_R) %in% filtered_samples_R_AMG_IPSI_A_CR]
