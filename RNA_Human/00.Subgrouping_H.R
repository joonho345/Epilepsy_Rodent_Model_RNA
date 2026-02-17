library(dplyr)

################################
#### BEFORE Sample Selection ####
################################

#### choose groups to compare #### 
## hippo: only hippo
coldata_df_hippo <- coldata_df %>%
  filter(Brain_Location == "HIPPO") %>%
  mutate(across(everything(), ~ replace(.x, .x == "", "N/A")))
filtered_samples_hippo <- rownames(coldata_df_hippo)
adjusted_df_hippo <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples_hippo]
# c("MTLEHS", "MTLE", "NL", "SCZD", "MDD", "BPD", "AUD")


################################
#### AFTER Sample Selection ####
################################
OUTLIER_samples <- c("SRR10867969", "SRR10868146", "SRR16522693SRR16522703", "SRR16522694SRR16522704", "SRR28007092", "SRR3628384", "SRR10868129", "SRR15406573", "SRR10868083", "SRR16522692SRR16522702")

## FILTERED_1_MTLEALL_NL - MTLEALL & NL (HIPPO ONLY)
# Ensure sample alignment between coldata and expression matrix
coldata_df_F <- coldata_df[coldata_df$Diagnosis %in% c("MTLEHS", "MTLE", "NL"), ]
coldata_df_F$Diagnosis_F <- ifelse(coldata_df_F$Diagnosis %in% c("MTLEHS", "MTLE"), "MTLEALL", "NL")
coldata_df_F$Diagnosis <- NULL
colnames(coldata_df_F)[which(names(coldata_df_F) == "Diagnosis_F")] <- "Diagnosis"
coldata_df_F <- coldata_df_F %>% 
  filter(Brain_Location == "HIPPO") %>% 
  mutate(across(everything(), ~ replace(.x, .x == "", "N/A"))) %>%
  mutate(Diagnosis_Sub = ifelse(Diagnosis_Sub == "TLE", "MTLE", Diagnosis_Sub))
filtered_samples_F <- rownames(coldata_df_F)
adjusted_df_F <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples_F]
# Filter out outliers
coldata_df_1 <- coldata_df_F %>%
  filter(!(rownames(coldata_df_F) %in% OUTLIER_samples))
filtered_samples_F_no_outliers <- rownames(coldata_df_1)
adjusted_df_1 <- adjusted_df_F[, filtered_samples_F_no_outliers, drop = FALSE]
# Ensure final alignment
if (!all(colnames(adjusted_df_1) == rownames(coldata_df_1))) {
  coldata_df_1 <- coldata_df_1[colnames(adjusted_df_1), , drop = FALSE]
}

## FILTERED_1_MTLEALL_NL_W - MTLEALL & NL (HIPPO & WHOLE ONLY)
# Ensure sample alignment between coldata and expression matrix
coldata_df_J <- coldata_df[coldata_df$Diagnosis %in% c("MTLEHS", "MTLE", "NL"), ]
coldata_df_J$Diagnosis_J <- ifelse(coldata_df_J$Diagnosis %in% c("MTLEHS", "MTLE"), "MTLEALL", "NL")
coldata_df_J$Diagnosis <- NULL
colnames(coldata_df_J)[which(names(coldata_df_J) == "Diagnosis_J")] <- "Diagnosis"
coldata_df_J <- coldata_df_J %>% 
  filter(Brain_Location == "HIPPO") %>% 
  filter(Brain_Location_Sub == "Whole") %>%
  mutate(across(everything(), ~ replace(.x, .x == "", "N/A")))
filtered_samples_J <- rownames(coldata_df_J)
adjusted_df_J <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples_J]
# Filter out outliers
coldata_df_1_W <- coldata_df_J %>%
  filter(!(rownames(coldata_df_J) %in% OUTLIER_samples))
filtered_samples_J_no_outliers_W <- rownames(coldata_df_1_W)
adjusted_df_1_W <- adjusted_df_J[, filtered_samples_J_no_outliers_W, drop = FALSE]
# Ensure final alignment
if (!all(colnames(adjusted_df_1_W) == rownames(coldata_df_1_W))) {
  coldata_df_1_W <- coldata_df_1_W[colnames(adjusted_df_1_W), , drop = FALSE]
}

## FILTERED_2_MTLEHS_NL - MTLEHS & NL (HIPPO ONLY)
# Ensure sample alignment between coldata and expression matrix
coldata_df_B <- coldata_df %>% 
  filter(Diagnosis == "MTLEHS" | Diagnosis == "NL") %>% 
  filter(Brain_Location == "HIPPO") %>% 
  mutate(across(everything(), ~ replace(.x, .x == "", "N/A")))
filtered_samples_B <- rownames(coldata_df_B)
adjusted_df_B <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples_B]
# Filter out outliers
coldata_df_2 <- coldata_df_B %>%
  filter(!(rownames(coldata_df_B) %in% OUTLIER_samples))
filtered_samples_B_no_outliers <- rownames(coldata_df_2)
adjusted_df_2 <- adjusted_df_B[, filtered_samples_B_no_outliers, drop = FALSE]
# Ensure final alignment
if (!all(colnames(adjusted_df_2) == rownames(coldata_df_2))) {
  coldata_df_2 <- coldata_df_2[colnames(adjusted_df_2), , drop = FALSE]
}

## FILTERED_3_MTLE_NL - MTLE & NL (HIPPO ONLY)
# Ensure sample alignment between coldata and expression matrix
coldata_df_H <- coldata_df %>% 
  filter(Diagnosis == "MTLE" | Diagnosis == "NL") %>%
  filter(Brain_Location == "HIPPO") %>% 
  mutate(across(everything(), ~ replace(.x, .x == "", "N/A")))
filtered_samples_H <- rownames(coldata_df_H)
adjusted_df_H <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples_H]
# Filter out outliers
coldata_df_3 <- coldata_df_H %>%
  filter(!(rownames(coldata_df_H) %in% OUTLIER_samples))
filtered_samples_H_no_outliers <- rownames(coldata_df_3)
adjusted_df_3 <- adjusted_df_H[, filtered_samples_H_no_outliers, drop = FALSE]
# Ensure final alignment
if (!all(colnames(adjusted_df_3) == rownames(coldata_df_3))) {
  coldata_df_3 <- coldata_df_3[colnames(adjusted_df_3), , drop = FALSE]
}
