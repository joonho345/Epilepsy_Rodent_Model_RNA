library(dplyr)

#### choose groups to compare #### 
## hippo: only hippo
coldata_df_hippo <- coldata_df %>%
  filter(Brain_Location == "HIPPO") %>%
  mutate(across(everything(), ~ replace(.x, .x == "", "N/A")))
filtered_samples_hippo <- rownames(coldata_df_hippo)
adjusted_df_hippo <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples_hippo]
# c("MTLEHS", "MTLE", "NL", "SCZD", "MDD", "BPD", "AUD")

########################
## A: MTLEHS & NL
coldata_df_A <- coldata_df %>% 
  filter(Diagnosis == "MTLEHS" | Diagnosis == "NL") %>% 
  mutate(across(everything(), ~ replace(.x, .x == "", "N/A")))
filtered_samples_A <- rownames(coldata_df_A)
adjusted_df_A <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples_A]

## B: MTLEHS & NL (HIPPO ONLY)
coldata_df_B <- coldata_df %>% 
  filter(Diagnosis == "MTLEHS" | Diagnosis == "NL") %>% 
  filter(Brain_Location == "HIPPO") %>% 
  mutate(across(everything(), ~ replace(.x, .x == "", "N/A")))
filtered_samples_B <- rownames(coldata_df_B)
adjusted_df_B <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples_B]

## C: MTLEHS & NL (HIPPO & WHOLE ONLY)
coldata_df_C <- coldata_df %>% 
  filter(Diagnosis == "MTLEHS" | Diagnosis == "NL") %>% 
  filter(Brain_Location == "HIPPO") %>% 
  filter(Brain_Location_Sub == "Whole") %>%
  mutate(across(everything(), ~ replace(.x, .x == "", "N/A")))
filtered_samples_C <- rownames(coldata_df_C)
adjusted_df_C <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples_C]

## D: MTLEHS & OtherEpilepsy
coldata_df_D <- coldata_df[coldata_df$Diagnosis %in% c("MTLEHS", "MTLE", "FCD", "TSC"), ]
coldata_df_D$Diagnosis_D <- ifelse(coldata_df_D$Diagnosis == "MTLEHS", "MTLEHS", "OtherEpilepsy")
coldata_df_D$Diagnosis <- NULL
colnames(coldata_df_D)[which(names(coldata_df_D) == "Diagnosis_D")] <- "Diagnosis"
filtered_samples_D <- rownames(coldata_df_D)
adjusted_df_D <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples_D]

## E: MTLEALL & NL
coldata_df_E <- coldata_df[coldata_df$Diagnosis %in% c("MTLEHS", "MTLE", "NL"), ]
coldata_df_E$Diagnosis_E <- ifelse(coldata_df_E$Diagnosis %in% c("MTLEHS", "MTLE"), "MTLEALL", "NL")
coldata_df_E$Diagnosis <- NULL
colnames(coldata_df_E)[which(names(coldata_df_E) == "Diagnosis_E")] <- "Diagnosis"
filtered_samples_E <- rownames(coldata_df_E)
adjusted_df_E <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples_E]

## F: MTLEALL & NL (HIPPO ONLY)
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

## G: MTLEALL & OtherEpilepsy
coldata_df_G <- coldata_df[coldata_df$Diagnosis %in% c("MTLEHS", "MTLE", "FCD", "TSC"), ]
coldata_df_G$Diagnosis_G <- ifelse(coldata_df_G$Diagnosis %in% c("MTLEHS", "MTLE"), "MTLEALL", "OtherEpilepsy")
coldata_df_G$Diagnosis <- NULL
colnames(coldata_df_G)[which(names(coldata_df_G) == "Diagnosis_G")] <- "Diagnosis"
filtered_samples_G <- rownames(coldata_df_G)
adjusted_df_G <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples_G]

## H: MTLE & NL (HIPPO)
coldata_df_H <- coldata_df %>% 
  filter(Diagnosis == "MTLE" | Diagnosis == "NL") %>%
  filter(Brain_Location == "HIPPO") %>% 
  mutate(across(everything(), ~ replace(.x, .x == "", "N/A")))
filtered_samples_H <- rownames(coldata_df_H)
adjusted_df_H <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples_H]

## I: MTLE & NL
coldata_df_I <- coldata_df %>% 
  filter(Diagnosis == "MTLE" | Diagnosis == "NL") %>%
  mutate(across(everything(), ~ replace(.x, .x == "", "N/A")))
filtered_samples_I <- rownames(coldata_df_I)
adjusted_df_I <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples_I]


############################
Target_comparison <- 'A_MTLEHS_NL'
Target_comparison <- 'B_MTLEHS_NL_HIPPO'
Target_comparison <- 'C_MTLEHS_NL_HIPPO_W'
Target_comparison <- 'D_MTLEHS_OtherEpi'
Target_comparison <- 'E_MTLEALL_NL'
Target_comparison <- 'F_MTLEALL_NL_HIPPO'
Target_comparison <- 'G_MTLEALL_OtherEpi'
Target_comparison <- 'H_MTLE_NL_HIPPO'
Target_comparison <- 'I_MTLE_NL'

group_mapping <- list(
  'A_MTLEHS_NL' = c("MTLEHS", "NL"),
  'B_MTLEHS_NL_HIPPO' = c("MTLEHS", "NL"),
  'C_MTLEHS_NL_HIPPO_W' = c("MTLEHS", "NL"),
  'D_MTLEHS_OtherEpi' = c("MTLEHS", "OtherEpilepsy"),
  'E_MTLEALL_NL' = c("MTLEALL", "NL"),
  'F_MTLEALL_NL_HIPPO' = c("MTLEALL", "NL"),
  'G_MTLEALL_OtherEpi' = c("MTLEALL", "OtherEpilepsy"),
  'H_MTLE_NL_HIPPO' = c("MTLE", "NL"),
  'I_MTLE_NL' = c("MTLE", "NL")
)

data_mapping <- list(
  'A_MTLEHS_NL' = list(Target_df = adjusted_df_A, Target_coldata = coldata_df_A),
  'B_MTLEHS_NL_HIPPO' = list(Target_df = adjusted_df_B, Target_coldata = coldata_df_B),
  'C_MTLEHS_NL_HIPPO_W' = list(Target_df = adjusted_df_C, Target_coldata = coldata_df_C),
  'D_MTLEHS_OtherEpi' = list(Target_df = adjusted_df_D, Target_coldata = coldata_df_D),
  'E_MTLEALL_NL' = list(Target_df = adjusted_df_E, Target_coldata = coldata_df_E),
  'F_MTLEALL_NL_HIPPO' = list(Target_df = adjusted_df_F, Target_coldata = coldata_df_F),
  'G_MTLEALL_OtherEpi' = list(Target_df = adjusted_df_G, Target_coldata = coldata_df_G),
  'H_MTLE_NL_HIPPO' = list(Target_df = adjusted_df_H, Target_coldata = coldata_df_H),
  'I_MTLE_NL' = list(Target_df = adjusted_df_I, Target_coldata = coldata_df_I)
)



################################
#### AFTER Sample Selection ####
################################
OUTLIER_samples <- c(
  "SRR15406573", "SRR28007131", "SRR28007092", "SRR8669936",  "SRR16522693", "SRR28007020")
#OUTLIER_samples <- c(
#  "SRR15406573", 
#  "SRR28007131", "SRR28007018", "SRR28007079",
#  "SRR16522705", "SRR28007092", "SRR28007081", "SRR28007133", "SRR16522707",
#  "SRR16522708", "SRR28007072", "SRR16522706", "SRR8669931", "SRR28007122",
#  "SRR28007034", "SRR28007155", "SRR16522692", "SRR28007132", "SRR16522690",
#  "SRR8669932", "SRR8669936", "SRR16522693", "SRR16522689", "SRR28007020",
#  "SRR16522691", "SRR16522694", "SRR28007089", "SRR8669933"
#)
#OUTLIER_samples <- c(
#  "SRR5241783", "SRR9733960", "SRR9733958", "SRR15406573", "SRR8669937",
#  "SRR5241781", "SRR8669940", "SRR28007131", "SRR28007018", "SRR28007079",
#  "SRR16522705", "SRR28007092", "SRR28007081", "SRR28007133", "SRR16522707",
#  "SRR16522708", "SRR28007072", "SRR16522706", "SRR8669931", "SRR28007122",
#  "SRR28007034", "SRR28007155", "SRR16522692", "SRR28007132", "SRR16522690",
#  "SRR8669932", "SRR8669936", "SRR16522693", "SRR16522689", "SRR28007020",
#  "SRR16522691", "SRR16522694", "SRR28007089", "SRR8669933"
#)

## F: MTLEALL & NL (HIPPO ONLY)
## FILTERED_1_MTLEALL_NL
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

coldata_df_1 <- coldata_df_F %>%
  filter(!(rownames(coldata_df_F) %in% OUTLIER_samples))
filtered_samples_F_no_outliers <- rownames(coldata_df_1)
adjusted_df_1 <- adjusted_df_F[, colnames(adjusted_df_F) %in% filtered_samples_F_no_outliers]

## B: MTLEHS & NL (HIPPO ONLY)
## FILTERED_2_MTLEHS_NL
coldata_df_B <- coldata_df %>% 
  filter(Diagnosis == "MTLEHS" | Diagnosis == "NL") %>% 
  filter(Brain_Location == "HIPPO") %>% 
  mutate(across(everything(), ~ replace(.x, .x == "", "N/A")))
filtered_samples_B <- rownames(coldata_df_B)
adjusted_df_B <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples_B]

coldata_df_2 <- coldata_df_B %>%
  filter(!(rownames(coldata_df_B) %in% OUTLIER_samples))
filtered_samples_B_no_outliers <- rownames(coldata_df_2)
adjusted_df_2 <- adjusted_df_B[, colnames(adjusted_df_B) %in% filtered_samples_B_no_outliers]

## H: MTLE & NL (HIPPO)
## FILTERED_3_MTLE_NL
coldata_df_H <- coldata_df %>% 
  filter(Diagnosis == "MTLE" | Diagnosis == "NL") %>%
  filter(Brain_Location == "HIPPO") %>% 
  mutate(across(everything(), ~ replace(.x, .x == "", "N/A")))
filtered_samples_H <- rownames(coldata_df_H)
adjusted_df_H <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples_H]

coldata_df_3 <- coldata_df_H %>%
  filter(!(rownames(coldata_df_H) %in% OUTLIER_samples))
filtered_samples_H_no_outliers <- rownames(coldata_df_3)
adjusted_df_3 <- adjusted_df_H[, colnames(adjusted_df_H) %in% filtered_samples_H_no_outliers]

############################
Target_comparison <- 'FILTERED_1_MTLEALL_NL'
Target_comparison <- 'FILTERED_2_MTLEHS_NL'
Target_comparison <- 'FILTERED_3_MTLE_NL'

group_mapping <- list(
  'FILTERED_1_MTLEALL_NL' = c("MTLEALL", "NL"),
  'FILTERED_2_MTLEHS_NL' = c("MTLEHS", "NL"),
  'FILTERED_3_MTLE_NL' = c("MTLE", "NL")
)

data_mapping <- list(
  'FILTERED_1_MTLEALL_NL' = list(Target_df = adjusted_df_1, Target_coldata = coldata_df_1),
  'FILTERED_2_MTLEHS_NL' = list(Target_df = adjusted_df_2, Target_coldata = coldata_df_2),
  'FILTERED_3_MTLE_NL' = list(Target_df = adjusted_df_3, Target_coldata = coldata_df_3)
)

