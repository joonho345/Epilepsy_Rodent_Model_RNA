library(dplyr)


### use already defined adjusted_df_1 (raw)
### use already defined coldata_df_1

## 
coldata_df_1_MTLE_NL <- coldata_df_1_MTLE_NL %>%   mutate(ip_group = case_when(
  Diagnosis != "NL" & IP %in% c('HT(Head Trauma)', 'TCPC(Tonic-Clonic Postictal Confusion)',
                                'CVA(Stroke)', 'FS(Febrile Seizure)',
                                'MEN(Meningitis)', 'GFS(Generalized Febrile Seizure)') ~ "Initial precipitating injury",
  Diagnosis != "NL" & IP == 'ND(No Data)' ~ "No event",
  TRUE ~ NA_character_
))

#grouping
coldata_df_HS <- coldata_df_1_MTLE_NL %>%
  filter(Diagnosis_Sub == "MTLEHS" | Diagnosis == "NL")
filtered_samples <- rownames(coldata_df_HS)
adjusted_df_HS <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples]

coldata_df_withoutHS <- coldata_df_1_MTLE_NL %>%
  filter(Diagnosis_Sub == "MTLE" | Diagnosis == "NL")
filtered_samples <- rownames(coldata_df_withoutHS)
adjusted_df_withoutHS <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples]

coldata_df_IP <- coldata_df_1_MTLE_NL %>%
  filter((Diagnosis == "MTLEALL" & ip_group == "Initial precipitating injury") | Diagnosis == "NL")
filtered_samples <- rownames(coldata_df_IP)
adjusted_df_IP <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples]

coldata_df_NO <- coldata_df_1_MTLE_NL %>%
  filter((Diagnosis == "MTLEALL" & ip_group == "No event") | Diagnosis == "NL")
filtered_samples <- rownames(coldata_df_NO)
adjusted_df_NO <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples]

coldata_df_IP_Sei <- coldata_df_1_MTLE_NL %>%
  filter((Diagnosis == "MTLEALL" & IP %in% c('TCPC(Tonic-Clonic Postictal Confusion)', 'FS(Febrile Seizure)', 
                                             'GFS(Generalized Febrile Seizure)')) | Diagnosis == "NL")
filtered_samples <- rownames(coldata_df_IP_Sei)
adjusted_df_IP_Sei <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples]

coldata_df_IP_Other <- coldata_df_1_MTLE_NL %>%
  filter((Diagnosis == "MTLEALL" & IP %in% c('HT(Head Trauma)', 'CVA(Stroke)', 'MEN(Meningitis)')) 
         | Diagnosis == "NL")
filtered_samples <- rownames(coldata_df_IP_Other)
adjusted_df_IP_Other <- adjusted_df[, colnames(adjusted_df) %in% filtered_samples]


