library(dplyr)
library(ggplot2)
library(ggsankey)
library(readxl)


### use already defined coldata_df_1
coldata_df_1_MTLE <- coldata_df_1 %>%
  filter(Diagnosis != "NL")

## 
coldata_df_1_MTLE <- coldata_df_1_MTLE %>% mutate(ip_group = case_when(
  IP == 'HT(Head Trauma)' ~ "Initial precipitating injury (Other)",
  IP == 'TCPC(Tonic-Clonic Postictal Confusion)' ~ "Initial precipitating injury (Seizure)",
  IP == 'CVA(Stroke)' ~ "Initial precipitating injury (Other)",
  IP == 'FS(Febrile Seizure)' ~ "Initial precipitating injury (Seizure)",
  IP == 'MEN(Meningitis)' ~ "Initial precipitating injury (Other)",
  IP == 'GFS(Generalized Febrile Seizure)' ~ "Initial precipitating injury (Seizure)",
  IP == 'ND(No Data)' ~ "No event",
  TRUE ~ NA_character_
))

#### plot
sankey_df <- coldata_df_1_MTLE %>% 
  mutate(
    Diagnosis = replace_na(as.character(Diagnosis_Sub), "Missing"),
    ip_group  = replace_na(as.character(ip_group), "Missing")
  ) %>% 
  # Convert your data from wide to long format.
  make_long(Diagnosis, ip_group)

node_colors <- c(
  "MTLEHS" = "#B4Eeb4", # darkseagreen2
  "MTLE" = "#CAE1FF",   # lightsteelblue1
  "Missing" = "#EEDFCC",      # antiquewhite2
  "No event" = "#EEDFCC",      # antiquewhite2
  "Initial precipitating injury (Seizure)" = "#EEDFCC",      # antiquewhite2
  "Initial precipitating injury (Other)" = "#EEDFCC"      # antiquewhite2
)

# sankey
plot_width <- 6
plot_height <- 3
plot_unit <- 'in'
plot_dpi <- 300
plot_target <- ggplot(sankey_df, 
                      aes(x = x, 
                          next_x = next_x, 
                          node = node, 
                          next_node = next_node, 
                          fill = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = "grey30") +
  geom_sankey_label(aes(label = node), size = 3, fill = "white") +
  scale_fill_manual(values = node_colors) +  # custom color scale
  theme_sankey(base_size = 14)
file_name <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/11.Sample_variables/00.sankey_plot/00.Sample_variable_sankey.png"
ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)
