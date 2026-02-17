library(ggplot2)

file_path <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/08.Enrichment/GSEA_03.Output/gsea_report_for_MTLEALL_GOBP_1000_new_revision.tsv"
name <- "_1000_new_revision"

df <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
# df$NOM.p.val[df$NOM.p.val == 0] <- 0.01
df$log_pval <- -log10(df$NOM.p.val)

df <- df[!grepl("^(ERK|CC|TLR|NR)", df$NAME), ]
df$Phase <- ifelse(grepl("^(NT|IC)_", df$NAME), "Acute",
                   ifelse(grepl("^(NI|ND|TNF)_", df$NAME), "Subacute",
                          ifelse(grepl("^(SP|GG|NG|MF)_", df$NAME), "Chronic", NA)))

df$Phase_sub <- ifelse(grepl("^NT_", df$NAME), "NT",
                       ifelse(grepl("^IC_", df$NAME), "IC",
                              ifelse(grepl("^CC_", df$NAME), "CC",
                                     ifelse(grepl("^NI_", df$NAME), "NI",
                                            ifelse(grepl("^ND_", df$NAME), "ND",
                                                   ifelse(grepl("^TLR_", df$NAME), "TLR",
                                                          ifelse(grepl("^TNF_", df$NAME), "TNF",
                                                                 ifelse(grepl("^NR_", df$NAME), "NR",
                                                                        ifelse(grepl("^NG_", df$NAME), "NG",
                                                                               ifelse(grepl("^GG_", df$NAME), "GG",
                                                                                      ifelse(grepl("^MF_", df$NAME), "MF",
                                                                                             ifelse(grepl("^SP_", df$NAME), "SP", NA))))))))))))

# Create a lookup table for GO_NAME
go_mapping <- data.frame(
  NAME = c(
    "GG_GO0042063_410", "GG_GO0014015_94", "GG_GO0014004_16", "GG_GO0048708_94", "GG_GO0014014_59",
    "NG_GO0050769_331", "NG_GO0002052_44", "NG_GO0050768_163", "NG_GO0007406_16", "NG_GO0050772_85", "NG_GO0050771_56",
    "SP_GO0048168_78", "SP_GO0048169_39", "SP_GO0048172_22",
    "NT_GO0051932_73", "NT_GO0035249_131", "NT_GO0008328_44", "NT_GO0032281_29", "NT_GO0032983_6", "NT_GO0017146_9",
    "IC_GO0005245_43", "IC_GO0005249_98", "IC_GO0005248_24",
    "NI_GO0006954_882", "NI_GO0150076_87", "NI_GO0048143_20", "NI_GO0001774_49", "NI_GO0050729_153", "NI_GO0050728_173",
    "ND_GO0051402_388", "ND_GO0110088_8", "ND_GO0043525_105", "ND_GO0043524_214",
    "MF_GO0098686_61", "MF_GO0048668_35",
    "TNF_GO0033209_103", "TNF_GO1903265_11", "TNF_GO0010804_27"
  ),
  GO_NAME = c(
    "Gliogenesis", "Positive Regulation of Gliogenesis", "Microglia Differentiation", "Astrocyte Differentiation", "Negative Regulation of Gliogenesis",
    "Positive Regulation of Neurogenesis", "Positive Regulation of Neuroblast Proliferation", "Negative Regulation of Neurogenesis", "Negative Regulation of Neuroblast Proliferation",
    "Positive Regulation of Axonogenesis", "Negative Regulation of Axonogenesis",
    "Regulation of Neuronal Synaptic Plasticity", "Regulation of Long-term Neuronal Synaptic Plasticity", "Regulation of Short-term Neuronal Synaptic Plasticity",
    "Synaptic Transmission, GABAergic", "Synaptic Transmission, Glutamatergic", "Ionotropic Glutamate Receptor Complex", "AMPA Glutamate Receptor Complex", 
    "Kainate Selective Glutamate Receptor Complex", "NMDA Selective Glutamate Receptor Complex",
    "Voltage-gated Calcium Channel Activity", "Voltage-gated Potassium Channel Activity", "Voltage-gated Sodium Channel Activity",
    "Inflammatory Response", "Neuroinflammatory Response", "Astrocyte Activation", "Microglial Cell Activation", "Positive Regulation of Inflammatory Response", 
    "Negative Regulation of Inflammatory Response",
    "Neuron Apoptotic Process", "Hippocampal Neuron Apoptotic Process", "Positive Regulation of Neuron Apoptotic Process", "Negative Regulation of Neuron Apoptotic Process",
    "Hippocampal Mossy Fiber to CA3 Synapse", "Collateral Sprouting",
    "TNF-mediated Signaling Pathway", "Positive Regulation of TNF-mediated Signaling Pathway", "Negative Regulation of TNF-mediated Signaling Pathway"
  ),
  stringsAsFactors = FALSE
)

df <- merge(df, go_mapping, by = "NAME", all.x = TRUE)
phase_order <- c("NT", "IC", "NI", "ND", "TNF", "NG", "GG", "MF", "SP")
df <- df[order(factor(df$Phase_sub, levels = phase_order), -df$NES), ]
df$GO_NAME <- factor(df$GO_NAME, levels = rev(unique(df$GO_NAME)))  # Reverse for top-to-bottom plotting
df <- df[grepl("_", df$NAME), ]

phase_colors <- c(
  "NT" = "#1a9850",
  "IC" = "#91cf60",
  "NI" = "#313695",
  "ND" = "#74add1",
  "TNF" = "#4575b4",
  "NG" = "#d73027",
  "GG" = "#f46d43",
  "MF" = "#fee08b",
  "SP" = "#fdae61"
)
df$ID_color <- phase_colors[df$Phase_sub]
df_filtered <- df %>% filter(!is.na(NOM.p.val))

df_filtered$Phase_sub_c <- ifelse(df_filtered$NES > 0, 
                                  paste0(df_filtered$Phase_sub, "_MTLEALL"), 
                                  paste0(df_filtered$Phase_sub, "_NL"))
df_filtered$Group <- ifelse(df_filtered$NES > 0, "MTLEALL", "NL")

######## combined p-value #########
# fisher_method
fisher_method <- function(p_values) {
  # Replace 0 or very small p-values with machine minimum to avoid log(0) = -Inf
  p_values <- pmax(p_values, .Machine$double.xmin)
  chi_sq <- -2 * sum(log(p_values))
  df <- 2 * length(p_values)  # Degrees of freedom
  p_combined <- pchisq(chi_sq, df, lower.tail = FALSE)  # Get combined p-value
  # If the result is too small, return the minimum representable value
  if (p_combined == 0 || is.na(p_combined)) {
    p_combined <- .Machine$double.xmin
  }
  return(p_combined)
}
combined_p_values <- df_filtered %>%
  group_by(Phase_sub_c) %>%
  summarise(
    combined_p = fisher_method(NOM.p.val)
  )

# Plotting the data
go_colors <- c(
  "NT" = "#1a9850",
  "IC" = "#91cf60",
  "NI" = "#313695",
  "ND" = "#74add1",
  "TNF" = "#4575b4",
  "NG" = "#d73027",
  "GG" = "#f46d43",
  "MF" = "#fee08b",
  "SP" = "#fdae61"
)

# Cap log_combined_p at maximum for visualization (to handle extremely small p-values
max_log_pval <- 5.0  # Cap at 5 (p < 1e-5)
combined_p_values <- combined_p_values %>%
  mutate(
    log_combined_p_raw = ifelse(grepl("MTLEALL", Phase_sub_c), -log10(combined_p), log10(combined_p)),
    log_combined_p = ifelse(grepl("MTLEALL", Phase_sub_c), 
                            pmin(log_combined_p_raw, max_log_pval),  # Cap positive values (MTLEALL)
                            pmax(log_combined_p_raw, -max_log_pval)), # Cap negative values (NL)
    Phase = factor(gsub("_MTLEALL|_NL", "", Phase_sub_c), levels = c("NT", "IC", "NI", "ND", "TNF", "NG", "GG", "MF", "SP")),
    alpha_value = ifelse(abs(log_combined_p_raw) >= -log10(0.05), 1, 0.5)  # Adjust alpha based on significance (use raw for alpha)
  )

axis_text_size <- 8
legend_text_size <- 8
title_text_size <- 12
title_face <- "bold"
font_family <- "Arial"
plot_target <- ggplot(combined_p_values, aes(x = Phase, y = log_combined_p, color = Phase, alpha = alpha_value)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -log10(0.05), ymax = log10(0.05)), fill = "grey96", inherit.aes = FALSE) +
  
  geom_segment(aes(x = Phase, xend = Phase, y = 0, yend = log_combined_p), color = "grey50", size = 1) +
  geom_hline(yintercept = c(-log10(0.05), log10(0.05)), linetype = "dashed", color = "gray", size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.7) +
  geom_point(size = 4) +
  
  scale_color_manual(values = go_colors) +
  scale_alpha_identity() +
  scale_y_continuous(
    breaks = c(-2, -1, 0, 1, 2, 3, 4, 5), 
    labels = c("-2", "-1", "0", "1", "2", "3", "4", "≥5") 
  ) +
  labs(y = "-Log10 Combined p-value") +
  
  theme_minimal() +
  theme(
    axis.title = element_text(size = title_text_size, family = font_family, face = title_face),
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.text.y = element_text(size = axis_text_size, family = font_family), 
    legend.position = "none",
    axis.line.y = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  ylim(-2.5, 5.5)

# Save the plot
plot_width <- 4
plot_height <- 4
plot_unit <- 'in'
plot_dpi <- 600
file_name <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/08.Enrichment/GSEA_03.Output/GOBP_MTLEALL_NL",name,"_combined.png")
ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)
