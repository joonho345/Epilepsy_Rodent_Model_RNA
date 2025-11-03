# Load required libraries
library(ggplot2)
library(dplyr)

# Define GO terms and colors
GO_to_analyze <- c("GG", "NG", "SP", "MF", "IC", "NT", "NI", "ND", "TNF")
GO_colors <- c(
  "NT" = "#1a9850",   # Deep green
  "IC" = "#91cf60",   # Vibrant green
  
  "NI" = "#313695",   # Deep indigo
  "ND" = "#74add1",   # Light blue
  "TNF" = "#4575b4",   # Rich blue
  
  "NG" = "#d73027",   # Bold red
  "GG" = "#f46d43",   # Distinct orange
  "MF" = "#fee08b",   # Light yellow-orange
  "SP" = "#fdae61"   # Bright orange-yellow
)

# Define file paths
file_paths <- list(
  MTLEALL = "/home/joonho345/3_RNA/RNA_Animal/07.Enrichment/GSEA_03.Output/gsea_report_for_MTLEALL_GOBP_1000_new.tsv",
  M_KAI_ST_BO_O_HA = "/home/joonho345/3_RNA/RNA_Animal/07.Enrichment/GSEA_03.Output/M_KAI_ST_BO_O_HA.tsv",
  M_KAI_ST_BO_O_AC = "/home/joonho345/3_RNA/RNA_Animal/07.Enrichment/GSEA_03.Output/M_KAI_ST_BO_O_AC.tsv",
  M_KAI_ST_BO_O_IM = "/home/joonho345/3_RNA/RNA_Animal/07.Enrichment/GSEA_03.Output/M_KAI_ST_BO_O_IM.tsv",
  M_KAI_ST_BO_O_CR = "/home/joonho345/3_RNA/RNA_Animal/07.Enrichment/GSEA_03.Output/M_KAI_ST_BO_O_CR.tsv"
)
Phase_level <- c("HA", "AC", "IM", "CR", "MTLEALL")
Treatment_type <- 'KAI_ST'

file_paths <- list(
  MTLEALL = "/home/joonho345/3_RNA/RNA_Animal/07.Enrichment/GSEA_03.Output/gsea_report_for_MTLEALL_GOBP_1000_new.tsv",
  M_PILO_Y_HA = "/home/joonho345/3_RNA/RNA_Animal/07.Enrichment/GSEA_03.Output/M_PILO_Y_HA.tsv",
  M_PILO_Y_AC = "/home/joonho345/3_RNA/RNA_Animal/07.Enrichment/GSEA_03.Output/M_PILO_Y_AC.tsv",
  M_PILO_Y_IM = "/home/joonho345/3_RNA/RNA_Animal/07.Enrichment/GSEA_03.Output/M_PILO_Y_IM.tsv"
)
Treatment_type <- 'PILO'
Phase_level <- c("HA", "AC", "IM", "MTLEALL")

file_paths <- list(
  MTLEALL = "/home/joonho345/3_RNA/RNA_Animal/07.Enrichment/GSEA_03.Output/gsea_report_for_MTLEALL_GOBP_1000_new.tsv",
  R_PPS_O_AC = "/home/joonho345/3_RNA/RNA_Animal/07.Enrichment/GSEA_03.Output/R_PPS_O_AC.tsv",
  R_PPS_O_IM = "/home/joonho345/3_RNA/RNA_Animal/07.Enrichment/GSEA_03.Output/R_PPS_O_IM.tsv",
  R_PPS_O_DOFS = "/home/joonho345/3_RNA/RNA_Animal/07.Enrichment/GSEA_03.Output/R_PPS_O_DOFS.tsv",
  R_PPS_O_CR = "/home/joonho345/3_RNA/RNA_Animal/07.Enrichment/GSEA_03.Output/R_PPS_O_CR.tsv"
)
Treatment_type <- 'PPS'
Phase_level <- c("AC", "DOFS", "IM", "CR", "MTLEALL")

# Initialize an empty list to store NES values
filtered_dfs <- list()

for (GO_name in GO_to_analyze) {
  for (id in names(file_paths)) {
    file_path <- file_paths[[id]]
    tsv_data <- read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "\"")
    go_data <- tsv_data %>% filter(NAME == GO_name)
    if (nrow(go_data) > 0) {
      go_data$Phase <- ifelse(grepl("HA", id), "HA",
                              ifelse(grepl("AC", id), "AC",
                                     ifelse(grepl("IM", id), "IM",
                                            ifelse(grepl("DOFS", id), "DOFS",
                                                   ifelse(grepl("CR", id), "CR", "MTLEALL")))))
      go_data$GO <- GO_name
      filtered_dfs[[paste(GO_name, id, sep = "_")]] <- go_data}}}

nes_df <- do.call(rbind, filtered_dfs)
nes_df$Phase <- factor(nes_df$Phase, levels = Phase_level)
legend_order <- c("NT", "IC", "NI", "ND", "TNF", "NG", "GG", "MF", "SP")
nes_df$GO <- factor(nes_df$GO, levels = legend_order)

GO_legend_labels <- c(
  "NT" = "Neurotransmission",
  "IC" = "Ion Channel",
  "NI" = "Neuroinflammation",
  "ND" = "Neuronal Death",
  "TNF" = "TNF Signaling",
  "NG" = "Neurogenesis",
  "GG" = "Gliogenesis",
  "MF" = "Mossy Fiber Sprouting",
  "SP" = "Synaptic Plasticity"
)
nes_df$GO <- factor(nes_df$GO, levels = names(GO_legend_labels))

# Create dot plot for all GO terms with reordered legend
axis_text_size <- 8
legend_text_size <- 8
title_text_size <- 12
title_face <- "bold"
font_family <- "Arial"
plot_target <- ggplot(nes_df, aes(x = Phase, y = NES, color = GO, group = GO)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", size = 0.6) +
  geom_hline(yintercept = c(-2, -1, 1, 2), linetype = "dashed", color = "gray", size = 0.2) +
  geom_vline(xintercept = which(levels(nes_df$Phase) == "MTLEALL") - 0.5, linetype = "dotted", color = "red", size = 0.8) +  
  geom_line(size = 1, alpha = 0.5) +  # Connect dots with lines
  geom_point(size = 3, alpha = 1) +  # Add dots
  scale_color_manual(values = GO_colors, labels = GO_legend_labels) +
  labs(
    y = "Normalized Enrichment Score (NES)",
    color = "Category"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_blank(),
    axis.title =  element_text(size = title_text_size, family = font_family, face = title_face),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, family = font_family, face = title_face, angle = 45, vjust = 1.0, hjust = 1.0),
    axis.text.y = element_text(size = axis_text_size, family = font_family), 
    legend.title = element_text(size = legend_text_size, family = font_family),
    legend.text = element_text(size = legend_text_size, family = font_family), 
    legend.position = "right",
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA), 
  )

# Save the plot
plot_width <- 6
plot_height <- 4
plot_unit <- 'in'
plot_dpi <- 300
file_name <- paste0("/home/joonho345/3_RNA/RNA_Animal/07.Enrichment/GSEA_05.NES_upper/NES_Plot_", Treatment_type, ".png")
ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)


##### without legend 
# Create dot plot for all GO terms with reordered legend
axis_text_size <- 8
legend_text_size <- 8
title_text_size <- 10
title_face <- "bold"
font_family <- "Arial"
plot_target <- ggplot(nes_df, aes(x = Phase, y = NES, color = GO, group = GO)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", size = 0.6) +
  geom_hline(yintercept = c(-2, -1, 1, 2), linetype = "dashed", color = "gray", size = 0.2) +
  geom_vline(xintercept = which(levels(nes_df$Phase) == "MTLEALL") - 0.5, linetype = "dotted", color = "red", size = 0.8) +  
  geom_line(size = 1, alpha = 0.5) +  # Connect dots with lines
  geom_point(size = 3, alpha = 1) +  # Add dots
  scale_color_manual(values = GO_colors) +
  labs(
    y = "Normalized Enrichment Score",
    color = "Category"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_blank(),
    axis.title =  element_text(size = title_text_size, family = font_family, face = title_face),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = axis_text_size, family = font_family, angle = 45, hjust = 1),
    axis.text.y = element_text(size = axis_text_size, family = font_family), 
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA), 
  )

# Save the plot
plot_width <- 5
plot_height <- 3
plot_unit <- 'in'
plot_dpi <- 300
file_name <- paste0("/home/joonho345/3_RNA/RNA_Animal/07.Enrichment/GSEA_05.NES_upper/NES_Plot_", Treatment_type, "_without.png")
ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)

file_name <- paste0("/home/joonho345/3_RNA/RNA_Animal/07.Enrichment/GSEA_05.NES_upper/NES_Plot_", Treatment_type, "_without.svg")
ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit, device = 'svg')


