library(ggplot2)

file_path <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/08.Enrichment/GSEA_03.Output/gsea_report_for_MTLEALL_GOBP_1000_new_revision.tsv"
file_name <- "_1000_new_revision"


df <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
# Create a flag for significant terms (p < 0.05) before modifying p-values
df$is_significant <- df$NOM.p.val < 0.05
# Handle p-value = 0 by replacing with machine minimum to avoid log(0) = Inf
df$NOM.p.val[df$NOM.p.val == 0] <- .Machine$double.xmin
df$log_pval <- -log10(df$NOM.p.val)
# Cap log_pval at a maximum for visualization (e.g., 5 or 10)
# This ensures extremely significant terms (p < 1e-5 or p < 1e-10) are still visible
max_log_pval <- 5  # Adjust this value as needed (5 = p < 1e-5, 10 = p < 1e-10)
df$log_pval_capped <- pmin(df$log_pval, max_log_pval)

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

# Create a color map for axis.text.y
color_map <- setNames(df$ID_color, df$GO_NAME)
# Set hline positions after 9th and 22nd items
hline_positions <- c(17, 30)
threshold <- -log10(0.05)

# Forest plot
axis_text_size <- 8
legend_text_size <- 8
title_text_size <- 12
title_face <- "bold"
font_family <- "Arial"
plot_target <- ggplot(df, aes(x = NES, y = GO_NAME, color = log_pval_capped)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "darkgray") + 
  geom_hline(yintercept = hline_positions[1] - 0.5, color = "gray40", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = hline_positions[2] - 0.5, color = "gray40", linetype = "dashed", linewidth = 0.5) +
  
  geom_segment(aes(x = 0, xend = NES, y = GO_NAME, yend = GO_NAME), color = "gray60", linewidth = 1.0) +
  geom_point(size = 5.0) + 
  # Add gray border for all significant terms (p < 0.05)
  geom_point(data = subset(df, is_significant), 
             shape = 21, size = 5.0, stroke = 1.5, color = "gray40", fill = NA) +

  scale_x_continuous(limits = c(-2, 2)) + 
  scale_color_gradient(low = "gray95", high = "firebrick", limits = c(0, max_log_pval), 
                       name = "-log10(p-value)",
                       guide = guide_colorbar(title.position = "top")) +
  
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway Name"
  ) +
  
  theme_minimal() + 
  theme(
    plot.title = element_blank(),
    axis.title =  element_blank(),
    axis.text.x = element_text(size = title_text_size, family = font_family, face = title_face),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10, family = font_family, color = 'black'),
    legend.title = element_text(size = legend_text_size, family = font_family),
    legend.text = element_text(size = legend_text_size, family = font_family), 
    legend.position = "right",
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA), 
    panel.grid.major.y = element_line(color = "gray90"),
  )

# Save the plot
plot_width <- 9.0
plot_height <- 8.7
plot_unit <- 'in'
plot_dpi <- 600
filename <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/08.Enrichment/GSEA_03.Output/GOBP_MTLEALL_NL", file_name, ".png")
ggsave(filename, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)
