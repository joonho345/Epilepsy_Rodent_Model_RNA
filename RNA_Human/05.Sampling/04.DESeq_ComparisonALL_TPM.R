library(ggplot2)
library(ggpubr)
library(scales)

###############################################################
##################### Use TPM Count matrix ####################
DATA <- 'TPM'
###############################################################

# vline PC1 -1 & hline PC2 -24
OUTLIER_samples <- c(
  "SRR5241783", "SRR9733960", "SRR9733958", "SRR15406573", "SRR8669937",
  "SRR5241781", "SRR8669940", "SRR28007131", "SRR28007018", "SRR28007079",
  "SRR16522705", "SRR28007092", "SRR28007081", "SRR28007133", "SRR16522707",
  "SRR16522708", "SRR28007072", "SRR16522706", "SRR8669931", "SRR28007122",
  "SRR28007034", "SRR28007155", "SRR16522692", "SRR28007132", "SRR16522690",
  "SRR8669932", "SRR8669936", "SRR16522693", "SRR16522689", "SRR28007020",
  "SRR16522691", "SRR16522694", "SRR28007089", "SRR8669933"
)
subset <- 'vh'

# vline PC1 -1
OUTLIER_samples <- c(
  "SRR15406573", 
  "SRR28007131", "SRR28007018", "SRR28007079",
  "SRR16522705", "SRR28007092", "SRR28007081", "SRR28007133", "SRR16522707",
  "SRR16522708", "SRR28007072", "SRR16522706", "SRR8669931", "SRR28007122",
  "SRR28007034", "SRR28007155", "SRR16522692", "SRR28007132", "SRR16522690",
  "SRR8669932", "SRR8669936", "SRR16522693", "SRR16522689", "SRR28007020",
  "SRR16522691", "SRR16522694", "SRR28007089", "SRR8669933"
)
subset <- 'v'

# select gene
Target_gene <- "GFAP"
y_intercept <- c(0, 2500, 5000, 7500, 10000)
y_coord <- c(0,10000)

Target_gene <- "NOTCH1"
y_intercept <- c(0, 25, 50, 75, 100)
y_coord <- c(0,100)


# choose target groups (hippo)
Target_matrix <- adjusted_df_hippo
Target_coldata <- coldata_df_hippo
specified_groups <- c("MTLEHS", "MTLE", "NL")
Target_groups <- 'HIPPO'
group_colors <- c("MTLEHS" = "darkseagreen2", "MTLE" = "lightsteelblue1", "NL" = "antiquewhite2")

# Ensure the sample names match between the expression data and coldata
Target_coldata <- Target_coldata[match(colnames(Target_matrix), rownames(Target_coldata)), ]

# Extract the expression data for the target gene
gene_expression <- Target_matrix[Target_gene, ]
gene_expression <- t(gene_expression)
gene_expression <- as.data.frame(gene_expression)
gene_expression$sample_name <- rownames(gene_expression)
colnames(gene_expression)[1] <- "target_gene_expression"

gene_expression_merged <- merge(gene_expression, Target_coldata, by.x = "sample_name", by.y = "Run")
gene_expression_filtered <- gene_expression_merged[gene_expression_merged$Diagnosis %in% specified_groups, ]
gene_expression_filtered$Diagnosis <- factor(gene_expression_filtered$Diagnosis, levels = specified_groups)
gene_expression_filtered$highlight <- ifelse(gene_expression_filtered$sample_name %in% OUTLIER_samples, "highlight", "normal")

# Prepare data for plotting
expression_data <- data.frame(
  Expression = gene_expression_filtered$target_gene_expression,
  Group = gene_expression_filtered$Diagnosis,
  Highlight = gene_expression_filtered$highlight
)


############## PLOT ###############
axis_text_size <- 8
legend_text_size <- 8
title_text_size <- 12
title_face <- "bold"
font_family <- "Arial"
plot_target <- ggplot(expression_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_hline(yintercept = y_intercept, linetype = "dashed", color = "gray", size = 0.3) +
  geom_boxplot(width = 0.7, outlier.shape = NA, position = position_dodge(width = 0.3)) +
  
  geom_jitter(
    data = subset(expression_data, Highlight == "normal"),
    aes(color = Highlight),
    width = 0.2, size = 0.8, alpha = 0.4
  ) +
  geom_jitter(
    data = subset(expression_data, Highlight == "highlight"),
    aes(color = Highlight),
    width = 0.2, size = 1.0, alpha = 1.0, shape = 4, stroke = 1.0
  ) +
  
  scale_fill_manual(values = group_colors) + 
  scale_color_manual(values = c("highlight" = "darkred", "normal" = "gray")) + 
  labs(title = paste0(DATA, " Expression of ", Target_gene),
       y = paste0(Target_gene," Normalized TPM Expression")) +
  
  coord_cartesian(ylim = y_coord) + 
  theme_minimal() + 
  theme(
    axis.title =  element_text(size = title_text_size, family = font_family, face = title_face),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = axis_text_size, family = font_family),
    legend.title = element_text(size = legend_text_size, family = font_family),
    legend.text = element_text(size = legend_text_size, family = font_family), 
    legend.position = "right",
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA), 
  ) 

plot_width <- 3
plot_height <- 5
plot_unit <- 'in'
plot_dpi <- 300
file_name <- paste("/home/joonho345/3_RNA/RNA_Human/05.Sampling/04.DESeq_ComparisonALL/Boxplot_", Target_gene, "_",
                   Target_groups, "_", subset, "_", DATA, ".png", sep="")
ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)



