library(dplyr)
library(pheatmap)
library(Cairo)
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(tidyr)

########### Correlation coefficient ##########
# MOUSE
input_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/04.CellType_Matrix/01.Correlation_Epileptogenesis_Allcelltype_1.txt")
correlation_matrix_CELLTYPE <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
df_list <- list(correlation_matrix_CELLTYPE) 
combined_df_M <- do.call(rbind, df_list)

rownames(combined_df_M) <- c("Excitatory Neuron", "Inhibitory Neuron", "Astrocyte", "Microglia", "Oligodendrocyte", "OPC", "Endothelial Cell")
colnames(combined_df_M) <- gsub("\\.", "_", colnames(combined_df_M))


# RAT
input_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/05.GO_Matrix/01.Correlation_WholeGO_upper_Rat_FILTERED.txt")
correlation_matrix_GO <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
df_list <- list(correlation_matrix_GO)
combined_df_R <- do.call(rbind, df_list)

rownames(combined_df_R) <- c("Neurotransmission", "Ion Channel", "Neuroinflammation", "Neuronal Death", "TNF Signaling",
                             "Neurogenesis", "Gliogenesis", "Mossy Fiber Sprouting", "Synaptic Plasticity")
colnames(combined_df_R) <- gsub("\\.", "_", colnames(combined_df_R))

additional_rows <- c("Excitatory Neuron", "Inhibitory Neuron", "Astrocyte", "Microglia", "Oligodendrocyte", "OPC", "Endothelial Cell")
empty_rows <- data.frame(matrix(NA, nrow = length(additional_rows), ncol = ncol(combined_df_R)))
rownames(empty_rows) <- additional_rows
colnames(empty_rows) <- colnames(combined_df_R)  
combined_df_R <- rbind(empty_rows)


# combine
correlation_matrix_W <- cbind(combined_df_M, combined_df_R)
Group_order <- c(
  "M_KAI_IH_IPSI_A_HA", "M_KAI_IH_IPSI_A_AC", "M_KAI_IH_IPSI_A_IM", "M_KAI_IH_IPSI_A_CR",
  "M_KAI_IH_CON_A_AC", "M_KAI_IH_CON_A_IM", "M_KAI_IH_CON_A_CR", 
  "M_KAI_IA_IPSI_A_AC", "M_KAI_IA_IPSI_A_IM", "M_KAI_IA_IPSI_A_CR",
  "M_KAI_IP_A_HA", "R_KAI_IP_A_CR", "R_KAI_SUB_I_IM", 
  "M_PILO_IP_A_HA", "M_PILO_IP_A_AC", "M_PILO_IP_A_IM", "M_PILO_IP_A_CR", 
  "R_PILO_IP_A_CR", "R_PPS_A_AC", "R_PPS_A_DOFS", "R_PPS_A_IM", "R_PPS_A_CR", 
  "R_AMG_IPSI_A_CR", "R_TBI_IPSI_A_CR"
)
Group_order <- Group_order[Group_order %in% colnames(correlation_matrix_W)]
correlation_matrix_W <- correlation_matrix_W[, Group_order, drop = FALSE]


################# p value #################
# MOUSE
#input_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human_old/10.Correlation/00.Correlation_Matrix_Mouse_FILTERED.txt")
#correlation_matrix_SAMPLE <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
#input_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human_old/12.Correlation_Matrix/04.Pvalue_WholeGO_upper_Mouse_FILTERED.txt")
#correlation_matrix_GO <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
input_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/04.CellType_Matrix/02.Pvalue_Epileptogenesis_Allcelltype_1.txt")
correlation_matrix_CELLTYPE <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
df_list <- list(correlation_matrix_CELLTYPE) 
combined_df_M <- do.call(rbind, df_list)

rownames(combined_df_M) <- c("Excitatory Neuron", "Inhibitory Neuron", "Astrocyte", "Microglia", "Oligodendrocyte", "OPC", "Endothelial Cell")
colnames(combined_df_M) <- gsub("\\.", "_", colnames(combined_df_M))

# RAT
input_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/05.GO_Matrix/02.Pvalue_WholeGO_upper_Rat_FILTERED.txt")
correlation_matrix_GO <- read.table(file = input_file, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
df_list <- list(correlation_matrix_GO)
combined_df_R <- do.call(rbind, df_list)

rownames(combined_df_R) <- c("Neurotransmission", "Ion Channel", "Neuroinflammation", "Neuronal Death", "TNF Signaling",
                             "Neurogenesis", "Gliogenesis", "Mossy Fiber Sprouting", "Synaptic Plasticity")
colnames(combined_df_R) <- gsub("\\.", "_", colnames(combined_df_R))

additional_rows <- c("Excitatory Neuron", "Inhibitory Neuron", "Astrocyte", "Microglia", "Oligodendrocyte", "OPC", "Endothelial Cell")
empty_rows <- data.frame(matrix(NA, nrow = length(additional_rows), ncol = ncol(combined_df_R)))
rownames(empty_rows) <- additional_rows
colnames(empty_rows) <- colnames(combined_df_R)  
combined_df_R <- rbind(empty_rows)

# combine
pvalue_matrix_W <- cbind(combined_df_M, combined_df_R)
Group_order <- c(
  "M_KAI_IH_IPSI_A_HA", "M_KAI_IH_IPSI_A_AC", "M_KAI_IH_IPSI_A_IM", "M_KAI_IH_IPSI_A_CR",
  "M_KAI_IH_CON_A_AC", "M_KAI_IH_CON_A_IM", "M_KAI_IH_CON_A_CR", 
  "M_KAI_IA_IPSI_A_AC", "M_KAI_IA_IPSI_A_IM", "M_KAI_IA_IPSI_A_CR",
  "M_KAI_IP_A_HA", "R_KAI_IP_A_CR", "R_KAI_SUB_I_IM", 
  "M_PILO_IP_A_HA", "M_PILO_IP_A_AC", "M_PILO_IP_A_IM", "M_PILO_IP_A_CR", 
  "R_PILO_IP_A_CR", "R_PPS_A_AC", "R_PPS_A_DOFS", "R_PPS_A_IM", "R_PPS_A_CR", 
  "R_AMG_IPSI_A_CR", "R_TBI_IPSI_A_CR"
)
Group_order <- Group_order[Group_order %in% colnames(pvalue_matrix_W)]
pvalue_matrix_W <- pvalue_matrix_W[, Group_order, drop = FALSE]

################ total ############
correlation_matrix_W
pvalue_matrix_W



# Ensure both matrices are properly formatted
correlation_matrix_W <- as.matrix(correlation_matrix_W)
pvalue_matrix_W <- as.matrix(pvalue_matrix_W)
cor_long <- as.data.frame(as.table(correlation_matrix_W))
pval_long <- as.data.frame(as.table(pvalue_matrix_W))
colnames(cor_long) <- c("Process", "Animal_Model", "Correlation")
colnames(pval_long) <- c("Process", "Animal_Model", "PValue")
plot_df <- left_join(cor_long, pval_long, by = c("Process", "Animal_Model"))

plot_df <- plot_df %>%
  mutate(LogP = -log10(PValue))
process_order <- c("Excitatory Neuron", "Inhibitory Neuron", "Astrocyte", "Microglia", "Oligodendrocyte", "OPC", "Endothelial Cell")
plot_df$Process <- factor(plot_df$Process, levels = rev(process_order)) 

max_logp <- max(plot_df$LogP, na.rm = TRUE) 
max_correlation <- max(plot_df$Correlation, na.rm = TRUE) 
min_logp <- min(plot_df$LogP, na.rm = TRUE) 
min_correlation <- min(plot_df$Correlation, na.rm = TRUE)  #


# Define dot plot
axis_text_size <- 8
legend_text_size <- 8
title_text_size <- 12
font_family <- "Arial"
plot_target <- ggplot(plot_df, aes(x = Animal_Model, y = Process, size = Correlation, color = LogP)) +
  geom_hline(aes(yintercept = as.numeric(Process)), linetype = "dotted", color = "gray70") +
  geom_vline(aes(xintercept = as.numeric(Animal_Model)), linetype = "dotted", color = "gray70") +
  geom_point() +
  scale_color_gradient(low = "gray95", high = "firebrick", limits = c(0, max_logp)) +  # Dynamic max
  scale_size_continuous(range = c(1, 10), limits = c(min_correlation, max_correlation)) +  # Dynamic max
  theme_minimal() +
  theme(
    plot.title = element_blank(),
    axis.title =  element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = legend_text_size, family = font_family),
    legend.text = element_text(size = legend_text_size, family = font_family), 
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 48, r = 20, b = 47, l = 20, unit = "pt")
  ) +
  labs(x = "Animal Model", y = "Biological Process", color = "-log10(PValue)", size = "Correlation")
plot_width <- 15
plot_height <- 3.5
plot_unit <- 'in'
plot_dpi <- 300
file_name <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/10.Correlation/04.CellType_Matrix/03.Dotplot_Celltype_Epileptogenesis.png")
ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)
