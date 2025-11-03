library(dplyr)
library(ggplot2)
library(reshape2)

### ONLY HIPPO ###
GSE <- "2019CA"
target <- paste0(GSE, "_scRNA_matrix_all_cluster_1")

GSE <- "GSE186538"
target <- paste0(GSE, "_scRNA_matrix_top1000_cluster_1")


#### Exclude outliers ####
OUTLIER_samples <- c("SRR15406573", "SRR28007131", "SRR28007092", "SRR8669936", "SRR16522693", "SRR28007020")
prop_coldata_df <- coldata_df_hippo[!rownames(coldata_df_hippo) %in% OUTLIER_samples, ]
Diagnosis_vector <- c("MTLEHS", "MTLE", "NL")
prop_coldata_df <- prop_coldata_df %>%
  filter(Diagnosis %in% Diagnosis_vector)

# Fraction from CIBERSORT
Fraction_path <- paste0("/home/joonho345/3_RNA/RNA_Human/08.Deconvolution/03.CIBERSORT_",
                        GSE, "_02Fraction/", target, "/CIBERSORTx_Results.txt")
prop_df <- read.csv(Fraction_path, sep="\t")
prop_df <- prop_df %>%
  filter(Mixture %in% prop_coldata_df$Run)
merged_df <- merge(prop_df, prop_coldata_df, by.x="Mixture", by.y="Run")

melted_df <- melt(merged_df, id.vars=c("Mixture", "Diagnosis"), 
                  measure.vars=c("ExN", "InN", "Oligo", "Astro", "Micro"))
melted_df$Diagnosis <- factor(melted_df$Diagnosis, levels=Diagnosis_vector)
melted_df$variable <- factor(melted_df$variable, 
                             levels = c("ExN", "InN", "Astro", "Oligo", "Micro"))

# Plot
colors_fill <- c("ExN" = "antiquewhite2", 
                 "InN" = "darkseagreen2", 
                 "Astro" = "steelblue", 
                 "Oligo" = "lightgoldenrod3", 
                 "Micro" = "coral")

boxplot <- ggplot(melted_df, aes(x=Diagnosis, y=value, fill=variable)) +
  geom_hline(yintercept = c(0.00, 0.25, 0.50, 0.75, 1.00), linetype = "dashed", color = "gray") +
  geom_jitter(aes(color = variable), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
              size = 0.6, alpha = 0.5) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outlier.shape = NA) +
  scale_fill_manual(values=colors_fill) +
  scale_color_manual(values=colors_fill) +
  
  theme_minimal() + 
  theme(
    axis.text.x = element_text(hjust = 1,family = "Arial"),
    legend.title = element_text(size = 8, family = "Arial"),
    legend.text = element_text(size = 6, family = "Arial"), 
    panel.grid = element_blank(), 
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA), 
    axis.line = element_line(color = "black"),
  ) + 
  
  labs(title="Proportion of Each Cell Type by Diagnosis", 
       x="Diagnosis", y="Proportion")

file_name <- paste("/home/joonho345/3_RNA/RNA_Human/08.Deconvolution/04.CIBER_02FractionPlot_final/",
                   target, "_Boxplot.png", sep="")
ggsave(file_name, plot = boxplot, width = 8, height = 7, dpi = 300)


###################################
### ONLY HIPPO ###
GSE <- "GSE186538"
target <- paste0(GSE, "_scRNA_matrix_top1000_cluster_2")

#### Exclude outliers ####
OUTLIER_samples <- c("SRR15406573", "SRR28007131", "SRR28007092", "SRR8669936", "SRR16522693", "SRR28007020")
prop_coldata_df <- coldata_df_hippo[!rownames(coldata_df_hippo) %in% OUTLIER_samples, ]
Diagnosis_vector <- c("MTLEHS", "MTLE", "NL")
prop_coldata_df <- prop_coldata_df %>%
  filter(Diagnosis %in% Diagnosis_vector)

# Fraction from CIBERSORT
Fraction_path <- paste0("/home/joonho345/3_RNA/RNA_Human/08.Deconvolution/03.CIBERSORT_",
                        GSE, "_02Fraction/", target, "/CIBERSORTx_Results.txt")
prop_df <- read.csv(Fraction_path, sep="\t")
prop_df <- prop_df %>%
  filter(Mixture %in% prop_coldata_df$Run)
merged_df <- merge(prop_df, prop_coldata_df, by.x="Mixture", by.y="Run")

melted_df <- melt(merged_df, id.vars=c("Mixture", "Diagnosis"), 
                  measure.vars=c("ExN_DG", "ExN_CA", "ExN_SUB"))
melted_df$Diagnosis <- factor(melted_df$Diagnosis, levels=Diagnosis_vector)
melted_df$variable <- factor(melted_df$variable, 
                             levels = c("ExN_DG", "ExN_CA", "ExN_SUB"))

# Plot
colors_fill <- c("ExN_DG" = "antiquewhite2", 
                 "ExN_CA" = "darkseagreen2", 
                 "ExN_SUB" = "steelblue")

boxplot <- ggplot(melted_df, aes(x=Diagnosis, y=value, fill=variable)) +
  geom_hline(yintercept = c(0.00, 0.25, 0.50, 0.75, 1.00), linetype = "dashed", color = "gray") +
  geom_jitter(aes(color = variable), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
              size = 0.6, alpha = 0.5) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outlier.shape = NA) +
  scale_fill_manual(values=colors_fill) +
  scale_color_manual(values=colors_fill) +
  
  theme_minimal() + 
  theme(
    axis.text.x = element_text(hjust = 1,family = "Arial"),
    legend.title = element_text(size = 8, family = "Arial"),
    legend.text = element_text(size = 6, family = "Arial"), 
    panel.grid = element_blank(), 
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA), 
    axis.line = element_line(color = "black"),
  ) + 
  
  labs(title="Proportion of Each Cell Type by Diagnosis", 
       x="Diagnosis", y="Proportion")

file_name <- paste("/home/joonho345/3_RNA/RNA_Human/08.Deconvolution/04.CIBER_02FractionPlot_final/",
                   target, "_Boxplot.png", sep="")
ggsave(file_name, plot = boxplot, width = 6, height = 7, dpi = 300)



