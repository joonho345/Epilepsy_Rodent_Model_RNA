library(dplyr)
library(ggplot2)
library(reshape2)

### ONLY HIPPO ###
GSE <- "GSE185862"
# ONLY HIPPO: Load coldata_df_M_hippo -> 00.Subgrouping_H.R
prop_coldata_df <- coldata_df_M_hippo
Treatment_Lateral_vector <- as.vector(unique(prop_coldata_df$Treatment_Lateral))

target <- paste0(GSE, "_scRNA_matrix_all_cluster_1")
target <- paste0(GSE, "_scRNA_matrix_top2000_cluster_1")
target <- paste0(GSE, "_scRNA_matrix_top1000_cluster_1")


###################################
# Fraction from CIBERSORT
Fraction_path <- paste0("/home/joonho345/3_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_",
                        GSE, "_02Fraction/", target, "/CIBERSORTx_Results.txt")
prop_df <- read.csv(Fraction_path, sep="\t")
merged_df <- merge(prop_df, prop_coldata_df, by.x="Mixture", by.y="Run")

# Melt the dataframe for easier plotting
melted_df <- melt(merged_df, id.vars=c("Mixture", "Treatment_Lateral"), 
                  measure.vars=c("ExN", "InN", "Astro", "Micro"))

# Set the order of the x-axis categories
melted_df$Treatment_Lateral <- factor(melted_df$Treatment_Lateral, levels=Treatment_Lateral_vector)

###################################
# Plot the proportions of each cell type by diagnosis group
boxplot <- ggplot(melted_df, aes(x=Treatment_Lateral, y=value, fill=variable)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales="free") +
  labs(title="Proportion of Each Cell Type by Treatment_Lateral", 
       x="Treatment_Lateral", y="Proportion") +
  theme(axis.text.x = element_text(angle=45, hjust=1))

file_name <- paste("/home/joonho345/3_RNA/RNA_Animal/06.Deconvolution/04.CIBER_02FractionPlot_cluster_1/",
                   target, "_Boxplot.png", sep="")
ggsave(file_name, plot = boxplot, width = 20, height = 15, dpi = 300)


# Plot the proportions of each cell type by diagnosis group
stacked_barplot <- ggplot(melted_df, aes(x=Treatment_Lateral, y=value, fill=variable)) +
  geom_bar(stat="identity", position="fill") +
  labs(title="Proportion of Each Cell Type by Treatment_Lateral", 
       x="Treatment_Lateral", y="Proportion") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_y_continuous(labels = scales::percent) +
  guides(fill=guide_legend(title="Cell Type"))

file_name <- paste("/home/joonho345/3_RNA/RNA_Animal/06.Deconvolution/04.CIBER_02FractionPlot_cluster_1/",
                   target, "_Stackplot.png", sep="")
ggsave(file_name, plot = stacked_barplot, width = 10, height = 10, dpi = 300)

