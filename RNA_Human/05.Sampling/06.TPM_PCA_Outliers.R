library(DESeq2)
library(ggplot2)
library(pasilla)

### all genes
markers <- 'all'

### Choose one matrix to analyze ###
### raw data
dds_data <- adjusted_df_F
dds_coldata <- coldata_df_F
Target_comparison <- 'F_MTLEALL_NL_HIPPO'
contrast_vector <- c("Diagnosis","MTLEALL","NL")


# Check if coldata and count matrix have the same order with respect to samples
if (!all(rownames(dds_coldata) == colnames(dds_data))) {
  stop("Rownames of coldata do not match colnames of count data")
}
## dds
dds <- DESeqDataSetFromMatrix(countData = dds_data,
                              colData = dds_coldata,
                              design = ~ Diagnosis)
## Prefiltering
keep <- rowSums(counts(dds)) >= 400
dds <- dds[keep,]
## Differential expression analysis
dds <- DESeq(dds)
resultsNames(dds)

## pcdDATA
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Diagnosis", "Brain_Location"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
dds_coldata_pca <- dds_coldata[, !colnames(dds_coldata) %in% c("Diagnosis", "Brain_Location")]
merged_pcaData <- cbind(pcaData, dds_coldata_pca)
merged_pcaData$Diagnosis_Sub[merged_pcaData$Diagnosis_Sub == "TLE"] <- "MTLE"
identical_values <- merged_pcaData$name == merged_pcaData$Run
merged_pcaData$Diagnosis_Sub <- factor(merged_pcaData$Diagnosis_Sub, levels = c('MTLEHS', 'MTLE', 'NL'))


########### mark outlier ###########
OUTLIER_samples <- ALL_outliers

## Outliers from TPM expr
ALL_outliers

## Shared outliers from TPM expr & DESeq PCA
shared_outliers <- intersect(OUTLIER_samples, ALL_outliers)
shared_outliers

merged_pcaData$Outlier <- ifelse(merged_pcaData$Run %in% shared_outliers, "Outlier", "Non-Outlier")

## PCA clustering plot (sub-Diagnosis)
axis_text_size <- 8
legend_text_size <- 8
title_text_size <- 12
title_face <- "bold"
font_family <- "Arial"
plot_target <- ggplot(merged_pcaData, aes(x=PC1, y=PC2, color=Diagnosis_Sub)) +
  geom_hline(yintercept = c(-40, -20, 0, 20), linetype = "dashed", color = "gray", size = 0.3) +
  geom_vline(xintercept = c(-20, 0, 20), linetype = "dashed", color = "gray", size = 0.3) +
  geom_vline(xintercept = c(-1), size = 1.2, linetype = "dashed", color = "darkred") +
  #geom_hline(yintercept = c(-24), size = 1.2, linetype = "dashed", color = "darkred") +
  geom_point(size=3) +
  
  geom_point(data = subset(merged_pcaData, Outlier == "Outlier"),
             aes(PC1, PC2),
             shape = 4, # Shape 4 is an "X"
             size = 3,  # Size of the X marker
             color = "darkred", # Color of the X
             stroke = 1.2) + # Thickness of the X
  
  scale_color_manual(
    values = c("MTLEHS" = "darkseagreen2", "MTLE" = "lightsteelblue1", "NL" = "antiquewhite2")
  ) +
  
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  
  theme_minimal() +
  theme(
    axis.title =  element_text(size = title_text_size, family = font_family, face = title_face),
    axis.text.x = element_text(size = axis_text_size, family = font_family),
    axis.text.y = element_text(size = axis_text_size, family = font_family), 
    legend.title = element_text(size = legend_text_size, family = font_family),
    legend.text = element_text(size = legend_text_size, family = font_family), 
    legend.position = "right",
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA)
  )

plot_width <- 7
plot_height <- 6
plot_unit <- 'in'
plot_dpi <- 300
file_name <- paste0("/home/joonho345/3_RNA/RNA_Human/05.Sampling/06.TPM_ComparisonALL/PCA_plot_",
                    Target_comparison, "_", markers, "_TPM_outlier.png")
#file_name <- paste0("/home/joonho345/3_RNA/RNA_Human/05.Sampling/03.DESeq_PCA_Outliers/PCA_plot_",
#                    Target_comparison, "_", markers, "_sub_vh_outlier.png")
ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)




