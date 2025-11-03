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
OUTLIER_samples <- c(
  "SRR15406573", 
  "SRR28007131", "SRR28007018", "SRR28007079",
  "SRR16522705", "SRR28007092", "SRR28007081", "SRR28007133", "SRR16522707",
  "SRR16522708", "SRR28007072", "SRR16522706", "SRR8669931", "SRR28007122",
  "SRR28007034", "SRR28007155", "SRR16522692", "SRR28007132", "SRR16522690",
  "SRR8669932", "SRR8669936", "SRR16522693", "SRR16522689", "SRR28007020",
  "SRR16522691", "SRR16522694", "SRR28007089", "SRR8669933"
)
#OUTLIER_samples <- c(
#  "SRR5241783", "SRR9733960", "SRR9733958", "SRR15406573", "SRR8669937",
#  "SRR5241781", "SRR8669940", "SRR28007131", "SRR28007018", "SRR28007079",
#  "SRR16522705", "SRR28007092", "SRR28007081", "SRR28007133", "SRR16522707",
#  "SRR16522708", "SRR28007072", "SRR16522706", "SRR8669931", "SRR28007122",
#  "SRR28007034", "SRR28007155", "SRR16522692", "SRR28007132", "SRR16522690",
#  "SRR8669932", "SRR8669936", "SRR16522693", "SRR16522689", "SRR28007020",
#  "SRR16522691", "SRR16522694", "SRR28007089", "SRR8669933"
#)

merged_pcaData$Outlier <- ifelse(merged_pcaData$Run %in% OUTLIER_samples, "Outlier", "Non-Outlier")

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
file_name <- paste0("/home/joonho345/3_RNA/RNA_Human/05.Sampling/03.DESeq_PCA_Outliers/PCA_plot_",
                    Target_comparison, "_", markers, "_sub_v_outlier.png")
#file_name <- paste0("/home/joonho345/3_RNA/RNA_Human/05.Sampling/03.DESeq_PCA_Outliers/PCA_plot_",
#                    Target_comparison, "_", markers, "_sub_vh_outlier.png")
ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)

