### Choose one matrix to analyze ###
merged_matrix <- read.table("/home/joonho345/3_RNA/RNA_Human/02.Quantification/merged_matrix.txt",
                          sep = "\t", header = TRUE, row.names = 1)
prcomp_data <- merged_matrix
prcomp_name <- "merged_matrix"

adjusted_merged_matrix_1 <- read.table("/home/joonho345/3_RNA/RNA_Human/02.Quantification/adjusted_merged_matrix_1.txt",
                            sep = "\t", header = TRUE, row.names = 1)
prcomp_data <- adjusted_merged_matrix_1
prcomp_name <- "adjusted_merged_matrix_1"

##
merged_matrix_TPM <- read.table("/home/joonho345/3_RNA/RNA_Human/04.Normalization/merged_matrix_TPM.txt",
                            sep = "\t", header = TRUE, row.names = 1)
prcomp_data <- merged_matrix_TPM
prcomp_name <- "merged_matrix_TPM"

adjusted_merged_matrix_1_TPM <- read.table("/home/joonho345/3_RNA/RNA_Human/04.Normalization/adjusted_merged_matrix_1_TPM.txt",
                                       sep = "\t", header = TRUE, row.names = 1)
prcomp_data <- adjusted_merged_matrix_1_TPM
prcomp_name <- "adjusted_merged_matrix_1_TPM"


#prcomp_data <- adjusted_merged_matrix
#prcomp_name <- "adjusted_merged_matrix"
#prcomp_data <- adjusted_merged_matrix_2
#prcomp_name <- "adjusted_merged_matrix_2"

coldata_df <- read.table("/home/joonho345/3_RNA/RNA_Human/02.Quantification/filtered_coldata.txt",
                         sep = "\t", header = TRUE, row.names = 1, fill = TRUE)
prcomp_coldata <- coldata_df


### PCA analyis using prcomp ###
Public_pca <- prcomp(
  prcomp_data,
  center = T,
  scale. = T
)

summary_df <- summary(Public_pca)$importance
Public_pca_df <- as.data.frame(Public_pca[2]$rotation)

Public_pca_df[,"PRJNA"] <- prcomp_coldata$PRJNA
Public_pca_df[,"Diagnosis"] <- prcomp_coldata$Diagnosis
Public_pca_df[,"Sequencing"] <- prcomp_coldata$Sequencing

Public_pca_per <- round(summary_df["Proportion of Variance",] * 100, 1)
Public_pca_per <- paste0(colnames(Public_pca_df)[which(substr(colnames(Public_pca_df),1,2) == "PC")]," (",paste0(as.character(Public_pca_per),"%",")", sep=""))

#choose color aes
cohort_colors <- c(
  "PRJNA1073977" = "#FDE3A7", # Lighter orange
  "PRJNA1077986" = "#FAB78B", # Lighter red-orange
  "PRJNA280563" = "#F58A72",  # Lighter deep red
  "PRJNA290212" = "#C48EDC",  # Lighter purple
  "PRJNA322318" = "#C8749F",  # Lighter magenta
  "PRJNA373909" = "#8CA6E5",  # Lighter blue
  "PRJNA525671" = "#85C4C9",  # Lighter teal
  "PRJNA556159" = "#FFC6D6",  # Lighter pink
  "PRJNA600414" = "#9A74B8",  # Lighter dark purple
  "PRJNA753504" = "#93DCB1",  # Lighter light green
  "PRJNA773413" = "#D8F296",  # Lighter lime green
  "PRJNA787219" = "#F5ED92"   # Lighter yellow
)
# exclude PRJNA787219
Public_pca_df <- subset(Public_pca_df, PRJNA != "PRJNA787219")

### ggplot ###
axis_text_size <- 8
legend_text_size <- 8
title_text_size <- 12
title_face <- "bold"
font_family <- "Arial"
plot_target <- ggplot(data = Public_pca_df, aes(x = PC1,  y = PC2)) +
  geom_point(aes(color = PRJNA), size = 2.5, alpha = 0.8) +
  scale_colour_manual(values = cohort_colors) +
  
  labs(x = Public_pca_per[1],
       y = Public_pca_per[2]) +

  theme_minimal() +
  theme(
    plot.title = element_blank(),
    axis.title =  element_text(size = title_text_size, family = font_family, face = title_face),
    axis.text.x = element_text(size = axis_text_size, family = font_family),
    axis.text.y = element_text(size = axis_text_size, family = font_family), 
    legend.title = element_text(size = legend_text_size, family = font_family),
    legend.text = element_text(size = legend_text_size, family = font_family), 
    legend.position = "right",
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA), 
    plot.background = element_rect(fill = "transparent", color = NA)
  )

plot_width <- 7
plot_height <- 4
plot_unit <- 'in'
plot_dpi <- 300
file_name <- paste0("/home/joonho345/3_RNA/RNA_Human/02.Quantification/05.clustering_PCA/02.PCA_plot_",prcomp_name,".png")
ggsave(file_name, plot = plot_target, width = plot_width, height = plot_height, dpi = plot_dpi, units = plot_unit)
