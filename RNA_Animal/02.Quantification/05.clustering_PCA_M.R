### Choose one matrix to analyze ###
# Mouse
prcomp_data <- merged_matrix_M
prcomp_name <- "merged_matrix_M"
prcomp_data <- adjusted_merged_matrix_M
prcomp_name <- "adjusted_merged_matrix_M"
prcomp_data <- adjusted_merged_matrix_M_1
prcomp_name <- "adjusted_merged_matrix_M_1"

prcomp_coldata <- coldata_df_M


### PCA analyis using prcomp ###
Public_pca <- prcomp(
  prcomp_data,
  center = T,
  scale. = T
)

summary_df <- summary(Public_pca)$importance
Public_pca_df <- as.data.frame(Public_pca[2]$rotation)

Public_pca_df[,"GEO"] <- prcomp_coldata$GEO
Public_pca_df[,"Treatment_Lateral"] <- prcomp_coldata$Treatment_Lateral
Public_pca_df[,"Sequencing"] <- prcomp_coldata$Sequencing

Public_pca_per <- round(summary_df["Proportion of Variance",] * 100, 1)
Public_pca_per <- paste0(colnames(Public_pca_df)[which(substr(colnames(Public_pca_df),1,2) == "PC")]," (",paste0(as.character(Public_pca_per),"%",")", sep=""))

# aes
cohort_colors_M <- c(
  "GSE72402"  = "#1F78B4",  # Blue
  "GSE99577"  = "#33A02C",  # Green
  # "GSE138100" = "#E31A1C",  # Red
  "GSE148028" = "#FF7F00",  # Orange
  # "GSE153976" = "#6A3D9A",  # Purple
  "GSE198498" = "#B15928",  # Brown
  "GSE205373" = "#A6CEE3",  # Light Blue
  # "GSE205454" = "#B2DF8A",  # Light Green
  "GSE213393" = "#FB9A99",  # Light Pink
  "GSE241219" = "#FDBF6F"   # Light Orange
)
treatment_shapes_M <- c(
  "CTL_SAL_IP"           = 20,  # Circle
  "PILO"                 = 17,  # Triangle
  "CTL_SAL_HIPPO_BOTH"   = 15,  # Square
  "KAI_HIPPO_BOTH"       = 18,  # Diamond
  "KAI_HIPPO_IPSI"       = 19,  # Solid Circle
  "CTL_SAL_HIPPO_IPSI"   = 8,   # Star
  "KAI_IP"               = 4,   # Cross
  "KAI_HIPPO_CONTRA"     = 3,   # Plus
  "CTL_SAL_HIPPO_CONTRA" = 7    # Empty Triangle
)
sequencing_shapes_M <- c(
  "GenomeAnalyzerII" = 20,  # Circle
  "NextSeq500"       = 17,  # Triangle
  # "HiSeq3000"        = 15,  # Square
  "HiSeq2500"        = 18,  # Diamond
  "NovaSeq6000"      = 15,  # Square
  "HiSeq1500"        = 8    # Star
)

### ggplot ###
######## Mouse #########
# cohort & treatment
name_shape <- 'treatment'
Public_pca_P <- ggplot(data = Public_pca_df,
                       aes(x = PC1,  y = PC2)) +
  geom_point(aes(color = GEO, shape = Treatment_Lateral), size = 4) +
  scale_colour_manual(values = cohort_colors_M) +
  scale_shape_manual(values = treatment_shapes_M) +
  theme_bw() +
  labs(title = "Expression pattern of dataset") +
  #eliminates background, gridlines, and chart border
  theme(axis.title       = element_text(size = 20,
                                        face = 'plain'),
        title            = element_text(size = 20,
                                        face = 'bold'),
        plot.title       = element_text(hjust = 0.5,
                                        margin = ggplot2::margin(5,20,20,20)),
        plot.background  = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(colour = "black",
                                        fill   = NA,
                                        linewidth   = 1,),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.line = element_line(color = 'black'),
        axis.text.x = element_text(size = 15,
                                   colour = "black",
                                   vjust = -0.5),
        axis.text.y = element_text(size = 15,
                                   colour = "black",
                                   hjust = -0.1),
        axis.title.x = element_text(size = 15,
                                    face = "plain",
                                    vjust = -1,
                                    margin = ggplot2::margin(15,0,20,0)),
        axis.title.y = element_text(size = 15,
                                    face = "plain",
                                    vjust= 4,
                                    margin = ggplot2::margin(0,15,0,20))) +
  # axis, main title
  labs(#title = "Top1000 before PCA",
    x     = Public_pca_per[1],
    y     = Public_pca_per[2])
filepath <- sprintf("/home/joonho345/3_RNA/RNA_Animal/02.Quantification/05.clustering_PCA/02.PCA_plot_%s_%s.png", 
                    prcomp_name, name_shape)
ggsave(filepath, plot = Public_pca_P, width = 14, height = 8, units = "in")

# cohort & sequencing
name_shape <- 'sequencing'
Public_pca_P <- ggplot(data = Public_pca_df,
                       aes(x = PC1,  y = PC2)) +
  geom_point(aes(color = GEO, shape = Sequencing), size = 4) +
  scale_colour_manual(values = cohort_colors_M) +
  scale_shape_manual(values = sequencing_shapes_M) +
  theme_bw() +
  labs(title = "Expression pattern of dataset") +
  #eliminates background, gridlines, and chart border
  theme(axis.title       = element_text(size = 20,
                                        face = 'plain'),
        title            = element_text(size = 20,
                                        face = 'bold'),
        plot.title       = element_text(hjust = 0.5,
                                        margin = ggplot2::margin(5,20,20,20)),
        plot.background  = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(colour = "black",
                                        fill   = NA,
                                        linewidth   = 1,),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.line = element_line(color = 'black'),
        axis.text.x = element_text(size = 15,
                                   colour = "black",
                                   vjust = -0.5),
        axis.text.y = element_text(size = 15,
                                   colour = "black",
                                   hjust = -0.1),
        axis.title.x = element_text(size = 15,
                                    face = "plain",
                                    vjust = -1,
                                    margin = ggplot2::margin(15,0,20,0)),
        axis.title.y = element_text(size = 15,
                                    face = "plain",
                                    vjust= 4,
                                    margin = ggplot2::margin(0,15,0,20))) +
  # axis, main title
  labs(#title = "Top1000 before PCA",
    x     = Public_pca_per[1],
    y     = Public_pca_per[2])
filepath <- sprintf("/home/joonho345/3_RNA/RNA_Animal/02.Quantification/05.clustering_PCA/02.PCA_plot_%s_%s.png", 
                    prcomp_name, name_shape)
ggsave(filepath, plot = Public_pca_P, width = 14, height = 8, units = "in")
