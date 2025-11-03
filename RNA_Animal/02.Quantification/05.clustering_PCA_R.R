### Choose one matrix to analyze ###
# Rat
prcomp_data <- merged_matrix_R
prcomp_name <- "merged_matrix_R"
prcomp_data <- adjusted_merged_matrix_R
prcomp_name <- "adjusted_merged_matrix_R"
prcomp_data <- adjusted_merged_matrix_R_1
prcomp_name <- "adjusted_merged_matrix_R_1"

prcomp_coldata <- coldata_df_R


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
cohort_colors_R <- c(
  "GSE50079"  = "#CAB2D6",  # Lavender
  "GSE75120"  = "#FDBF6F",  # Light Orange
  "GSE75402"  = "#FF6F61",  # Coral
  "GSE136913" = "#F0E442",  # Yellow
  "GSE137473" = "#6A3D9A",  # Deep Purple
  "GSE143555" = "#377EB8"   # Dark Blue
)
treatment_shapes_R <- c(
  "PILO"             = 20,  # Circle
  "CTL_SAL_IP"       = 17,  # Triangle
  "CTL_TBI_IPSI"     = 15,  # Square
  "TBI_IPSI"         = 18,  # Diamond
  "CTL_AMG_STI_IPSI" = 19,  # Solid Circle
  "AMG_STI_IPSI"     = 8,   # Star
  "KAI_SUB"          = 4,   # Cross
  "CTL_SAL_SUB"      = 3,   # Plus
  "PPS"              = 7,   # Empty Triangle
  "CTL_PPS"          = 9,   # Diamond with Plus
  "KAI_IP"           = 10   # Empty Square
)
sequencing_shapes_R <- c(
  "GenomeAnalyzerIIx" = 16,  # Circle
  "HiSeq2500"         = 17,  # Triangle
  "HiSeq4000"         = 15   # Square
)

### ggplot ###
######## Rat #########
# cohort & treatment
name_shape <- 'treatment'
Public_pca_P <- ggplot(data = Public_pca_df,
                       aes(x = PC1,  y = PC2)) +
  geom_point(aes(color = GEO, shape = Treatment_Lateral), size = 4) +
  scale_colour_manual(values = cohort_colors_R) +
  scale_shape_manual(values = treatment_shapes_R) +
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
  scale_colour_manual(values = cohort_colors_R) +
  scale_shape_manual(values = sequencing_shapes_R) +
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

