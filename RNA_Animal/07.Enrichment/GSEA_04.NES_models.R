library(ggplot2)
library(dplyr)

# Define GO term
GO_name <- "IONOTROPIC_GLUTAMATE_KAINATE_RECEPTORS"

# File paths
file_paths <- list(
  MTLEALL = "/home/joonho345/3_RNA/RNA_Animal/05.Enrichment/GSEA_03.Output/gsea_report_for_MTLEALL_GOBP_1000_1.tsv",
  M_KAI_ST_BO_O_HA = "/home/joonho345/3_RNA/RNA_Animal/05.Enrichment/GSEA_03.Output/M_KAI_ST_BO_O_HA.tsv",
  M_KAI_ST_BO_O_AC = "/home/joonho345/3_RNA/RNA_Animal/05.Enrichment/GSEA_03.Output/M_KAI_ST_BO_O_AC.tsv",
  M_KAI_ST_BO_O_IM = "/home/joonho345/3_RNA/RNA_Animal/05.Enrichment/GSEA_03.Output/M_KAI_ST_BO_O_IM.tsv",
  M_KAI_ST_BO_O_CR = "/home/joonho345/3_RNA/RNA_Animal/05.Enrichment/GSEA_03.Output/M_KAI_ST_BO_O_CR.tsv",
  M_PILO_Y_HA = "/home/joonho345/3_RNA/RNA_Animal/05.Enrichment/GSEA_03.Output/M_PILO_Y_HA.tsv",
  M_PILO_Y_AC = "/home/joonho345/3_RNA/RNA_Animal/05.Enrichment/GSEA_03.Output/M_PILO_Y_AC.tsv",
  M_PILO_Y_IM = "/home/joonho345/3_RNA/RNA_Animal/05.Enrichment/GSEA_03.Output/M_PILO_Y_IM.tsv",
  R_PPS_O_AC = "/home/joonho345/3_RNA/RNA_Animal/05.Enrichment/GSEA_03.Output/R_PPS_O_AC.tsv",
  R_PPS_O_IM = "/home/joonho345/3_RNA/RNA_Animal/05.Enrichment/GSEA_03.Output/R_PPS_O_IM.tsv",
  R_PPS_O_DOFS = "/home/joonho345/3_RNA/RNA_Animal/05.Enrichment/GSEA_03.Output/R_PPS_O_DOFS.tsv",
  R_PPS_O_CR = "/home/joonho345/3_RNA/RNA_Animal/05.Enrichment/GSEA_03.Output/R_PPS_O_CR.tsv"
)

# Initialize an empty list to store NES values
nes_values <- list()

# Loop through each file and extract NES for the specified GO term
for (name in names(file_paths)) {
  file_path <- file_paths[[name]]
  
  # Read the file
  tsv_data <- read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "\"")
  
  # Filter the data for the specified GO term
  go_data <- tsv_data %>% filter(NAME == GO_name)
  
  # Extract NES if the GO term is found
  if (nrow(go_data) > 0) {
    nes_values[[name]] <- go_data$NES
  } else {
    nes_values[[name]] <- NA  # Set NA if the GO term is not found
  }
}

# Convert the NES values into a DataFrame
nes_df <- data.frame(
  Sample = names(nes_values),
  NES = unlist(nes_values)
)

# Plot NES values
nes_plot <- ggplot(nes_df, aes(x = Sample, y = NES)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = paste("NES for", GO_name), x = "Sample", y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# Save the plot
output_plot <- paste0("/home/joonho345/3_RNA/RNA_Animal/05.Enrichment/GSEA_03.Output/NES_Plot_", GO_name, ".png")
ggsave(output_plot, plot = nes_plot, width = 8, height = 6)

