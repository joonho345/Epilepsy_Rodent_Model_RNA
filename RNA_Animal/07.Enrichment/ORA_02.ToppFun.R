library(dplyr)
library(ggplot2)

# Define list of data names to process
data_names <- c(
  # '02.' prefix - intersection and unique between M_KAI_IH and R_PPS
  'intersection_M_KAI_IH_IPSI_R_PPS_A_up',
  'unique_M_KAI_IH_IPSI_A_up',
  'unique_R_PPS_A_up',
  'intersection_M_KAI_IH_IPSI_R_PPS_A_down',
  'unique_M_KAI_IH_IPSI_A_down',
  'unique_R_PPS_A_down',
  # '03.' prefix - intersection and unique between M_KAI_IH and MTLE
  'intersection_M_KAI_IH_IPSI_MTLE_up',
  'unique_M_KAI_IH_IPSI_A_MTLE_up',
  'unique_MTLE_M_KAI_IH_IPSI_up',
  'intersection_M_KAI_IH_IPSI_MTLE_down',
  'unique_M_KAI_IH_IPSI_A_MTLE_down',
  'unique_MTLE_M_KAI_IH_IPSI_down',
  # '04.' prefix - intersection and unique between R_PPS and MTLE
  'intersection_R_PPS_MTLE_up',
  'unique_R_PPS_A_MTLE_up',
  'unique_MTLE_R_PPS_up',
  'intersection_R_PPS_MTLE_down',
  'unique_R_PPS_A_MTLE_down',
  'unique_MTLE_R_PPS_down'
)

# Settings
GO_type <- "GO: Biological Process"
GO_name <- "Biological"
HitCountGenome <- 1000

# Loop through each data name
for (DataName in data_names) {
  cat("Processing", DataName, "...\n")
  
  # Determine prefix based on data name (check more specific patterns first)
  if (grepl("^intersection_M_KAI_IH_IPSI_MTLE_|^unique_M_KAI_IH_IPSI_A_MTLE_|^unique_MTLE_M_KAI_IH_IPSI_", DataName)) {
    prefix <- "03."
  } else if (grepl("^intersection_R_PPS_MTLE_|^unique_R_PPS_A_MTLE_|^unique_MTLE_R_PPS_", DataName)) {
    prefix <- "04."
  } else if (grepl("^intersection_M_KAI_IH_IPSI_R_PPS_A_|^unique_M_KAI_IH_IPSI_A_|^unique_R_PPS_A_", DataName)) {
    prefix <- "02."
  } else {
    prefix <- "02."  # Default
  }
  
  # Construct input file path
  InputFile <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_03.ToppFun/", prefix, "ToppFun_", DataName, ".txt")
  
  # Check if input file exists
  if (!file.exists(InputFile)) {
    cat("Warning: Input file", InputFile, "does not exist. Skipping...\n")
    next
  }
  
  # Read and filter data
  data <- read.table(InputFile, header = TRUE, sep = "\t", quote = "", comment.char = "")
  filtered_data <- data %>% filter(Category == GO_type)
  
  # Check if filtered data is empty
  if (nrow(filtered_data) == 0) {
    cat("Warning: No data found for", DataName, "after filtering. Skipping...\n")
    next
  }
  
  # Filter and process data
  top_data <- filtered_data %>% 
    filter(Hit.Count.in.Genome <= HitCountGenome) %>%
    filter(q.value.FDR.B.H < 0.05) %>%
    arrange(q.value.FDR.B.H)
  
  # Check if top_data is empty after filtering
  if (nrow(top_data) == 0) {
    cat("Warning: No significant results found for", DataName, "after filtering. Skipping...\n")
    next
  }
  
  top_data <- top_data %>%
    mutate(log_q_value = -log10(`q.value.FDR.B.H`),
           gene_ratio = `Hit.Count.in.Query.List` / `Hit.Count.in.Genome`)
  
  # Export files
  filename <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_03.ToppFun/", prefix, "ToppFun_", DataName, "_final.txt")
  write.table(top_data, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("Saved", filename, "with", nrow(top_data), "rows\n")
}



