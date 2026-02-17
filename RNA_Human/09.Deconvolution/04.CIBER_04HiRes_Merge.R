library(dplyr)

cohort <- 'GSE160189_186538_Integration'
target <- paste0(cohort, "_scRNA_matrix_DEGs_CellType1")

celltype <- "ExN_H"
celltype <- "InN_H"
celltype <- "Astro_H"
celltype <- "Micro_H"
celltype <- "Oligo_H"
celltype <- "OPC_H"
celltype <- "Endo_H"

window <- 'Window28'

# Set the working directory
setwd(paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/03.CIBERSORT_", cohort, "_04HiRes/", target))
subdirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
data_list <- list()

# Loop through each subdirectory and read the data files
for (dir in subdirs) {
  file_path <- file.path(dir, paste0("CIBERSORTxHiRes_NA_", celltype, "_", window, ".txt"))
  if (file.exists(file_path)) {
    data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    rownames(data) <- NULL
    data_list[[dir]] <- data
  }
}

# Bind all data frames by rows
merged_data <- do.call(rbind, data_list)
merged_data <- merged_data[complete.cases(merged_data), ]

#### Group by GeneSymbol and calculate the mean for duplicated genes ####
# There are some duplicated genes across different gene sets #
merged_data <- merged_data %>%
  group_by(GeneSymbol) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  ungroup()

# unique genes -> 2295
# unique genes without NaN -> 2150

# Set GeneSymbol as rownames and remove the first column
merged_data <- as.data.frame(merged_data)
rownames(merged_data) <- merged_data$GeneSymbol
merged_data <- merged_data[, -1]

# Save the merged data to a file
write.table(merged_data, file = paste0("Merged_", celltype, ".txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

