# Install and load dplyr if not already installed
# install.packages("dplyr")
library(dplyr)

GSE <- 'GSE185862'
target <- paste0(GSE, "_10X_scRNA_matrix_DEGs_CellType1")

celltype <- "ExN-M"
celltype <- "InN-M"
celltype <- "Astro-M"
celltype <- "Micro-M"
celltype <- "Oligo-M"
celltype <- "OPC-M"
celltype <- "Endo-M"

window <- 'Window28'

# Set the working directory
setwd(paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_", GSE, "_04HiRes/", target))
subdirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
data_list <- list()

# Convert celltype from hyphen format (e.g., "Endo-M") to underscore format (e.g., "Endo_M") for filename
celltype_filename <- gsub("-", "_", celltype)

# Loop through each subdirectory and read the data files
for (dir in subdirs) {
  file_path <- file.path(dir, paste0("CIBERSORTxHiRes_NA_", celltype_filename, "_", window, ".txt"))
  if (file.exists(file_path)) {
    data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # Ensure GeneSymbol is a column (it should already be from header=TRUE)
    # If GeneSymbol is in rownames, convert it to a column
    if (!"GeneSymbol" %in% colnames(data) && !is.null(rownames(data))) {
      data$GeneSymbol <- rownames(data)
      rownames(data) <- NULL
    }
    data_list[[dir]] <- data
  }
}

# Check if any data was loaded
if (length(data_list) == 0) {
  stop("ERROR: No data files found. Check that files exist with pattern: CIBERSORTxHiRes_NA_", celltype_filename, "_", window, ".txt")
}

# Bind all data frames by rows
merged_data <- do.call(rbind, data_list)

# Check if merged_data has rows and columns before using complete.cases
if (nrow(merged_data) > 0 && ncol(merged_data) > 0) {
  merged_data <- merged_data[complete.cases(merged_data), ]
  # Check if we still have data after filtering
  if (nrow(merged_data) == 0) {
    stop("ERROR: All rows were removed after filtering complete cases. Check for missing values in input files.")
  }
} else {
  stop("ERROR: Merged data has no rows or columns. Check input files.")
}

# Group by GeneSymbol and calculate the mean for duplicated genes
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
