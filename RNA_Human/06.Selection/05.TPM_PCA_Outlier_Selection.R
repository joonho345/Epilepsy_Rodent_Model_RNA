#############################
# TPM-PCA Common Outlier Selection Script
# Finds outliers common to both TPM-based and PCA-based methods
#############################

library(dplyr)

###############################################################
##################### Define Helper Functions ####################
###############################################################

# Function to get cohort info for a list of samples
get_cohort_info <- function(sample_list, coldata) {
  if (length(sample_list) == 0) {
    return(data.frame())
  }
  
  # Handle concatenated sample names (e.g., "SRR16522693SRR16522703")
  # Split and check each part
  all_samples <- character()
  for (sample in sample_list) {
    # Try to find exact match first
    if (sample %in% rownames(coldata)) {
      all_samples <- c(all_samples, sample)
    } else {
      # Try to split and find matches (for concatenated names)
      # Look for patterns like SRR followed by numbers
      parts <- strsplit(sample, "(?<=[0-9])(?=SRR)", perl = TRUE)[[1]]
      for (part in parts) {
        if (part %in% rownames(coldata)) {
          all_samples <- c(all_samples, part)
        }
      }
    }
  }
  
  if (length(all_samples) == 0) {
    return(data.frame())
  }
  
  # Get unique samples
  all_samples <- unique(all_samples)
  
  # Extract cohort information (include PRJNA if available)
  cohort_cols <- c("Diagnosis", "Diagnosis_Sub", "Brain_Location", "Brain_Location_Sub")
  if ("PRJNA" %in% colnames(coldata)) {
    cohort_cols <- c(cohort_cols, "PRJNA")
  }
  
  cohort_info <- coldata[all_samples, cohort_cols, drop = FALSE]
  cohort_info$Sample <- rownames(cohort_info)
  rownames(cohort_info) <- NULL
  
  return(cohort_info)
}

###############################################################
##################### Read TPM Outliers ####################
###############################################################

tpm_outlier_dir <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/06.Selection/01.TPM_Outliers"

# Read TPM outliers
tpm_all_outliers_file <- paste0(tpm_outlier_dir, "/01.ALL_Outliers_All_Genes.txt")
tpm_mtle_outliers_file <- paste0(tpm_outlier_dir, "/01.MTLE_Outliers_All_Genes.txt")
tpm_nl_outliers_file <- paste0(tpm_outlier_dir, "/01.NL_Outliers_All_Genes.txt")

tpm_all_outliers <- if (file.exists(tpm_all_outliers_file) && file.info(tpm_all_outliers_file)$size > 0) {
  tryCatch({
    tmp <- trimws(read.table(tpm_all_outliers_file, stringsAsFactors = FALSE)$V1)
    tmp[tmp != ""]
  }, error = function(e) {
    character(0)
  })
} else {
  character(0)
}

tpm_mtle_outliers <- if (file.exists(tpm_mtle_outliers_file) && file.info(tpm_mtle_outliers_file)$size > 0) {
  tryCatch({
    tmp <- trimws(read.table(tpm_mtle_outliers_file, stringsAsFactors = FALSE)$V1)
    tmp[tmp != ""]
  }, error = function(e) {
    character(0)
  })
} else {
  character(0)
}

tpm_nl_outliers <- if (file.exists(tpm_nl_outliers_file) && file.info(tpm_nl_outliers_file)$size > 0) {
  tryCatch({
    tmp <- trimws(read.table(tpm_nl_outliers_file, stringsAsFactors = FALSE)$V1)
    tmp[tmp != ""]
  }, error = function(e) {
    character(0)
  })
} else {
  character(0)
}

###############################################################
##################### Read PCA Outliers ####################
###############################################################

pca_outlier_dir <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/06.Selection/03.PCA_Outliers"

# Read PCA outliers
pca_all_outliers_file <- paste0(pca_outlier_dir, "/03.ALL_Outliers_PCA.txt")
pca_mtle_outliers_file <- paste0(pca_outlier_dir, "/03.MTLE_Outliers_PCA.txt")
pca_nl_outliers_file <- paste0(pca_outlier_dir, "/03.NL_Outliers_PCA.txt")

pca_all_outliers <- if (file.exists(pca_all_outliers_file) && file.info(pca_all_outliers_file)$size > 0) {
  tryCatch({
    tmp <- trimws(read.table(pca_all_outliers_file, stringsAsFactors = FALSE)$V1)
    tmp[tmp != ""]
  }, error = function(e) {
    character(0)
  })
} else {
  character(0)
}

pca_mtle_outliers <- if (file.exists(pca_mtle_outliers_file) && file.info(pca_mtle_outliers_file)$size > 0) {
  tryCatch({
    tmp <- trimws(read.table(pca_mtle_outliers_file, stringsAsFactors = FALSE)$V1)
    tmp[tmp != ""]
  }, error = function(e) {
    character(0)
  })
} else {
  character(0)
}

pca_nl_outliers <- if (file.exists(pca_nl_outliers_file) && file.info(pca_nl_outliers_file)$size > 0) {
  tryCatch({
    tmp <- trimws(read.table(pca_nl_outliers_file, stringsAsFactors = FALSE)$V1)
    tmp[tmp != ""]
  }, error = function(e) {
    character(0)
  })
} else {
  character(0)
}

###############################################################
##################### Find Common Outliers ####################
###############################################################

# Common outliers (intersection of TPM and PCA)
common_all_outliers <- intersect(tpm_all_outliers, pca_all_outliers)
common_mtle_outliers <- intersect(tpm_mtle_outliers, pca_mtle_outliers)
common_nl_outliers <- intersect(tpm_nl_outliers, pca_nl_outliers)

###############################################################
##################### Save Results ####################
###############################################################

# Create output directory
output_dir <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/06.Selection/05.TPM_PCA_Outliers"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save common outliers
if (length(common_all_outliers) > 0) {
  file_name <- paste0(output_dir, "/05.ALL_Outliers_TPM_PCA_Common.txt")
  write.table(common_all_outliers, file = file_name, 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

if (length(common_mtle_outliers) > 0) {
  file_name <- paste0(output_dir, "/05.MTLE_Outliers_TPM_PCA_Common.txt")
  write.table(common_mtle_outliers, file = file_name, 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

if (length(common_nl_outliers) > 0) {
  file_name <- paste0(output_dir, "/05.NL_Outliers_TPM_PCA_Common.txt")
  write.table(common_nl_outliers, file = file_name, 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

###############################################################
##################### Get Cohort Information ####################
###############################################################

# Read coldata to get cohort information
coldata_df <- read.table("/home/joonho345/1_Epilepsy_RNA/RNA_Human/02.Quantification/filtered_coldata.txt",
                         sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE, 
                         fill = TRUE, quote = "", comment.char = "")

cat("\n=== Cohort Information for Common Outliers ===\n")
if (length(common_all_outliers) > 0) {
  cohort_info <- get_cohort_info(common_all_outliers, coldata_df)
  if (nrow(cohort_info) > 0) {
    cat("\nCohort breakdown:\n")
    for (i in 1:nrow(cohort_info)) {
      cat("  ", cohort_info$Sample[i], ": Diagnosis =", cohort_info$Diagnosis[i], 
          ", Diagnosis_Sub =", cohort_info$Diagnosis_Sub[i],
          ", Brain_Location =", cohort_info$Brain_Location[i],
          ", Brain_Location_Sub =", cohort_info$Brain_Location_Sub[i])
      if ("PRJNA" %in% colnames(cohort_info)) {
        cat(", PRJNA =", cohort_info$PRJNA[i])
      }
      cat("\n")
    }
    
    # Summary by Diagnosis_Sub
    cat("\nSummary by Diagnosis_Sub:\n")
    if ("Diagnosis_Sub" %in% colnames(cohort_info)) {
      diag_sub_table <- table(cohort_info$Diagnosis_Sub)
      for (diag in names(diag_sub_table)) {
        cat("  ", diag, ":", diag_sub_table[diag], "samples\n")
      }
    }
    
    # Summary by Brain_Location_Sub
    cat("\nSummary by Brain_Location_Sub:\n")
    if ("Brain_Location_Sub" %in% colnames(cohort_info)) {
      location_sub_table <- table(cohort_info$Brain_Location_Sub)
      for (loc in names(location_sub_table)) {
        cat("  ", loc, ":", location_sub_table[loc], "samples\n")
      }
    }
    
    # Summary by PRJNA
    if ("PRJNA" %in% colnames(cohort_info)) {
      cat("\nSummary by PRJNA:\n")
      prjna_table <- table(cohort_info$PRJNA)
      for (prjna in names(prjna_table)) {
        cat("  ", prjna, ":", prjna_table[prjna], "samples\n")
      }
    }
  } else {
    cat("  ⚠ Could not find cohort information for common outliers\n")
  }
} else {
  cat("  No common outliers to summarize\n")
}
sink()

