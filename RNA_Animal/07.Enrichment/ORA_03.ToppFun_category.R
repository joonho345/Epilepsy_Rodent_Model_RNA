### to check each terms which category they fit well

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

# Define gene sets (same for all data names)
gene_sets <- list(
  H_gene_set_IC  = H_gene_set_IC,
  H_gene_set_NT  = H_gene_set_NT,
  H_gene_set_ND  = H_gene_set_ND,
  H_gene_set_NI  = H_gene_set_NI,
  H_gene_set_TNF = H_gene_set_TNF,
  H_gene_set_MF  = H_gene_set_MF,
  H_gene_set_NG  = H_gene_set_NG,
  H_gene_set_GG  = H_gene_set_GG,
  H_gene_set_SP  = H_gene_set_SP,
  H_gene_set_ERK = H_gene_set_ERK
)

# Loop through each data name
for (DataName in data_names) {
  cat("\n", "=", rep("=", 60), "\n", sep = "")
  cat("Processing:", DataName, "\n")
  cat("=", rep("=", 60), "\n", sep = "")
  
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
    arrange(p.value)
  
  # Check if top_data is empty after filtering
  if (nrow(top_data) == 0) {
    cat("Warning: No significant results found for", DataName, "after filtering. Skipping...\n")
    next
  }
  
  top_data <- top_data %>%
    mutate(log_q_value = -log10(`q.value.FDR.B.H`),
           gene_ratio = `Hit.Count.in.Query.List` / `Hit.Count.in.Genome`)
  
  # Add Color column if GO_term_colors exists, otherwise skip
  if (exists("GO_term_colors")) {
    top_data <- top_data %>%
      mutate(Color = GO_term_colors[ID])
  }
  
  # Get top 40 rows
  top_data <- top_data %>% head(40)
  
  # Check if we have enough rows
  if (nrow(top_data) == 0) {
    cat("Warning: No rows available for", DataName, ". Skipping...\n")
    next
  }
  
  # Extract gene lists from each row
  rows <- list()
  for (i in 1:min(40, nrow(top_data))) {
    rows[[i]] <- top_data[i, 'Hit.in.Query.List']
  }
  
  # Create a list to store the results for all rows
  all_membership_results <- list()
  
  # Loop through each row and calculate membership
  for (i in seq_along(rows)) {
    # Skip if row is NA or empty
    if (is.na(rows[[i]]) || length(rows[[i]]) == 0 || rows[[i]] == "") {
      next
    }
    
    # Input genes as a comma-separated string for the current row
    genes_input <- rows[[i]]
    genes_vector <- strsplit(genes_input, ",\\s*")[[1]]
    
    # Remove empty strings
    genes_vector <- genes_vector[genes_vector != ""]
    
    if (length(genes_vector) == 0) {
      next
    }
    
    # Check which gene sets each gene belongs to
    gene_set_membership <- lapply(genes_vector, function(gene) {
      sets <- names(gene_sets)[sapply(gene_sets, function(set) gene %in% set)]
      if (length(sets) == 0) {
        return("None")
      } else {
        return(sets)
      }
    })
    
    # Flatten the list of memberships to a single vector
    all_gene_sets <- unlist(gene_set_membership)
    
    # Remove "None" from counts if present
    all_gene_sets <- all_gene_sets[all_gene_sets != "None"]
    
    # Count the occurrences of each gene set
    if (length(all_gene_sets) > 0) {
      gene_set_counts <- table(all_gene_sets)
      
      # Find the most frequently occurring gene set(s)
      max_count <- max(gene_set_counts)
      most_common_gene_sets <- names(gene_set_counts[gene_set_counts == max_count])
      
      # Get GO term name and ID for this row
      go_term_name <- ifelse("Name" %in% colnames(top_data), as.character(top_data[i, "Name"]), "")
      go_term_id <- ifelse("ID" %in% colnames(top_data), as.character(top_data[i, "ID"]), "")
      
      # Store the result for the current row
      all_membership_results[[paste0("row_", i)]] <- list(
        Membership = data.frame(
          Gene = genes_vector,
          GeneSets = sapply(gene_set_membership, function(sets) paste(sets, collapse = ", "))
        ),
        MostCommonGeneSets = most_common_gene_sets,
        GeneSetCounts = gene_set_counts,
        GOTermName = go_term_name,
        GOTermID = go_term_id
      )
    }
  }
  
  # Prepare output data for writing
  output_lines <- c()
  output_lines <- c(output_lines, paste0("Processing: ", DataName))
  output_lines <- c(output_lines, "")
  
  # Calculate summary for top 10 rows
  if (length(all_membership_results) > 0) {
    output_lines <- c(output_lines, "Summary: Top 10 GO Terms and Their Most Common Gene Sets")
    output_lines <- c(output_lines, "")
    
    # Collect gene sets for counting
    top10_gene_sets <- c()
    
    for (i in 1:min(10, length(all_membership_results))) {
      row_key <- paste0("row_", i)
      if (row_key %in% names(all_membership_results)) {
        result <- all_membership_results[[row_key]]
        
        # Get GO term info
        go_term_id <- result$GOTermID
        go_term_name <- result$GOTermName
        most_common <- result$MostCommonGeneSets
        
        # Collect gene sets for summary count
        if (length(most_common) > 0) {
          top10_gene_sets <- c(top10_gene_sets, most_common)
        }
        
        # Format the line
        if (nchar(go_term_name) > 0) {
          if (length(most_common) > 0) {
            gene_sets_str <- paste(most_common, collapse = ", ")
            output_lines <- c(output_lines, paste0("  ", i, ". ", go_term_name, " (", go_term_id, ") -> ", gene_sets_str))
          } else {
            output_lines <- c(output_lines, paste0("  ", i, ". ", go_term_name, " (", go_term_id, ") -> None"))
          }
        } else if (nchar(go_term_id) > 0) {
          if (length(most_common) > 0) {
            gene_sets_str <- paste(most_common, collapse = ", ")
            output_lines <- c(output_lines, paste0("  ", i, ". ", go_term_id, " -> ", gene_sets_str))
          } else {
            output_lines <- c(output_lines, paste0("  ", i, ". ", go_term_id, " -> None"))
          }
        } else {
          if (length(most_common) > 0) {
            gene_sets_str <- paste(most_common, collapse = ", ")
            output_lines <- c(output_lines, paste0("  ", i, ". Row ", i, " -> ", gene_sets_str))
          } else {
            output_lines <- c(output_lines, paste0("  ", i, ". Row ", i, " -> None"))
          }
        }
      }
    }
    
    output_lines <- c(output_lines, "")
    
    # Add summary count of GO terms per category
    if (length(top10_gene_sets) > 0) {
      top10_summary <- table(top10_gene_sets)
      top10_summary <- sort(top10_summary, decreasing = TRUE)
      
      output_lines <- c(output_lines, "Summary: Number of GO Terms per Category (Top 10)")
      output_lines <- c(output_lines, "")
      for (j in 1:length(top10_summary)) {
        output_lines <- c(output_lines, paste0("  ", names(top10_summary)[j], ": ", top10_summary[j]))
      }
      output_lines <- c(output_lines, "")
    }
    
    output_lines <- c(output_lines, paste(rep("=", 60), collapse = ""))
    output_lines <- c(output_lines, "")
  }
  
  # Print and collect the most common gene sets for each row
  if (length(all_membership_results) > 0) {
    for (i in seq_along(all_membership_results)) {
      row_key <- paste0("row_", i)
      if (row_key %in% names(all_membership_results)) {
        result <- all_membership_results[[row_key]]
        
        # Print to console
        cat("\n", DataName, " - Row_", i, ":\n", sep = "")
        cat("Most Common Gene Sets:\n")
        print(result$MostCommonGeneSets)
        cat("\nGene Set Counts:\n")
        print(result$GeneSetCounts)
        
        # Add to output lines
        output_lines <- c(output_lines, paste0(DataName, " - Row_", i, ":"))
        if (nchar(result$GOTermID) > 0) {
          output_lines <- c(output_lines, paste0("GO Term ID: ", result$GOTermID))
        }
        if (nchar(result$GOTermName) > 0) {
          output_lines <- c(output_lines, paste0("GO Term Name: ", result$GOTermName))
        }
        output_lines <- c(output_lines, "Most Common Gene Sets:")
        if (length(result$MostCommonGeneSets) > 0) {
          output_lines <- c(output_lines, paste(result$MostCommonGeneSets, collapse = ", "))
        } else {
          output_lines <- c(output_lines, "None")
        }
        output_lines <- c(output_lines, "")
        output_lines <- c(output_lines, "Gene Set Counts:")
        if (length(result$GeneSetCounts) > 0) {
          counts_df <- data.frame(
            GeneSet = names(result$GeneSetCounts),
            Count = as.numeric(result$GeneSetCounts)
          )
          for (j in 1:nrow(counts_df)) {
            output_lines <- c(output_lines, paste0("  ", counts_df$GeneSet[j], ": ", counts_df$Count[j]))
          }
        } else {
          output_lines <- c(output_lines, "  None")
        }
        output_lines <- c(output_lines, "")
      }
    }
  } else {
    cat("No results found for", DataName, "\n")
    output_lines <- c(output_lines, paste0("No results found for ", DataName))
  }
  
  # Write results to file
  output_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/Gene_03.ToppFun/", prefix, "ToppFun_", DataName, "_Category.txt")
  writeLines(output_lines, output_file)
  cat("\nResults written to:", output_file, "\n")
}

