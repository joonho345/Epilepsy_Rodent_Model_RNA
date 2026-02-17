# Load necessary libraries
library(dplyr)

Target_comparison <- 'FILTERED_1_MTLEALL_NL'

# Define the DESeq file and cutoff values
DESeq_file <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/07.DESeq/01.DESeq_files/DESeq_", Target_comparison,".txt")


##########################
cutoff_log2FoldChange <- 1
cutoff_padj <- 0.05

# Load data
data <- read.table(DESeq_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Filter genes based on the criteria using the cutoff variables
upregulated_genes <- rownames(data)[data$log2FoldChange > cutoff_log2FoldChange & data$padj < cutoff_padj]
downregulated_genes <- rownames(data)[data$log2FoldChange < -cutoff_log2FoldChange & data$padj < cutoff_padj]

# Define output file names using the variables
output_file_1 <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/08.Enrichment/GO_01.Matrix_DESeq/DESeq_", 
                        cutoff_log2FoldChange, "_", cutoff_padj, "_", Target_comparison, "_Up.txt")
output_file_2 <- paste0("/home/joonho345/1_Epilepsy_RNA/RNA_Human/08.Enrichment/GO_01.Matrix_DESeq/DESeq_", 
                        cutoff_log2FoldChange, "_", cutoff_padj, "_", Target_comparison,"_Down.txt")

# Export the gene lists to text files
write.table(upregulated_genes, file = output_file_1, sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(downregulated_genes, file = output_file_2, sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
