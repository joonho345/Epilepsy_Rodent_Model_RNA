library(dplyr)


#### coldata ####
coldata_all <- read.csv("/home/joonho345/3_RNA/RNA_Human/Raw_Data/Patients_sheet_FTP.csv", header = TRUE)
coldata <- coldata_all %>% filter(Included == 1) # filtered 677 Samples
coldata <- coldata %>% filter(Brain_Location != "Blood") # filtered 661 Samples
# coldata <- coldata %>% filter(PRJNA != "PRJNA556159") # filtered 625 Samples
rownames(coldata) <- coldata$Run
coldata_df <- as.data.frame(coldata, stringsAsFactors = FALSE)


#### matrix ####
# Define the path to your files and list all .htseq.count.txt files
path <- "/data/project/HS/RNA_Human/02.Quantification/01.HTseq"
files <- list.files(path, pattern = "*.htseq.count.txt", full.names = TRUE)
# Read each file and store it in a list of data frames
dfs <- lapply(files, function(file) {
  df <- read.table(file, sep = "\t", header = FALSE, col.names = c("Gene", gsub(".htseq.count.txt", "", basename(file))))
  return(df)
})
# Merge all data frames by "Gene"
merged_df <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), dfs)
merged_df <- merged_df[-(1:5), ]
# Remove genes whose names start with 'ENSG'
merged_df <- merged_df[!grepl("^ENSG", merged_df$Gene), ]
rownames(merged_df) <- merged_df$Gene
merged_df <- merged_df[, -1]


#### Order of samples ####
# coldata & expression matrix -> 661 Samples 
sorted_colnames <- sort(colnames(merged_df))
sorted_rownames <- sort(rownames(coldata_df))
merged_df <- merged_df[, sorted_colnames]
coldata_df <- coldata_df[sorted_rownames, ]
filtered_samples <- rownames(coldata_df)
merged_df <- merged_df[, colnames(merged_df) %in% filtered_samples]


#### Gene_name: Save the matrices and coldata ####
merged_matrix <- as.matrix(merged_df)
write.table(merged_matrix, file = "/home/joonho345/3_RNA/RNA_Human/02.Quantification/merged_matrix.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
coldata_matrix <- as.matrix(coldata_df)
write.table(coldata_matrix, file = "/home/joonho345/3_RNA/RNA_Human/02.Quantification/filtered_coldata.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


#### Change the gene name into gene id ####
merged_df_id <- merged_df
merged_df_id$gene_name <- rownames(merged_df_id)
merged_df_id <- merge(merged_df_id, gtf_table[, c("gene_name", "gene_id")], by = "gene_name", all.x = TRUE)
if(any(is.na(merged_df_id$gene_id))) {
  warning("There are NA values in gene_id. These rows will be removed.")
  merged_df_id <- merged_df_id[!is.na(merged_df_id$gene_id), ]
}
merged_df_id <- merged_df_id %>% dplyr::select(-gene_name)
rownames(merged_df_id) <- merged_df_id$gene_id
merged_df_id <- merged_df_id %>% dplyr::select(-gene_id)


#### Gene_id: Save the matrices and coldata ####
merged_matrix_id <- as.matrix(merged_df_id)
write.table(merged_matrix_id, file = "/home/joonho345/3_RNA/RNA_Human/02.Quantification/merged_matrix_id.txt", 
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

