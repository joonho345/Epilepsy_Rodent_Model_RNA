library(dplyr)
library(ggplot2)

DataName <- 'DESeq_1_0.05_FILTERED_1_MTLEALL_NL_Up'
DataName <- 'DESeq_1_0.05_FILTERED_1_MTLEALL_NL_Down'


InputFile <- sprintf("/home/joonho345/1_Epilepsy_RNA/RNA_Human/08.Enrichment/GO_02.ToppFun/ToppFun_%s.txt", DataName)
GO_type <- "GO: Biological Process"
GO_name <- "Biological"

data <- read.table(InputFile, header = TRUE, sep = "\t", quote = "", comment.char = "")
filtered_data <- data %>% filter(Category == GO_type)

#### data filter
HitCountGenome <- 1000

# all
top_data <- filtered_data %>% 
  filter(Hit.Count.in.Genome <= HitCountGenome) %>%
  filter(q.value.FDR.B.H < 0.05) %>%
  arrange(q.value.FDR.B.H)
top_data <- top_data %>%
  mutate(log_q_value = -log10(`q.value.FDR.B.H`),
         gene_ratio = `Hit.Count.in.Query.List` / `Hit.Count.in.Genome`)


#### Export files
filename <- sprintf("/home/joonho345/1_Epilepsy_RNA/RNA_Human/08.Enrichment/GO_02.ToppFun/%s_%s_top_final.txt", GO_name, DataName)
write.table(top_data, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)


