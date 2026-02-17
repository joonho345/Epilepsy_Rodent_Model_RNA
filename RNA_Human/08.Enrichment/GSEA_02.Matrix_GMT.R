gene_sets <- list(
  H_gene_set_ERK_GO0007265_122 = H_gene_set_ERK_GO0007265_122,
  H_gene_set_ERK_GO0051057_94 = H_gene_set_ERK_GO0051057_94,
  H_gene_set_ERK_GO0051058_64 = H_gene_set_ERK_GO0051058_64,
  H_gene_set_ERK_GO0070372_335 = H_gene_set_ERK_GO0070372_335,
  H_gene_set_ERK_GO0070373_83 = H_gene_set_ERK_GO0070373_83,
  H_gene_set_ERK_GO0070374_227 = H_gene_set_ERK_GO0070374_227,
  H_gene_set_ERK = H_gene_set_ERK,
  H_gene_set_GG_GO0014004_16 = H_gene_set_GG_GO0014004_16,
  H_gene_set_GG_GO0014014_59 = H_gene_set_GG_GO0014014_59,
  H_gene_set_GG_GO0014015_94 = H_gene_set_GG_GO0014015_94,
  H_gene_set_GG_GO0042063_410 = H_gene_set_GG_GO0042063_410,
  H_gene_set_GG_GO0048708_94 = H_gene_set_GG_GO0048708_94,
  H_gene_set_GG = H_gene_set_GG,
  H_gene_set_IC_GO0005245_43 = H_gene_set_IC_GO0005245_43,
  H_gene_set_IC_GO0005248_24 = H_gene_set_IC_GO0005248_24,
  H_gene_set_IC_GO0005249_98 = H_gene_set_IC_GO0005249_98,
  H_gene_set_IC = H_gene_set_IC,
  H_gene_set_MF_GO0048668_35 = H_gene_set_MF_GO0048668_35,
  H_gene_set_MF_GO0098686_61 = H_gene_set_MF_GO0098686_61,
  H_gene_set_MF = H_gene_set_MF,
  H_gene_set_ND_GO0043524_214 = H_gene_set_ND_GO0043524_214, 
  H_gene_set_ND_GO0043525_105 = H_gene_set_ND_GO0043525_105, 
  H_gene_set_ND_GO0051402_388 = H_gene_set_ND_GO0051402_388,
  H_gene_set_ND_GO0110088_8 = H_gene_set_ND_GO0110088_8,
  H_gene_set_ND = H_gene_set_ND,
  H_gene_set_NG_GO0002052_44 = H_gene_set_NG_GO0002052_44, 
  H_gene_set_NG_GO0007406_16 = H_gene_set_NG_GO0007406_16,
  H_gene_set_NG_GO0050768_163 = H_gene_set_NG_GO0050768_163, 
  H_gene_set_NG_GO0050769_331 = H_gene_set_NG_GO0050769_331,
  H_gene_set_NG_GO0050772_85 = H_gene_set_NG_GO0050772_85,
  H_gene_set_NG_GO0050771_56 = H_gene_set_NG_GO0050771_56,
  H_gene_set_NG = H_gene_set_NG,
  H_gene_set_NI_GO0001774_49 = H_gene_set_NI_GO0001774_49, 
  H_gene_set_NI_GO0006954_882 = H_gene_set_NI_GO0006954_882,
  H_gene_set_NI_GO0048143_20 = H_gene_set_NI_GO0048143_20,
  H_gene_set_NI_GO0050728_173 = H_gene_set_NI_GO0050728_173,
  H_gene_set_NI_GO0050729_153 = H_gene_set_NI_GO0050729_153,
  H_gene_set_NI_GO0150076_87 = H_gene_set_NI_GO0150076_87, 
  H_gene_set_NI = H_gene_set_NI,
  H_gene_set_NT_GO0008328_44 = H_gene_set_NT_GO0008328_44,
  H_gene_set_NT_GO0017146_9 = H_gene_set_NT_GO0017146_9, 
  H_gene_set_NT_GO0032281_29 = H_gene_set_NT_GO0032281_29,
  H_gene_set_NT_GO0032983_6 = H_gene_set_NT_GO0032983_6,
  H_gene_set_NT_GO0035249_131 = H_gene_set_NT_GO0035249_131,
  H_gene_set_NT_GO0051932_73 = H_gene_set_NT_GO0051932_73,
  H_gene_set_NT = H_gene_set_NT,
  H_gene_set_SP_GO0048168_78 = H_gene_set_SP_GO0048168_78,
  H_gene_set_SP_GO0048169_39 = H_gene_set_SP_GO0048169_39,
  H_gene_set_SP_GO0048172_22 = H_gene_set_SP_GO0048172_22,
  H_gene_set_SP = H_gene_set_SP,
  H_gene_set_TNF_GO0010804_27 = H_gene_set_TNF_GO0010804_27,
  H_gene_set_TNF_GO0033209_103 = H_gene_set_TNF_GO0033209_103,
  H_gene_set_TNF_GO1903265_11 = H_gene_set_TNF_GO1903265_11,
  H_gene_set_TNF = H_gene_set_TNF
)


# Define the output GMT file path
gmt_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Human/08.Enrichment/GSEA_02.Matrix_GSEA/output_gene_sets_new.gmt"

# Open a file connection for writing
file_conn <- file(gmt_file, "w")

# Process and write each gene set
for (name in names(gene_sets)) {
  # Extract the gene set name without the prefix
  gene_set_name <- sub("H_gene_set_", "", name)
  
  # Combine the gene set name, 'na', and the existing genes into a single line
  gmt_line <- paste(c(gene_set_name, "na", gene_sets[[name]]), collapse = "\t")
  print(gmt_line)
  # Write the line to the file
  writeLines(gmt_line, file_conn)
}

# Close the file connection
close(file_conn)
