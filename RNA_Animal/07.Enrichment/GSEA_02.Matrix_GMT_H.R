gene_sets <- list(
  M_gene_set_GG_GO0014004_16 = M_gene_set_GG_GO0014004_16,
  M_gene_set_GG_GO0014014_59 = M_gene_set_GG_GO0014014_59,
  M_gene_set_GG_GO0014015_94 = M_gene_set_GG_GO0014015_94,
  M_gene_set_GG_GO0042063_410 = M_gene_set_GG_GO0042063_410,
  M_gene_set_GG_GO0048708_94 = M_gene_set_GG_GO0048708_94,
  M_gene_set_GG = M_gene_set_GG,
  M_gene_set_IC_GO0005245_43 = M_gene_set_IC_GO0005245_43,
  M_gene_set_IC_GO0005248_24 = M_gene_set_IC_GO0005248_24,
  M_gene_set_IC_GO0005249_98 = M_gene_set_IC_GO0005249_98,
  M_gene_set_IC = M_gene_set_IC,
  M_gene_set_MF_GO0048668_35 = M_gene_set_MF_GO0048668_35,
  M_gene_set_MF_GO0098686_61 = M_gene_set_MF_GO0098686_61,
  M_gene_set_MF = M_gene_set_MF,
  M_gene_set_ND_GO0043524_214 = M_gene_set_ND_GO0043524_214, 
  M_gene_set_ND_GO0043525_105 = M_gene_set_ND_GO0043525_105, 
  M_gene_set_ND_GO0051402_388 = M_gene_set_ND_GO0051402_388,
  M_gene_set_ND_GO0110088_8 = M_gene_set_ND_GO0110088_8,
  M_gene_set_ND = M_gene_set_ND,
  M_gene_set_NG_GO0002052_44 = M_gene_set_NG_GO0002052_44, 
  M_gene_set_NG_GO0007406_16 = M_gene_set_NG_GO0007406_16,
  M_gene_set_NG_GO0050768_163 = M_gene_set_NG_GO0050768_163, 
  M_gene_set_NG_GO0050769_331 = M_gene_set_NG_GO0050769_331,
  M_gene_set_NG_GO0050772_85 = M_gene_set_NG_GO0050772_85,
  M_gene_set_NG_GO0050771_56 = M_gene_set_NG_GO0050771_56,
  M_gene_set_NG = M_gene_set_NG,
  M_gene_set_NI_GO0001774_49 = M_gene_set_NI_GO0001774_49, 
  M_gene_set_NI_GO0006954_882 = M_gene_set_NI_GO0006954_882,
  M_gene_set_NI_GO0048143_20 = M_gene_set_NI_GO0048143_20,
  M_gene_set_NI_GO0050728_173 = M_gene_set_NI_GO0050728_173,
  M_gene_set_NI_GO0050729_153 = M_gene_set_NI_GO0050729_153,
  M_gene_set_NI_GO0150076_87 = M_gene_set_NI_GO0150076_87, 
  M_gene_set_NI = M_gene_set_NI,
  M_gene_set_NT_GO0008328_44 = M_gene_set_NT_GO0008328_44,
  M_gene_set_NT_GO0017146_9 = M_gene_set_NT_GO0017146_9, 
  M_gene_set_NT_GO0032281_29 = M_gene_set_NT_GO0032281_29,
  M_gene_set_NT_GO0032983_6 = M_gene_set_NT_GO0032983_6,
  M_gene_set_NT_GO0035249_131 = M_gene_set_NT_GO0035249_131,
  M_gene_set_NT_GO0051932_73 = M_gene_set_NT_GO0051932_73,
  M_gene_set_NT = M_gene_set_NT,
  M_gene_set_SP_GO0048168_78 = M_gene_set_SP_GO0048168_78,
  M_gene_set_SP_GO0048169_39 = M_gene_set_SP_GO0048169_39,
  M_gene_set_SP_GO0048172_22 = M_gene_set_SP_GO0048172_22,
  M_gene_set_SP = M_gene_set_SP,
  M_gene_set_TNF_GO0010804_27 = M_gene_set_TNF_GO0010804_27,
  M_gene_set_TNF_GO0033209_103 = M_gene_set_TNF_GO0033209_103,
  M_gene_set_TNF_GO1903265_11 = M_gene_set_TNF_GO1903265_11,
  M_gene_set_TNF = M_gene_set_TNF
)

# Define the output GMT file path
gmt_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/GSEA_02.Matrix_GSEA/output_gene_sets_M_new.gmt"

# Open a file connection for writing
file_conn <- file(gmt_file, "w")

# Process and write each gene set
for (name in names(gene_sets)) {
  # Extract the gene set name without the prefix
  gene_set_name <- sub("M_gene_set_", "", name)
  
  # Combine the gene set name, 'na', and the existing genes into a single line
  gmt_line <- paste(c(gene_set_name, "na", gene_sets[[name]]), collapse = "\t")
  # Write the line to the file
  writeLines(gmt_line, file_conn)
}

# Close the file connection
close(file_conn)



######################################################
gene_sets <- list(
  R_gene_set_GG_GO0014004_16 = R_gene_set_GG_GO0014004_16,
  R_gene_set_GG_GO0014014_59 = R_gene_set_GG_GO0014014_59,
  R_gene_set_GG_GO0014015_94 = R_gene_set_GG_GO0014015_94,
  R_gene_set_GG_GO0042063_410 = R_gene_set_GG_GO0042063_410,
  R_gene_set_GG_GO0048708_94 = R_gene_set_GG_GO0048708_94,
  R_gene_set_GG = R_gene_set_GG,
  R_gene_set_IC_GO0005245_43 = R_gene_set_IC_GO0005245_43,
  R_gene_set_IC_GO0005248_24 = R_gene_set_IC_GO0005248_24,
  R_gene_set_IC_GO0005249_98 = R_gene_set_IC_GO0005249_98,
  R_gene_set_IC = R_gene_set_IC,
  R_gene_set_MF_GO0048668_35 = R_gene_set_MF_GO0048668_35,
  R_gene_set_MF_GO0098686_61 = R_gene_set_MF_GO0098686_61,
  R_gene_set_MF = R_gene_set_MF,
  R_gene_set_ND_GO0043524_214 = R_gene_set_ND_GO0043524_214, 
  R_gene_set_ND_GO0043525_105 = R_gene_set_ND_GO0043525_105, 
  R_gene_set_ND_GO0051402_388 = R_gene_set_ND_GO0051402_388,
  R_gene_set_ND_GO0110088_8 = R_gene_set_ND_GO0110088_8,
  R_gene_set_ND = R_gene_set_ND,
  R_gene_set_NG_GO0002052_44 = R_gene_set_NG_GO0002052_44, 
  R_gene_set_NG_GO0007406_16 = R_gene_set_NG_GO0007406_16,
  R_gene_set_NG_GO0050768_163 = R_gene_set_NG_GO0050768_163, 
  R_gene_set_NG_GO0050769_331 = R_gene_set_NG_GO0050769_331,
  R_gene_set_NG_GO0050772_85 = R_gene_set_NG_GO0050772_85,
  R_gene_set_NG_GO0050771_56 = R_gene_set_NG_GO0050771_56,
  R_gene_set_NG = R_gene_set_NG,
  R_gene_set_NI_GO0001774_49 = R_gene_set_NI_GO0001774_49, 
  R_gene_set_NI_GO0006954_882 = R_gene_set_NI_GO0006954_882,
  R_gene_set_NI_GO0048143_20 = R_gene_set_NI_GO0048143_20,
  R_gene_set_NI_GO0050728_173 = R_gene_set_NI_GO0050728_173,
  R_gene_set_NI_GO0050729_153 = R_gene_set_NI_GO0050729_153,
  R_gene_set_NI_GO0150076_87 = R_gene_set_NI_GO0150076_87, 
  R_gene_set_NI = R_gene_set_NI,
  R_gene_set_NT_GO0008328_44 = R_gene_set_NT_GO0008328_44,
  R_gene_set_NT_GO0017146_9 = R_gene_set_NT_GO0017146_9, 
  R_gene_set_NT_GO0032281_29 = R_gene_set_NT_GO0032281_29,
  R_gene_set_NT_GO0032983_6 = R_gene_set_NT_GO0032983_6,
  R_gene_set_NT_GO0035249_131 = R_gene_set_NT_GO0035249_131,
  R_gene_set_NT_GO0051932_73 = R_gene_set_NT_GO0051932_73,
  R_gene_set_NT = R_gene_set_NT,
  R_gene_set_SP_GO0048168_78 = R_gene_set_SP_GO0048168_78,
  R_gene_set_SP_GO0048169_39 = R_gene_set_SP_GO0048169_39,
  R_gene_set_SP_GO0048172_22 = R_gene_set_SP_GO0048172_22,
  R_gene_set_SP = R_gene_set_SP,
  R_gene_set_TNF_GO0010804_27 = R_gene_set_TNF_GO0010804_27,
  R_gene_set_TNF_GO0033209_103 = R_gene_set_TNF_GO0033209_103,
  R_gene_set_TNF_GO1903265_11 = R_gene_set_TNF_GO1903265_11,
  R_gene_set_TNF = R_gene_set_TNF
)


# Define the output GMT file path
gmt_file <- "/home/joonho345/1_Epilepsy_RNA/RNA_Animal/08.Enrichment/GSEA_02.Matrix_GSEA/output_gene_sets_R_new.gmt"

# Open a file connection for writing
file_conn <- file(gmt_file, "w")

# Process and write each gene set
for (name in names(gene_sets)) {
  # Extract the gene set name without the prefix
  gene_set_name <- sub("R_gene_set_", "", name)
  
  # Combine the gene set name, 'na', and the existing genes into a single line
  gmt_line <- paste(c(gene_set_name, "na", gene_sets[[name]]), collapse = "\t")
  # Write the line to the file
  writeLines(gmt_line, file_conn)
}

# Close the file connection
close(file_conn)
