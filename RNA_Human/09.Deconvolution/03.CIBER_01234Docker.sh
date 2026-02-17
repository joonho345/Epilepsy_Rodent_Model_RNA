##################################
# pull docker image into singularity
singularity pull docker://cibersortx/fractions
singularity pull docker://cibersortx/hires
singularity pull docker://cibersortx/gep

GSE="GSE160189_186538_Integration"
mixture="CIBERSORTx_adjusted_merged_matrix_1_TPM.txt"
GOIs=(
  "H_gene_set_ERK"
  "H_gene_set_IC"
  "H_gene_set_MF"
  "H_gene_set_ND"
  "H_gene_set_NG"
  "H_gene_set_NI"
  "H_gene_set_NT"
  "H_gene_set_SP"
  "H_gene_set_GG"
  "H_gene_set_TNF"
)
target="${GSE}_scRNA_matrix_DEGs_CellType1"

# Use scRNA data from scRNA_Human/02.${GSE}_Seurat_Deconv

########## cibersortx/fractions ##########
#### Run CIBERSORTx for signature matrix ####
cd /home/joonho345/resources/CIBERSORTx
scRNA_input=/home/joonho345/1_Epilepsy_RNA/scRNA_Human/02.${GSE}_Seurat_Deconv/
sig_output=/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/03.CIBERSORT_${GSE}_01Sig/${target}

singularity exec \
  --bind ${scRNA_input}:/src/data \
  --bind ${sig_output}:/src/outdir \
  /home/joonho345/resources/CIBERSORTx/fractions_latest.sif \
  /src/CIBERSORTxFractions \
  --username joonho345@yuhs.ac \
  --token 3feb8102722435d25be86816493ebe18 \
  --single_cell TRUE \
  --refsample /src/data/${target}.txt \
  --outdir /src/outdir

#### Run CIBERSORTx for fraction calculation ####
cd /home/joonho345/resources/CIBERSORTx
sig_output=/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/03.CIBERSORT_${GSE}_01Sig/${target}
fraction_output=/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/03.CIBERSORT_${GSE}_02Fraction/${target}
## move mixture file to sig_output ##
mixture_source=/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/01.IMPORT_MIXTURE/${mixture}
cp ${mixture_source} ${sig_output}/

singularity exec \
  --bind ${sig_output}:/src/data \
  --bind ${fraction_output}:/src/outdir \
  /home/joonho345/resources/CIBERSORTx/fractions_latest.sif \
  /src/CIBERSORTxFractions \
  --username joonho345@yuhs.ac \
  --token 3feb8102722435d25be86816493ebe18 \
  --sigmatrix /src/data/CIBERSORTx_${target}_inferred_phenoclasses.CIBERSORTx_${target}_inferred_refsample.bm.K999.txt \
  --mixture /src/data/${mixture} \
  --outdir /src/outdir 


########## cibersortx/hires ##########
# Run CIBERSORTx for HiRes - loop through each GOI
for GOI in "${GOIs[@]}"; do
  echo "  Processing GOI: ${GOI}"
  cd /home/joonho345/resources/CIBERSORTx
  sig_output=/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/03.CIBERSORT_${GSE}_01Sig/${target}
  hires_output=/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/03.CIBERSORT_${GSE}_04HiRes/${target}/${GOI}
  ## move mixture file to sig_output ##
  ## move GOI file to sig_output ##
  mixture_source=/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/01.IMPORT_MIXTURE/${mixture}
  cp ${mixture_source} ${sig_output}/
  goi_source=/home/joonho345/1_Epilepsy_RNA/RNA_Human/09.Deconvolution/02.IMPORT_GO/${GOI}.txt
  cp ${goi_source} ${sig_output}/

  singularity exec \
    --bind ${sig_output}:/src/data \
    --bind ${hires_output}:/src/outdir \
    /home/joonho345/resources/CIBERSORTx/hires_latest.sif \
    /src/CIBERSORTxHiRes \
    --username joonho345@yuhs.ac \
    --token 3feb8102722435d25be86816493ebe18 \
    --sigmatrix /src/data/CIBERSORTx_${target}_inferred_phenoclasses.CIBERSORTx_${target}_inferred_refsample.bm.K999.txt \
    --mixture /src/data/${mixture} \
    --subsetgenes /src/data/${GOI}.txt \
    --outdir /src/outdir
  
  echo "  Completed HiRes for ${GOI}"
done
echo "All HiRes analyses completed!"

