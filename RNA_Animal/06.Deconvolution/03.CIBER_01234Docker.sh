##################################
# pull docker image into singularity
singularity pull docker://cibersortx/fractions
singularity pull docker://cibersortx/hires
singularity pull docker://cibersortx/gep

GSE="GSE185862"

mixture="CIBERSORTx_adjusted_merged_matrix_TPM_M.txt"

# Define array of GOIs (Genes of Interest)
GOIs=("M_gene_set_ERK" "M_gene_set_IC" "M_gene_set_MF" "M_gene_set_ND" "M_gene_set_NG" "M_gene_set_NI" "M_gene_set_NT" "M_gene_set_SP" "M_gene_set_GG" "M_gene_set_TNF")

target="${GSE}_10X_scRNA_matrix_DEGs_CellType1"


function make_dir {
        if [ ! -d ${1}/ ];then
                mkdir -p ${1}/
        fi
}

########## cibersortx/fractions ##########
#### Run CIBERSORTx for signature matrix ####
cd /home/joonho345/resources/CIBERSORTx
scRNA_input=/home/joonho345/1_Epilepsy_RNA/scRNA_Animal/02.${GSE}_10X_Deconv/
sig_output=/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_${GSE}_01Sig/${target}
make_dir ${sig_output}

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
sig_output=/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_${GSE}_01Sig/${target}
fraction_output=/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_${GSE}_02Fraction/${target}
make_dir ${fraction_output}
## move mixture file to sig_output ##
mixture_source=/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/01.IMPORT_MIXTURE/${mixture}
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

#--rmbatchSmode TRUE \
#--refsample /src/data/${target}.txt \


########## cibersortx/hires ##########
# Run CIBERSORTx for HiRes - loop through each GOI
for GOI in "${GOIs[@]}"; do
  echo "  Processing GOI: ${GOI}"
  cd /home/joonho345/resources/CIBERSORTx
  sig_output=/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_${GSE}_01Sig/${target}
  hires_output=/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_${GSE}_04HiRes/${target}/${GOI}
  make_dir ${hires_output}
  ## move mixture file to sig_output ##
  ## move GOI file to sig_output ##
  mixture_source=/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/01.IMPORT_MIXTURE/${mixture}
  cp ${mixture_source} ${sig_output}/
  goi_source=/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/02.IMPORT_GO/${GOI}.txt
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


#####################################################
########## cibersortx/GEP ##########
# Run CIBERSORTx for Group Level GEPs
cd /home/joonho345/resources/CIBERSORTx
sig_output=/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_${GSE}_01Sig/${target}
gep_output=/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_${GSE}_03GEP/${target}
make_dir ${gep_output}
## move mixture file to sig_output ##
mixture_source=/home/joonho345/1_Epilepsy_RNA/RNA_Animal/06.Deconvolution/01.IMPORT_MIXTURE/${mixture}
cp ${mixture_source} ${sig_output}/

singularity exec \
  --bind ${sig_output}:/src/data \
  --bind ${gep_output}:/src/outdir \
  /home/joonho345/resources/CIBERSORTx/gep_latest.sif \
  /src/CIBERSORTxGep \
  --username joonho345@yuhs.ac \
  --token 3feb8102722435d25be86816493ebe18 \
  --sigmatrix /src/data/CIBERSORTx_${target}_inferred_phenoclasses.CIBERSORTx_${target}_inferred_refsample.bm.K999.txt \
  --mixture /src/data/${mixture} \
  --outdir /src/outdir
  #######################################################