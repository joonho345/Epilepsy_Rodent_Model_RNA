##################################
function make_dir {
        if [ ! -d ${1}/ ];then
                mkdir -p ${1}/
        fi
}

GSE="GSE185862"

mixture="CIBERSORTx_adjusted_merged_matrix_TPM_M.txt"

GOI="M_gene_set_CC"
GOI="M_gene_set_ERK"
GOI="M_gene_set_IC"
GOI="M_gene_set_MF"
GOI="M_gene_set_ND"
GOI="M_gene_set_NG"
GOI="M_gene_set_NI"
GOI="M_gene_set_NR"
GOI="M_gene_set_NT"
GOI="M_gene_set_SP"
GOI="M_gene_set_TLR"
GOI="M_gene_set_TNF"

target="GSE185862_scRNA_matrix_all_cluster_label"
target="GSE185862_scRNA_matrix_all_cluster_1"
target="GSE185862_scRNA_matrix_all_cluster_2"

target="GSE185862_scRNA_matrix_top2000_cluster_label"
target="GSE185862_scRNA_matrix_top2000_cluster_1"
target="GSE185862_scRNA_matrix_top2000_cluster_2"

target="GSE185862_scRNA_matrix_top1000_cluster_label"
target="GSE185862_scRNA_matrix_top1000_cluster_1"
target="GSE185862_scRNA_matrix_top1000_cluster_2"

# pull docker image into singularity 
singularity pull docker://cibersortx/fractions
singularity pull docker://cibersortx/gep
singularity pull docker://cibersortx/hires


########## cibersortx/fractions ##########
#### Run CIBERSORTx for signature matrix ####
cd /home/joonho345/CIBERSORTx
scRNA_input=/home/joonho345/3_RNA/scRNA_Animal/02.${GSE}_Seurat_Deconv/
sig_output=/home/joonho345/3_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_${GSE}_01Sig/${target}
make_dir ${sig_output}

singularity exec \
  --bind ${scRNA_input}:/src/data \
  --bind ${sig_output}:/src/outdir \
  /home/joonho345/CIBERSORTx/fractions_latest.sif \
  /src/CIBERSORTxFractions \
  --username joonho345@yuhs.ac \
  --token ff1b2ed2472db4a219523baa9103ab02 \
  --single_cell TRUE \
  --refsample /src/data/${target}.txt \
  --outdir /src/outdir

#### Run CIBERSORTx for fraction calculation ####
cd /home/joonho345/CIBERSORTx
sig_output=/home/joonho345/3_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_${GSE}_01Sig/${target}
fraction_output=/home/joonho345/3_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_${GSE}_02Fraction/${target}
make_dir ${fraction_output}
## move mixture file to sig_output ##

singularity exec \
  --bind ${sig_output}:/src/data \
  --bind ${fraction_output}:/src/outdir \
  /home/joonho345/CIBERSORTx/fractions_latest.sif \
  /src/CIBERSORTxFractions \
  --username joonho345@yuhs.ac \
  --token ff1b2ed2472db4a219523baa9103ab02 \
  --sigmatrix /src/data/CIBERSORTx_${target}_inferred_phenoclasses.CIBERSORTx_${target}_inferred_refsample.bm.K999.txt \
  --mixture /src/data/${mixture} \
  --outdir /src/outdir 

#--rmbatchSmode TRUE \
#--refsample /src/data/${target}.txt \


########## cibersortx/hires ##########
# Run CIBERSORTx for HiRes
cd /home/joonho345/CIBERSORTx
sig_output=/home/joonho345/3_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_${GSE}_01Sig/${target}
hires_output=/home/joonho345/3_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_${GSE}_04HiRes/${target}/${GOI}
make_dir ${hires_output}
## move mixture file to sig_output ##
## move GOI file to sig_output ##

singularity exec \
  --bind ${sig_output}:/src/data \
  --bind ${hires_output}:/src/outdir \
  /home/joonho345/CIBERSORTx/hires_latest.sif \
  /src/CIBERSORTxHiRes \
  --username joonho345@yuhs.ac \
  --token ff1b2ed2472db4a219523baa9103ab02 \
  --sigmatrix /src/data/CIBERSORTx_${target}_inferred_phenoclasses.CIBERSORTx_${target}_inferred_refsample.bm.K999.txt \
  --mixture /src/data/${mixture} \
  --subsetgenes /src/data/${GOI}.txt \
  --outdir /src/outdir


#####################################################
########## cibersortx/GEP ##########
# Run CIBERSORTx for Group Level GEPs
cd /home/joonho345/CIBERSORTx
sig_output=/home/joonho345/3_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_${GSE}_01Sig/${target}
gep_output=/home/joonho345/3_RNA/RNA_Animal/06.Deconvolution/03.CIBERSORT_${GSE}_03GEP/${target}
make_dir ${gep_output}
## move mixture file to sig_output ##

singularity exec \
  --bind ${sig_output}:/src/data \
  --bind ${gep_output}:/src/outdir \
  /home/joonho345/CIBERSORTx/gep_latest.sif \
  /src/CIBERSORTxGep \
  --username joonho345@yuhs.ac \
  --token ff1b2ed2472db4a219523baa9103ab02 \
  --sigmatrix /src/data/CIBERSORTx_${target}_inferred_phenoclasses.CIBERSORTx_${target}_inferred_refsample.bm.K999.txt \
  --mixture /src/data/${mixture} \
  --outdir /src/outdir
  #######################################################