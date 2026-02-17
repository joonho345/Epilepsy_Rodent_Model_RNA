#.Project Name
PROJECT_R="1_Epilepsy_RNA"
DIR_H="/home/joonho345/${PROJECT_R}/"
#### RNA_Human ####
DIR_H_H="/home/joonho345/${PROJECT_R}/RNA_Human/"
DIR_H="/data/project/${PROJECT_R}/RNA_Human/"
script_H=${DIR_H_H}script/
logPath_H=${DIR_H_H}out/
temp_H=${DIR_H}temp/
#### RNA_Animal ####
DIR_H_A="/home/joonho345/${PROJECT_R}/RNA_Animal/"
DIR_A="/data/project/${PROJECT_R}/RNA_Animal/"
script_A=${DIR_H_A}script/
logPath_A=${DIR_H_A}out/
temp_A=${DIR_A}temp/
#### scRNA_Human ####
DIR_H_S="/home/joonho345/${PROJECT_R}/scRNA_Human/"
DIR_S="/data/project/${PROJECT_R}/scRNA_Human/"
script_S=${DIR_H_S}script/
logPath_S=${DIR_H_S}out/
temp_S=${DIR_S}temp/
#### scRNA_Animal ####
DIR_H_C="/home/joonho345/${PROJECT_R}/scRNA_Animal/"
DIR_C="/data/project/${PROJECT_R}/scRNA_Animal/"
script_C=${DIR_H_C}script/
logPath_C=${DIR_H_C}out/
temp_C=${DIR_C}temp/


#5.Reference Path
GTF_GRCh38_p14_112="/home/joonho345/resources/Reference/Homo_sapiens.GRCh38.112.gtf"
GTF_GRCh38_p14_v46="/home/joonho345/resources/Reference/gencode.v46.annotation.gtf"
REF_GRCh38_p14="/home/joonho345/resources/Reference/GRCh38.p14.genome.fa"
GTF_GRCh38_p13_104="/home/joonho345/resources/Reference/Homo_sapiens.GRCh38.104.gtf"
GTF_GRCh38_p13_v38="/home/joonho345/resources/Reference/gencode.v38.annotation.gtf"
REF_GRCh38_p13="/home/joonho345/resources/Reference/GRCh38.p13.genome.fa"
GTF_mRatBN7_2_112="/home/joonho345/resources/Reference/Rattus_norvegicus.mRatBN7.2.112.gtf"
GFF_mRatBN7_2_112="/home/joonho345/resources/Reference/Rattus_norvegicus.mRatBN7.2.112.gff3"
REF_mRatBN7_2="/home/joonho345/resources/Reference/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa"
GTF_GRCm39_112="/home/joonho345/resources/Reference/Mus_musculus.GRCm39.112.gtf"
GFF_GRCm39_112="/home/joonho345/resources/Reference/Mus_musculus.GRCm39.112.gff3"
REF_GRCm39="/home/joonho345/resources/Reference/Mus_musculus.GRCm39.dna.primary_assembly.fa"


##Function
function make_dir {
        if [ ! -d ${1}/ ];then
                mkdir -p ${1}/
        fi
        #chown -R $(whoami):fdn ${1} >&/dev/null
}

##################################################################
###################### RNA_Human ##########################
#####Sample_List#####
SRR_ALL_LIST_PATH=${DIR_H_H}Raw_Data/SRR_Acc_List_ALL.txt
SRR_PAIR_LIST_PATH=${DIR_H_H}Raw_Data/SRR_Acc_List_Paired.txt
SRR_SINGLE_LIST_PATH=${DIR_H_H}Raw_Data/SRR_Acc_List_Single.txt
mapfile -t SRR_ALL_LIST < ${SRR_ALL_LIST_PATH}
mapfile -t SRR_PAIR_LIST < ${SRR_PAIR_LIST_PATH}
mapfile -t SRR_SINGLE_LIST < ${SRR_SINGLE_LIST_PATH}
SRR_GSE71058_LIST_PATH=${DIR_H_H}Raw_Data/Human_GSE71058/SRR_Acc_List_final.txt
mapfile -t SRR_GSE71058_LIST < ${SRR_GSE71058_LIST_PATH}
#####1.Pre-processing#####
#00.received
OutPath_SRA=${DIR_H}01.Pre-processing/00.recieved/
#01.raw
OutPath_SRA_raw=${DIR_H}01.Pre-processing/01.raw/
OutPath_SRA_raw_unzip=${DIR_H}01.Pre-processing/01.raw_single
FastqPath=${DIR_H}01.Pre-processing/01.raw/
FastqPath_unzip=${DIR_H}01.Pre-processing/01.raw_single/
FASTQ0_pre_LIST=$(ls ${FastqPath}| egrep "*_0.fq.gz")
FASTQ1_pre_LIST=$(ls ${FastqPath}| egrep "*_1.fq.gz")
FASTQ2_pre_LIST=$(ls ${FastqPath}| egrep "*_2.fq.gz")
FASTQ0_LIST=(${FASTQ0_pre_LIST// / })
FASTQ1_LIST=(${FASTQ1_pre_LIST// / })
FASTQ2_LIST=(${FASTQ2_pre_LIST// / })
FASTQ0="${ID}_0.fq.gz"
FASTQ0_unzip="${ID}_0.fq"
FASTQ1="${ID}_1.fq.gz"
FASTQ2="${ID}_2.fq.gz"
#03.fastp
FastpPath=${DIR_H}01.Pre-processing/03.fastp/${ID}/
FASTP0=${ID}.FP_R0.fq.gz
FASTP0_unzip=${ID}.FP_R0.fq
FASTP1=${ID}.FP_R1.fq.gz
FASTP2=${ID}.FP_R2.fq.gz
#04.Aligned_Index
IndexPath_50=${DIR_H}01.Pre-processing/04.aligned_index/01.Len_50/
SRR_50_LIST_PATH=${DIR_H_H}Raw_Data/SRR_Acc_List_Single.txt
mapfile -t SRR_50_LIST < ${SRR_50_LIST_PATH}
IndexPath_100=${DIR_H}01.Pre-processing/04.aligned_index/01.Len_100/
SRR_100_LIST_PATH=${DIR_H_H}Raw_Data/SRR_Acc_List_100.txt
mapfile -t SRR_100_LIST < ${SRR_100_LIST_PATH}
IndexPath_125=${DIR_H}01.Pre-processing/04.aligned_index/01.Len_125/
SRR_125_LIST_PATH=${DIR_H_H}Raw_Data/SRR_Acc_List_125.txt
mapfile -t SRR_125_LIST < ${SRR_125_LIST_PATH}
IndexPath_150=${DIR_H}01.Pre-processing/04.aligned_index/01.Len_150/
SRR_150_LIST_PATH=${DIR_H_H}Raw_Data/SRR_Acc_List_150.txt
mapfile -t SRR_150_LIST < ${SRR_150_LIST_PATH}
#04.Aligned
AlignPath=${DIR_H}01.Pre-processing/04.aligned/${ID}/
AlignAllPath=${DIR_H}01.Pre-processing/04.aligned/ALL/

#####2.Quantification#####
BAMAll_pre_LIST=$(ls ${AlignAllPath}| egrep "*.bam"| grep -v ".bai")
BAMAll_LIST=(${BAMAll_pre_LIST// / })
#01.HTseq
HTseqPath=${DIR_H}02.Quantification/01.HTseq/


##################################################################
###################### RNA_Animal ##########################
#####1.Pre-processing#####
##00.received
OutPath_SRA_A=${DIR_A}01.Pre-processing/00.recieved/
##01.raw
OutPath_SRA_raw_A=${DIR_A}01.Pre-processing/01.raw/
OutPath_SRA_raw_run_A=${DIR_A}01.Pre-processing/01.raw_run/
OutPath_SRA_raw_unzip_A=${DIR_A}01.Pre-processing/01.raw_single
FastqPath_A=${DIR_A}01.Pre-processing/01.raw_run/
FastqPath_unzip_A=${DIR_A}01.Pre-processing/01.raw_single/
FASTQ0_pre_LIST_A=$(ls ${FastqPath_A}| egrep "*_0.fq.gz")
FASTQ1_pre_LIST_A=$(ls ${FastqPath_A}| egrep "*_1.fq.gz")
FASTQ2_pre_LIST_A=$(ls ${FastqPath_A}| egrep "*_2.fq.gz")
FASTQ0_LIST_A=(${FASTQ0_pre_LIST_A// / })
FASTQ1_LIST_A=(${FASTQ1_pre_LIST_A// / })
FASTQ2_LIST_A=(${FASTQ2_pre_LIST_A// / })
FASTQ0_A="${ID}_0.fq.gz"
FASTQ0_unzip_A="${ID}_0.fq"
FASTQ1_A="${ID}_1.fq.gz"
FASTQ2_A="${ID}_2.fq.gz"
##02.fastqc
FastqcPath_A=${DIR_A}01.Pre-processing/02.fastqc/
##03.fastp
FastpPath_A=${DIR_A}01.Pre-processing/03.fastp/${ID}/
FASTP0_A=${ID}.FP_R0.fq.gz
FASTP0_unzip_A=${ID}.FP_R0.fq
FASTP1_A=${ID}.FP_R1.fq.gz
FASTP2_A=${ID}.FP_R2.fq.gz
SRR_35_LIST_SINGLE_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_Single_35.txt
mapfile -t SRR_35_LIST_SINGLE_A < ${SRR_35_LIST_SINGLE_PATH_A}
SRR_50_LIST_SINGLE_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_Single_50.txt
mapfile -t SRR_50_LIST_SINGLE_A < ${SRR_50_LIST_SINGLE_PATH_A}
SRR_75_LIST_SINGLE_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_Single_75.txt
mapfile -t SRR_75_LIST_SINGLE_A < ${SRR_75_LIST_SINGLE_PATH_A}
SRR_35_LIST_PAIRED_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_Paired_35.txt
mapfile -t SRR_35_LIST_PAIRED_A < ${SRR_35_LIST_PAIRED_PATH_A}
SRR_50_LIST_PAIRED_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_Paired_50.txt
mapfile -t SRR_50_LIST_PAIRED_A < ${SRR_50_LIST_PAIRED_PATH_A}
SRR_75_LIST_PAIRED_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_Paired_75.txt
mapfile -t SRR_75_LIST_PAIRED_A < ${SRR_75_LIST_PAIRED_PATH_A}
SRR_100_LIST_PAIRED_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_Paired_100.txt
mapfile -t SRR_100_LIST_PAIRED_A < ${SRR_100_LIST_PAIRED_PATH_A}
SRR_125_LIST_PAIRED_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_Paired_125.txt
mapfile -t SRR_125_LIST_PAIRED_A < ${SRR_125_LIST_PAIRED_PATH_A}
SRR_150_LIST_PAIRED_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_Paired_150.txt
mapfile -t SRR_150_LIST_PAIRED_A < ${SRR_150_LIST_PAIRED_PATH_A}
##04.Aligned_Index
#Mouse
Mouse_IndexPath_35_A=${DIR_A}01.Pre-processing/04.aligned_index/01.Mouse_Len_35/
Mouse_IndexPath_50_A=${DIR_A}01.Pre-processing/04.aligned_index/01.Mouse_Len_50/
Mouse_IndexPath_75_A=${DIR_A}01.Pre-processing/04.aligned_index/01.Mouse_Len_75/
Mouse_IndexPath_100_A=${DIR_A}01.Pre-processing/04.aligned_index/01.Mouse_Len_100/
Mouse_IndexPath_125_A=${DIR_A}01.Pre-processing/04.aligned_index/01.Mouse_Len_125/
Mouse_IndexPath_150_A=${DIR_A}01.Pre-processing/04.aligned_index/01.Mouse_Len_150/
#Rat
Rat_IndexPath_35_A=${DIR_A}01.Pre-processing/04.aligned_index/01.Rat_Len_35/
Rat_IndexPath_50_A=${DIR_A}01.Pre-processing/04.aligned_index/01.Rat_Len_50/
Rat_IndexPath_75_A=${DIR_A}01.Pre-processing/04.aligned_index/01.Rat_Len_75/
Rat_IndexPath_100_A=${DIR_A}01.Pre-processing/04.aligned_index/01.Rat_Len_100/
Rat_IndexPath_125_A=${DIR_A}01.Pre-processing/04.aligned_index/01.Rat_Len_125/
Rat_IndexPath_150_A=${DIR_A}01.Pre-processing/04.aligned_index/01.Rat_Len_150/
##04.Aligned
AlignPath_A=${DIR_A}01.Pre-processing/04.aligned/${ID}/
AlignAllPath_A=${DIR_A}01.Pre-processing/04.aligned/ALL/
#Mouse
SRR_M_PAIRED_35_LIST_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_M_Paired_35.txt
mapfile -t SRR_M_PAIRED_35_LIST_A < ${SRR_M_PAIRED_35_LIST_PATH_A}
SRR_M_PAIRED_50_LIST_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_M_Paired_50.txt
mapfile -t SRR_M_PAIRED_50_LIST_A < ${SRR_M_PAIRED_50_LIST_PATH_A}
SRR_M_PAIRED_100_LIST_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_M_Paired_100.txt
mapfile -t SRR_M_PAIRED_100_LIST_A < ${SRR_M_PAIRED_100_LIST_PATH_A}
SRR_M_PAIRED_125_LIST_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_M_Paired_125.txt
mapfile -t SRR_M_PAIRED_125_LIST_A < ${SRR_M_PAIRED_125_LIST_PATH_A}
SRR_M_PAIRED_150_LIST_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_M_Paired_150.txt
mapfile -t SRR_M_PAIRED_150_LIST_A < ${SRR_M_PAIRED_150_LIST_PATH_A}
SRR_M_SINGLE_50_LIST_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_M_Single_50.txt
mapfile -t SRR_M_SINGLE_50_LIST_A < ${SRR_M_SINGLE_50_LIST_PATH_A}
SRR_M_SINGLE_75_LIST_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_M_Single_75.txt
mapfile -t SRR_M_SINGLE_75_LIST_A < ${SRR_M_SINGLE_75_LIST_PATH_A}
#Rat
SRR_R_PAIRED_75_LIST_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_R_Paired_75.txt
mapfile -t SRR_R_PAIRED_75_LIST_A < ${SRR_R_PAIRED_75_LIST_PATH_A}
SRR_R_PAIRED_125_LIST_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_R_Paired_125.txt
mapfile -t SRR_R_PAIRED_125_LIST_A < ${SRR_R_PAIRED_125_LIST_PATH_A}
SRR_R_PAIRED_150_LIST_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_R_Paired_150.txt
mapfile -t SRR_R_PAIRED_150_LIST_A < ${SRR_R_PAIRED_150_LIST_PATH_A}
SRR_R_SINGLE_35_LIST_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_R_Single_35.txt
mapfile -t SRR_R_SINGLE_35_LIST_A < ${SRR_R_SINGLE_35_LIST_PATH_A}

#####2.Quantification#####
BAMAll_pre_LIST_A=$(ls ${AlignAllPath_A}| egrep "*.bam"| grep -v ".bai")
BAMAll_LIST_A=(${BAMAll_pre_LIST_A// / })
SRR_M_LIST_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_M.txt
mapfile -t SRR_M_LIST_A < ${SRR_M_LIST_PATH_A}
SRR_R_LIST_PATH_A=${DIR_H_A}Raw_Data/SRR_Acc_List_R.txt
mapfile -t SRR_R_LIST_A < ${SRR_R_LIST_PATH_A}
#01.HTseq
HTseqPath_A=${DIR_A}02.Quantification/01.HTseq/
