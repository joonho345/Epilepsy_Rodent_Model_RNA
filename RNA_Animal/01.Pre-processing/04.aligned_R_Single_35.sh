#!/bin/bash
#$ -cwd
#$ -S /bin/bash

ID=$1
ConfigFile=$2

source ${ConfigFile}

#####
IndexPath=${Rat_IndexPath_35_A}
DataPath=${FastpPath_A}
FASTQ0=${FASTP0_A} # this script is for FASTQ filtered by FASTP
OutPath=${AlignPath_A}
AllOutPath=${AlignAllPath_A}
Read_len=35

#####
make_dir ${OutPath}
make_dir ${AllOutPath}

#####
STAR \
--runThreadN 10 \
--runMode alignReads \
--genomeDir ${IndexPath} \
--readFilesIn ${DataPath}${FASTQ0} \
--sjdbOverhang ${Read_len} \
--sjdbGTFfile ${GTF_mRatBN7_2_112} \
--readFilesCommand gunzip -c \
--genomeLoad NoSharedMemory \
--twopassMode Basic \
--outFileNamePrefix ${OutPath}/${ID}_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes All

# --quantMode TranscriptomeSAM \
# --outSAMunmapped Within \

##### Indexing
samtools index ${OutPath}${ID}_Aligned.sortedByCoord.out.bam ${OutPath}${ID}_Aligned.sortedByCoord.out.bam.bai

##### Link
ln -s ${OutPath}${ID}_Aligned.sortedByCoord.out.bam ${AllOutPath}${ID}.bam
ln -s ${OutPath}${ID}_Aligned.sortedByCoord.out.bam.bai ${AllOutPath}${ID}.bam.bai
