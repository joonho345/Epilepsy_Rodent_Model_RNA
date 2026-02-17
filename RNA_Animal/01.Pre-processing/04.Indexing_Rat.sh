#!/bin/bash
#$ -cwd
#$ -S /bin/bash

ConfigFile="/home/joonho345/1_Epilepsy_RNA/script/Scratch_settings.sh"
source ${ConfigFile}

#ConfigFile=$1
#source ${ConfigFile}

#####
OutPath_35=${Rat_IndexPath_35_A}
OutPath_50=${Rat_IndexPath_50_A}
OutPath_75=${Rat_IndexPath_75_A}
OutPath_100=${Rat_IndexPath_100_A}
OutPath_125=${Rat_IndexPath_125_A}
OutPath_150=${Rat_IndexPath_150_A}

#####
make_dir ${OutPath_35}
make_dir ${OutPath_50}
make_dir ${OutPath_75}
make_dir ${OutPath_100}
make_dir ${OutPath_125}
make_dir ${OutPath_150}

#####
STAR \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir ${OutPath_35} \
--genomeFastaFiles ${REF_mRatBN7_2} \
--sjdbGTFfile ${GTF_mRatBN7_2_112} \
--sjdbOverhang 35
echo "This is Indexing of ${OutPath_35}"

STAR \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir ${OutPath_50} \
--genomeFastaFiles ${REF_mRatBN7_2} \
--sjdbGTFfile ${GTF_mRatBN7_2_112} \
--sjdbOverhang 50
echo "This is Indexing of ${OutPath_50}"

STAR \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir ${OutPath_75} \
--genomeFastaFiles ${REF_mRatBN7_2} \
--sjdbGTFfile ${GTF_mRatBN7_2_112} \
--sjdbOverhang 75
echo "This is Indexing of ${OutPath_75}"

STAR \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir ${OutPath_100} \
--genomeFastaFiles ${REF_mRatBN7_2} \
--sjdbGTFfile ${GTF_mRatBN7_2_112} \
--sjdbOverhang 100
echo "This is Indexing of ${OutPath_100}"

STAR \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir ${OutPath_125} \
--genomeFastaFiles ${REF_mRatBN7_2} \
--sjdbGTFfile ${GTF_mRatBN7_2_112} \
--sjdbOverhang 125
echo "This is Indexing of ${OutPath_125}"

STAR \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir ${OutPath_150} \
--genomeFastaFiles ${REF_mRatBN7_2} \
--sjdbGTFfile ${GTF_mRatBN7_2_112} \
--sjdbOverhang 150
echo "This is Indexing of ${OutPath_150}"

