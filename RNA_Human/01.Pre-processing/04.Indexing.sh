#!/bin/bash
#$ -cwd
#$ -S /bin/bash

ConfigFile="/home/joonho345/1_Epilepsy_RNA/script/Scratch_settings.sh"
source ${ConfigFile}

#ConfigFile=$1
#source ${ConfigFile}

#####
OutPath_50=${IndexPath_50}
OutPath_100=${IndexPath_100}
OutPath_125=${IndexPath_125}
OutPath_150=${IndexPath_150}

#####
STAR \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir ${OutPath_50} \
--genomeFastaFiles ${REF_GRCh38_p14} \
--sjdbGTFfile ${GTF_GRCh38_p14_112} \
--sjdbOverhang 50
echo "This is Indexing of ${OutPath_50}"

STAR \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir ${OutPath_100} \
--genomeFastaFiles ${REF_GRCh38_p14} \
--sjdbGTFfile ${GTF_GRCh38_p14_112} \
--sjdbOverhang 100
echo "This is Indexing of ${OutPath_100}"

STAR \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir ${OutPath_125} \
--genomeFastaFiles ${REF_GRCh38_p14} \
--sjdbGTFfile ${GTF_GRCh38_p14_112} \
--sjdbOverhang 125
echo "This is Indexing of ${OutPath_125}"

STAR \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir ${OutPath_150} \
--genomeFastaFiles ${REF_GRCh38_p14} \
--sjdbGTFfile ${GTF_GRCh38_p14_112} \
--sjdbOverhang 150
echo "This is Indexing of ${OutPath_150}"

