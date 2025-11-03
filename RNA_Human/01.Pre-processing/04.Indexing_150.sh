#!/bin/bash
#$ -cwd
#$ -S /bin/bash

ConfigFile="/home/joonho345/3_RNA/script/Scratch_settings.sh"
source ${ConfigFile}

#ConfigFile=$1
#source ${ConfigFile}

#####
OutPath_150=${IndexPath_150_MS}

#####
make_dir ${OutPath_150}

#####
STAR \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir ${OutPath_150} \
--genomeFastaFiles ${REF_GRCm39} \
--sjdbGTFfile ${GTF_GRCm39_112} \
--sjdbOverhang 150
echo "This is Indexing of ${OutPath_150}"

