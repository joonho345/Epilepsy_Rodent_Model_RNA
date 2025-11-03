#!/bin/bash
#$ -cwd
#$ -S /bin/bash

ConfigFile="/home/joonho345/3_RNA/script/Scratch_settings.sh"
source ${ConfigFile}

#ConfigFile=$1
#source ${ConfigFile}

#####
OutPath_50=${IndexPath_50}
OutPath_100=${IndexPath_100}
OutPath_125=${IndexPath_125}
OutPath_150=${IndexPath_150}

#####
make_dir ${OutPath_50}
make_dir ${OutPath_100}
make_dir ${OutPath_125}
make_dir ${OutPath_150}

#####
STAR \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir ${OutPath_125} \
--genomeFastaFiles ${REF_GRCh38_p14} \
--sjdbGTFfile ${GTF_GRCh38_p14_112} \
--sjdbOverhang 125
echo "This is Indexing of ${OutPath_125}"
