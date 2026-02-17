#!/bin/bash
#$ -cwd
#$ -S /bin/bash

ID=$1
ConfigFile=$2
source ${ConfigFile}

#####
DataPath=${AlignAllPath_A}
BAMFile="${ID}.bam"
OutPath=${HTseqPath_A}

#####
make_dir ${OutPath}

#### 01.HTseq
htseq-count \
-f bam \
-r name \
-s no \
-a 10 \
-t exon \
-m intersection-strict \
${DataPath}/${BAMFile} \
${GTF_GRCm39_112} > ${OutPath}/${ID}.htseq.count.txt

# -i gene_name \