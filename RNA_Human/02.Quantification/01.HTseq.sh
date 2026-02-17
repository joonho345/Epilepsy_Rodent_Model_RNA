#!/bin/bash
#$ -cwd
#$ -S /bin/bash

ID=$1
ConfigFile=$2
source ${ConfigFile}

#####
DataPath=${AlignAllPath}
BAMFile="${ID}.bam"
OutPath=${HTseqPath}

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
-i gene_name \
${DataPath}/${BAMFile} \
${GTF_GRCh38_p14_v46} > ${OutPath}/${ID}.htseq.count.txt
