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

# GTF_GRCh38_p14_112 is not working
# Error occured when processing GFF file (line 8 of file /home/joonho345/resources/Reference/Homo_sapiens.GRCh38.112.gtf):
# Feature ENSG00000228037 does not contain a 'gene_name' attribute