#!/bin/bash
#$ -cwd
#$ -S /bin/bash

ID=$1
ConfigFile=$2
source ${ConfigFile}

#####
DataPath=${FastqPath_unzip}
OutPath=${FastpPath}

#####
make_dir ${OutPath}

#####
fastp \
	--in1 ${DataPath}${FASTQ0_unzip} \
	--out1 ${OutPath}${FASTP0} \
	--html ${OutPath}${ID}.report.html \
	--json ${OutPath}${ID}.report.json \
	--length_required 25 \
	--low_complexity_filter \
	--overrepresentation_analysis -P 20

##########################
# https://github.com/OpenGene/fastp
# for single end data (not compressed)
# the output will be gzip-compressed if its file name ends with .gz / FASTP0_unzip=${ID}.FP_R0.fq??

#--trim_poly_x \
#--correction \
#p: --overrepresentation_analysis
#w: thread
#y:  --low_complexity_filter
#h: html result saved path
