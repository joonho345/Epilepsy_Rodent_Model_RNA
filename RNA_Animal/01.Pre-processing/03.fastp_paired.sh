#!/bin/bash
#$ -cwd
#$ -S /bin/bash

ID=$1
ConfigFile=$2
source ${ConfigFile}

#####
DataPath=${FastqPath_A}
OutPath=${FastpPath_A}

#####
make_dir ${OutPath}

#####
fastp \
	--in1 ${DataPath}${FASTQ1_A} \
	--in2 ${DataPath}${FASTQ2_A} \
	--out1 ${OutPath}${FASTP1_A} \
	--out2 ${OutPath}${FASTP2_A} \
	--html ${OutPath}${ID}.report.html \
	--json ${OutPath}${ID}.report.json \
	--length_required 50 \
	--low_complexity_filter \
	--detect_adapter_for_pe \
	--overrepresentation_analysis -P 20

# --length_required 25 \ for 35bp & 50bp reads

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
