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
