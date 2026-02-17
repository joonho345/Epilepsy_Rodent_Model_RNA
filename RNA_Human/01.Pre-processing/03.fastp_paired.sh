#!/bin/bash
#$ -cwd
#$ -S /bin/bash

ID=$1
ConfigFile=$2
source ${ConfigFile}

#####
DataPath=${FastqPath}
OutPath=${FastpPath}

#####
make_dir ${OutPath}

#####
fastp \
	--in1 ${DataPath}${FASTQ1} \
	--in2 ${DataPath}${FASTQ2} \
	--out1 ${OutPath}${FASTP1} \
	--out2 ${OutPath}${FASTP2} \
	--html ${OutPath}${ID}.report.html \
	--json ${OutPath}${ID}.report.json \
	--length_required 50 \
	--low_complexity_filter \
	--detect_adapter_for_pe \
	--overrepresentation_analysis -P 20

