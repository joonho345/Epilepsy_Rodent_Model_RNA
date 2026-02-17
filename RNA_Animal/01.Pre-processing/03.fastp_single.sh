#!/bin/bash
#$ -cwd
#$ -S /bin/bash

ID=$1
ConfigFile=$2
source ${ConfigFile}

#####
DataPath=${FastqPath_unzip_A}
OutPath=${FastpPath_A}

#####
make_dir ${OutPath}

#####
fastp \
	--in1 ${DataPath}${FASTQ0_unzip_A} \
	--out1 ${OutPath}${FASTP0_A} \
	--html ${OutPath}${ID}.report.html \
	--json ${OutPath}${ID}.report.json \
	--length_required 25 \
	--low_complexity_filter \
	--overrepresentation_analysis -P 20
